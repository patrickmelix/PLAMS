import pytest

from scm.plams.core.jobrunner import JobRunner
from scm.plams.unit_tests.test_basejob import DummySingleJob, log_call


class LoggedJobRunner(JobRunner):
    """
    Test job runner which logs call timings.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.call_log = {}

    def add_call_log_entry(self, entry):
        if entry.method in self.call_log:
            if isinstance(self.call_log[entry.method], list):
                self.call_log[entry.method].append(entry)
            else:
                self.call_log[entry.method] = [self.call_log[entry.method], entry]
        else:
            self.call_log[entry.method] = entry

    @log_call
    def call(self, runscript, workdir, out, err, runflags):
        return super().call(runscript, workdir, out, err, runflags)


class TestJobRunner:

    def test_serial_execution(self, config):
        # Given a job runner with parallel set to False in c'tor or via setter
        jobrunner1 = LoggedJobRunner()
        jobrunner2 = LoggedJobRunner(parallel=True)
        jobrunner2.parallel = False

        # When run jobs
        # Then verify run in serial
        self.verify_serial(jobrunner1, config)
        self.verify_serial(jobrunner2, config)

    @pytest.mark.parametrize(
        "maxjobs,maxthreads,use_setter",
        [
            (2, 256, False),
            (3, 256, True),
            (4, 256, False),
            (5, 256, True),
            (0, 2, False),
            (0, 3, True),
            (0, 4, False),
            (0, 5, True),
            (0, 10, False),
            (4, 2, False),
            (2, 4, False),
        ],
    )
    def test_parallel_execution_with_limit(self, maxjobs, maxthreads, use_setter, config):
        # Given a job runner with parallel constraint in c'tor or via setter
        if use_setter:
            jobrunner = LoggedJobRunner(parallel=True)
            jobrunner.maxjobs = maxjobs
            jobrunner.maxthreads = maxthreads
        else:
            jobrunner = LoggedJobRunner(parallel=True, maxjobs=maxjobs, maxthreads=maxthreads)

        # When run jobs
        # Then verify run in parallel with constraints
        self.verify_parallel_limit(jobrunner, config, maxjobs, maxthreads)

    def test_serial_to_parallel_toggle(self, config):
        # Given a serial job runner
        jobrunner = LoggedJobRunner(maxjobs=2)

        # When run jobs
        # Then runs in serial
        self.verify_serial(jobrunner, config)

        # When set to parallel and maxjobs increased
        jobrunner.parallel = True
        jobrunner.maxjobs = 0

        # When run jobs
        # Then runs in parallel
        self.verify_parallel_limit(jobrunner, config, 0, 256)

    def test_limit_setters(self, config):
        # Given jobrunner
        jobrunner = LoggedJobRunner(parallel=True)

        # When set maxjobs to <0 or 1
        # Then get error
        with pytest.raises(ValueError):
            jobrunner.maxjobs = -1
        with pytest.raises(ValueError):
            jobrunner.maxjobs = 1

        # When set maxjobs to 0
        # Then no semaphore
        jobrunner.maxjobs = 0
        assert jobrunner.maxjobs == 0
        assert jobrunner._job_limit is None

        # When set maxjobs to a new value
        # Then semaphore set
        jobrunner.maxjobs = 10
        assert jobrunner.maxjobs == 10
        assert jobrunner._job_limit is not None
        assert jobrunner._job_limit.max_value == 10

        # When set maxjobs to 0
        # Then no semaphore
        jobrunner.maxjobs = 0
        assert jobrunner.maxjobs == 0
        assert jobrunner._job_limit is None

        # When set maxthreads to <0 or 1
        # Then get error
        with pytest.raises(ValueError):
            jobrunner.maxthreads = -1
        with pytest.raises(ValueError):
            jobrunner.maxthreads = 1

        # When set maxthreads to 0
        # Then no semaphore
        jobrunner.maxthreads = 0
        assert jobrunner.maxthreads == 0
        assert jobrunner._jobthread_limit is None

        # When set maxthreads to a new value
        # Then semaphore set
        jobrunner.maxthreads = 10
        assert jobrunner.maxthreads == 10
        assert jobrunner._jobthread_limit is not None
        assert jobrunner._jobthread_limit.max_value == 10

        # When set maxthreads to 0
        # Then no semaphore
        jobrunner.maxthreads = 0
        assert jobrunner.maxthreads == 0
        assert jobrunner._job_limit is None

    def verify_serial(self, jobrunner, config):
        # When run 5 jobs
        jobs = []
        for i in range(5):
            job = DummySingleJob(wait=0.1)
            jobs.append(job)
            jobrunner._run_job(job, config.default_jobmanager)

        # Then check either that successive jobs are only executed after the previous has finished
        for i, j in enumerate(jobs):
            j.results.wait()
            if i > 0:
                # Then successive jobs are only executed after the previous has finished
                assert jobs[i].call_log["_execute"].timestamp >= jobs[i - 1].call_log["_finalize"].timestamp

    def verify_parallel_limit(self, jobrunner, config, maxjobs, maxthreads):
        jobs = []
        thread_limited = maxthreads < maxjobs or maxjobs == 0
        limit = maxthreads if thread_limited else maxjobs

        # When run 5 jobs in parallel with limit
        for i in range(5):
            job = DummySingleJob(wait=0.2)
            jobs.append(job)
            jobrunner._run_job(job, config.default_jobmanager)

        times = []
        for i, j in enumerate(jobs):
            j.results.wait()
            if thread_limited:
                start = jobs[i].call_log["_prepare"].timestamp
                end = jobs[i].call_log["_finalize"].timestamp
            else:
                start = jobrunner.call_log["call"][i].timestamp
                end = jobs[i].call_log["_finalize"].timestamp
            times.append((start, end))

        events = [e for s, e in times for e in [(s, 1), (e, -1)]]
        events.sort(key=lambda x: (x[0], x[1]))

        max_parallel_jobs = 0
        current_jobs = 0
        for _, change in events:
            current_jobs += change
            max_parallel_jobs = max(max_parallel_jobs, current_jobs)

        # Then number of jobs running in parallel is within limit
        assert max_parallel_jobs <= limit
