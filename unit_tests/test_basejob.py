import pytest
import uuid
from unittest.mock import patch
from datetime import datetime

from scm.plams.core.settings import Settings
from scm.plams.core.basejob import SingleJob
from scm.plams.core.errors import PlamsError, FileError
from scm.plams.core.jobrunner import JobRunner
from scm.plams.core.functions import add_to_instance
from scm.plams.core.enums import JobStatus


class DummySingleJob(SingleJob):
    """
    Single Job for testing.
    By default, creates an input file with a unique id of the form: 'Dummy input <id>'.
    When the job is run, creates an output file with the contents: 'Dummy output <id>'.
    There is also an option to add an arbitrary wait period before the job runs to simulate a longer-running calculation.
    """

    def __init__(self, inp=None, cmd=None, wait=0.0, **kwargs):
        super().__init__(**kwargs)
        self.calls = []
        self.id = uuid.uuid4()  # Ensure input unique
        self.input = f"Dummy input {self.id}" if inp is None else inp.replace("%ID%", str(self.id))
        self.command = "sed 's/input/output/g'" if cmd is None else cmd
        self.wait = wait

    def prerun(self):
        self.calls.append((self.prerun.__name__, datetime.now()))

    def postrun(self):
        self.calls.append((self.postrun.__name__, datetime.now()))

    def get_input(self):
        return self.input

    def get_runscript(self):
        self.calls.append((self.get_runscript.__name__, datetime.now()))
        return f"sleep {self.wait} && {self.command} {self._filename('inp')}"

    def check(self):
        x = self.results.read_file(self._filename("err"))
        return x == ""


class TestSingleJob:
    """
    Test suite for the Single Job.
    Not truly independent as relies upon the job runner/manager and results components.
    But this suite focuses on testing the methods on the job class itself.
    """

    def test_full_runscript_applies_pre_and_post_runscript_settings(self):
        # Given runscript settings
        s = Settings()
        s.runscript.shebang = "#!/bin/sh"
        s.runscript.pre = "# Add pre line"
        s.runscript.post = "\n# Add post line"

        # When get the runscript from the job
        job = DummySingleJob(settings=s)
        runscript = job.full_runscript()

        # Then the script is composed as expected
        assert (
            runscript
            == """#!/bin/sh

# Add pre line

sleep 0.0 && sed 's/input/output/g' plamsjob.in
# Add post line

""".replace(
                "\r\n", "\n"
            )
        )

    def test_run_wrapped_in_pre_and_postrun(self):
        # When run a job
        job = DummySingleJob()
        job.run()

        # Then makes calls to pre- and post-run
        assert [c for c, _ in job.calls] == ["prerun", "get_runscript", "postrun"]

    def test_run_writes_output_file_on_success(self):
        # When run a successful job
        job = DummySingleJob()
        results = job.run()

        # Then job check passes and output file written
        assert job.check()
        assert results.read_file("$JN.in") == f"Dummy input {job.id}"
        assert results.read_file("$JN.out") == f"Dummy output {job.id}"
        assert results.read_file("$JN.err") == ""

    def test_run_writes_error_file_on_failure(self):
        # When run an erroring job
        job = DummySingleJob(cmd="not_a_cmd")
        results = job.run()

        # Then job check fails and error file written
        assert not job.check()
        assert results.read_file("$JN.in") == f"Dummy input {job.id}"
        assert results.read_file("$JN.out") == ""
        assert results.read_file("$JN.err") != ""

    @pytest.mark.parametrize(
        "mode,expected",
        [
            [None, {(1, 1): None, (1, 2): None, (1, 3): None, (2, 2): None, (2, 3): None, (3, 3): None}],
            ["input", {(1, 1): True, (1, 2): False, (1, 3): True, (2, 2): True, (2, 3): False, (3, 3): True}],
            ["runscript", {(1, 1): True, (1, 2): True, (1, 3): False, (2, 2): True, (2, 3): False, (3, 3): True}],
            [
                "input+runscript",
                {(1, 1): True, (1, 2): False, (1, 3): False, (2, 2): True, (2, 3): False, (3, 3): True},
            ],
            ["not_a_mode", None],
        ],
        ids=["no_hashing", "input_hashing", "runscript_hashing", "both_hashing", "invalid_hashing"],
    )
    def test_hash_respects_mode(self, mode, expected, config):
        # Given jobs with different inputs and/or runscripts
        with patch("scm.plams.core.basejob.config", config):
            config.jobmanager.hashing = mode
            s = Settings()
            s.runscript.shebang = "#!/bin/sh"
            job1 = DummySingleJob(settings=s)
            job2 = DummySingleJob(inp=job1.input.replace("input", "inputx"), settings=s)
            job3 = DummySingleJob(inp=job1.input, cmd="echo 'foo' && sed 's/input/output/g'", settings=s)
            jobs = [job1, job2, job3]

            # When call hash with different modes
            # Then hashes match as expected
            if expected is None:
                with pytest.raises(PlamsError):
                    job1.hash()
            else:
                hashes = [
                    (i + 1, j + 1, job_i.hash(), job_j.hash())
                    for i, job_i in enumerate(jobs)
                    for j, job_j in enumerate(jobs)
                    if j >= i
                ]
                matches = {(i, j): None if h_i is None and h_j is None else h_i == h_j for i, j, h_i, h_j in hashes}
                assert matches == expected

    def test_run_multiple_independent_jobs_in_parallel(self):
        # Given parallel job runner
        runner = JobRunner(parallel=True, maxjobs=2)

        # When set up two jobs with no dependencies
        job1 = DummySingleJob(wait=0.2)
        job2 = DummySingleJob(wait=0.01)
        results = [job1.run(runner), job2.run(runner)]

        # Then both run in parallel
        # Shorter job finishes first even though started second
        # Pre-run call of second job is made before post-run call of first job
        results[1].wait()
        assert job2.status == JobStatus.SUCCESSFUL
        assert job1.status == JobStatus.RUNNING

        results[0].wait()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.calls[0][1] < job1.calls[2][1]

    def test_run_multiple_dependent_jobs_in_serial(self):
        # Given parallel job runner
        runner = JobRunner(parallel=True, maxjobs=2)

        # When set up two jobs with dependency
        job1 = DummySingleJob(wait=0.2)
        job2 = DummySingleJob(wait=0.01, depend=[job1])
        results = [job1.run(runner), job2.run(runner)]

        # Then run in serial
        # Second job finishes second even though shorter
        # Pre-run call of second job is made after post-run call of first job
        results[0].wait()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.status in [JobStatus.REGISTERED, JobStatus.STARTED, JobStatus.RUNNING]

        results[1].wait()
        assert job2.status == JobStatus.SUCCESSFUL
        assert job2.calls[0][1] >= job1.calls[2][1]

    def test_run_multiple_prerun_dependent_jobs_in_serial(self):
        # Given parallel job runner
        runner = JobRunner(parallel=True, maxjobs=2)

        # When set up two jobs with dependency via prerun
        job1 = DummySingleJob(wait=0.2)
        job2 = DummySingleJob(wait=0.01)

        @add_to_instance(job2)
        def prerun(s):
            s.calls.append((s.prerun.__name__, datetime.now()))
            job1.results.wait()

        results = [job1.run(runner), job2.run(runner)]

        # Then run in serial
        # Second job finishes second even though shorter
        # Get runscript call of second job is made after post-run call of first job
        results[0].wait()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.status in [JobStatus.REGISTERED, JobStatus.STARTED, JobStatus.RUNNING]

        results[1].wait()
        assert job2.status == JobStatus.SUCCESSFUL
        assert job2.calls[1][1] >= job1.calls[2][1]

    def test_ok_waits_on_results_and_checks_status(self):
        # Given job and a copy
        job1 = DummySingleJob(wait=0.1)
        job2 = DummySingleJob(inp=job1.input)
        job3 = DummySingleJob(cmd="not_a_cmd")

        # Then not ok before being started
        assert not job1.ok()

        # When call run and check ok
        job1.run()
        job2.run()
        job3.run()

        # Then waits for job to finish and checks status
        assert job1.ok()
        assert job1.status == JobStatus.SUCCESSFUL
        assert job2.ok()
        assert job2.status == JobStatus.COPIED
        assert not job3.ok()
        assert job3.status == JobStatus.CRASHED

    def test_load_non_existent_job_errors(self):
        # Given path that does not point to a job
        # When call load, then gives error
        with pytest.raises(FileError):
            DummySingleJob.load("not_a_path")

    def test_load_previously_run_job_succeeds(self):
        # Given successful previously run job
        job1 = DummySingleJob()
        job1.run()

        # When call load
        job2 = DummySingleJob.load(job1.path)

        # Then job loaded successfully
        assert job1.name == job2.name
        assert job1.id == job2.id
        assert job1.path == job2.path
        assert job1.settings == job2.settings
        assert job1._filenames == job2._filenames
