import pytest
import os
import uuid

from scm.plams.core.jobmanager import JobManager
from scm.plams.core.settings import JobManagerSettings
from scm.plams.core.errors import PlamsError
from scm.plams.core.functions import use_subdir
from scm.plams.unit_tests.test_basejob import DummySingleJob


class TestJobManager:

    def test_lazy_workdir(self):
        # Given job manager
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)

        # When first initialised
        # Then workdir does not exist
        assert not os.path.exists(job_manager._workdir)

        # When access workdir for the first time
        # Then workdir is created
        workdir = job_manager.workdir
        assert os.path.exists(workdir)
        assert os.path.exists(job_manager._workdir)

        # When access subsequent time
        # Then same workdir is returned
        assert job_manager.workdir == workdir

        os.rmdir(job_manager.workdir)

    def test_load_and_clean_do_not_create_workdir(self):
        # Given job manager
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)

        # When load job
        job = DummySingleJob()
        job.run()
        job.results.wait()
        job_manager.load_job(f"{job.path}/{job.name}.dill")

        # Then workdir not created
        assert not os.path.exists(job_manager._workdir)

        # When clean the jobmanager
        job_manager._clean()

        # Then workdir not created
        assert not os.path.exists(job_manager._workdir)

    @pytest.mark.parametrize(
        "path_exists,folder_exists,use_existing_folder,expected_workdir",
        [
            (True, False, False, "./{}/{}"),
            (True, True, False, "./{}/{}.002"),
            (True, False, True, "./{}/{}"),
            (True, True, True, "./{}/{}"),
            (False, False, False, None),
        ],
        ids=[
            "path_exists_new_folder",
            "path_exists_folder_renamed",
            "path_exists_new_folder_with_use_existing",
            "path_exists_reuse_folder_with_use_existing",
            "path_not_exists_errors",
        ],
    )
    def test_workdir_location(self, path_exists, folder_exists, use_existing_folder, expected_workdir):
        # Given path and folder which may already exist
        path = str(uuid.uuid4())
        folder = str(uuid.uuid4())
        expected_workdir = expected_workdir.format(path, folder) if expected_workdir else None
        if path_exists:
            os.mkdir(path)
            if folder_exists:
                os.mkdir(f"{path}/{folder}")

        if expected_workdir is None:
            # When create jobmanager where path does not exist
            # Then raises error
            with pytest.raises(PlamsError):
                job_manager = JobManager(
                    settings=JobManagerSettings(), path=path, folder=folder, use_existing_folder=use_existing_folder
                )
        else:
            # When create jobmanager where path and folder may exist
            job_manager = JobManager(
                settings=JobManagerSettings(), path=path, folder=folder, use_existing_folder=use_existing_folder
            )

            # Then workdir is created
            assert os.path.abspath(expected_workdir) == job_manager.workdir
            assert os.path.exists(job_manager.workdir)

            job_manager._clean()
            if os.path.exists(job_manager.workdir):
                os.rmdir(job_manager.workdir)

        if os.path.exists(f"{path}/{folder}"):
            os.rmdir(f"{path}/{folder}")
        if os.path.exists(path):
            os.rmdir(path)

    def test_job_registration_load_deletion(self):
        # Given job manager
        folder = str(uuid.uuid4())
        job_manager = JobManager(settings=JobManagerSettings(), folder=folder)

        # When register jobs (implicitly when running)
        base_name = "test_jobreg"
        job1 = DummySingleJob(name=base_name)
        job2 = DummySingleJob(name=base_name)
        job3 = DummySingleJob(name=base_name)
        job4 = DummySingleJob(name=base_name)
        jobs = [job1, job2, job3, job4]
        job_manager._register(job1)
        job_manager._register(job2)
        with use_subdir("foo"):
            job_manager._register(job3)
            with use_subdir("bar"):
                job_manager._register(job4)

        # Then jobs registered as expected
        def verify_job_registration(job, jm, expected_name, expected_subdir=None):
            # Verify job manager set on job
            assert job.jobmanager == jm
            # Verify job name has postfix if duplicate run
            assert job.name == expected_name
            # Verify job path in workdir/subdir/job
            if expected_subdir:
                assert job.path == f"{job_manager.workdir}/{expected_subdir}/{expected_name}"
            else:
                assert job.path == f"{job_manager.workdir}/{expected_name}"
            # Verify job status is registered
            assert job.status == "registered"

        verify_job_registration(job1, job_manager, base_name)
        verify_job_registration(job2, job_manager, f"{base_name}.002")
        verify_job_registration(job3, job_manager, base_name, expected_subdir="foo")
        verify_job_registration(job4, job_manager, base_name, expected_subdir="foo/bar")
        assert job_manager.jobs == jobs
        assert job_manager.names == {"foo/bar/test_jobreg": 1, "foo/test_jobreg": 1, "test_jobreg": 2}

        # Given job manager and saved jobs
        folder2 = str(uuid.uuid4())
        job_manager2 = JobManager(settings=JobManagerSettings(), folder=folder2)
        for job in jobs:
            job.pickle()

        # When load jobs
        loaded_jobs = []
        for job in jobs:
            loaded_jobs.append(job_manager2.load_job(f"{job.path}/{job.name}.dill"))

        # Then jobs loaded correctly
        verify_job_registration(loaded_jobs[0], job_manager2, base_name)
        verify_job_registration(loaded_jobs[1], job_manager2, f"{base_name}.002")
        verify_job_registration(loaded_jobs[2], job_manager2, base_name, expected_subdir="foo")
        verify_job_registration(loaded_jobs[3], job_manager2, base_name, expected_subdir="foo/bar")

        # Given same jobs
        # When removed
        for job in jobs:
            job_manager.remove_job(job)

        # Then jobs removed from job manager
        assert job_manager.jobs == []

        job_manager._clean()
        if os.path.exists(job_manager.workdir):
            os.rmdir(job_manager.workdir)
