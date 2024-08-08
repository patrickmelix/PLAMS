import pytest
from unittest.mock import patch, MagicMock, mock_open
import time
import threading
import os

from scm.plams.core.functions import init, finish
from scm.plams.core.settings import Settings
from scm.plams.unit_tests.test_helpers import assert_config_as_expected


class TestInitAndFinish:
    """
    Test suite for the init() and finish() public function
    """

    @pytest.fixture(autouse=True)
    def mock_jobmanager(self):
        """
        Mock out the job manager to avoid creating run directories etc.
        """
        with patch("scm.plams.core.jobmanager.JobManager") as mock_jobmanager:
            yield mock_jobmanager

    @pytest.fixture(autouse=True)
    def mock_file(self):
        """
        Mock out file for logging.
        """
        with patch("builtins.open", new_callable=mock_open) as mock_file:
            yield mock_file

    def test_init_sets_config_settings_default_values(self, config):
        # When call init
        init()

        # Then config settings set to default values
        assert_config_as_expected(config)

    def test_init_updates_config_settings(self, config):
        settings = Settings()
        settings.preview = True
        settings.job.runscript.stdout_redirect = True
        settings.log.stdout = 5

        # When call init with additional settings
        init(config_settings=settings)

        # Then config settings set to default values with overrides
        assert_config_as_expected(config, preview=True, stdout_redirect=True, stdout=5)

    def test_init_passes_args_to_default_jobmanager(self, config, mock_jobmanager):
        # When call init with job manager args
        init(path="foo/bar", folder="test_folder", use_existing_folder=True)

        # Then default jobmanager initialised with the passed args
        assert mock_jobmanager.call_args_list[0].args == (
            {"counter_len": 3, "hashing": "input", "remove_empty_directories": True},
            "foo/bar",
            "test_folder",
            True,
        )

    def test_init_updates_init_flag_on_config_settings(self, config):
        # When call init with job manager args
        init()

        # Then config is marked as initialised
        assert config.init

    def test_init_idempotent(self, config, mock_jobmanager):
        # When call init twice, the second time with updated settings
        settings = Settings()
        settings.preview = True
        settings.job.runscript.stdout_redirect = True
        settings.log.stdout = 5
        init()
        init(config_settings=settings)

        # Then the second call is a no-op
        assert_config_as_expected(config)
        assert mock_jobmanager.call_count == 1

    def test_init_initialises_again_when_config_init_flag_reset(self, config, mock_jobmanager):
        # When call init twice, the second time with updated settings
        settings = Settings()
        settings.preview = True
        settings.job.runscript.stdout_redirect = True
        settings.log.stdout = 5
        init()
        config.init = False
        init(config_settings=settings)

        # Then the second call is effective
        assert_config_as_expected(config, preview=True, stdout_redirect=True, stdout=5)
        assert mock_jobmanager.call_count == 2

    @pytest.mark.parametrize(
        "code,version,tasks,should_error",
        [
            (0, "20.10", "16", False),
            (0, "21.10", "16", False),
            (1, "20.11", "16", True),
            (0, "10.11", "16", True),
            (0, "21_11", "16", True),
            (0, "21.10", None, True),
            (0, "20.10", "xyz", True),
        ],
        ids=[
            "happy_v20.10",
            "happy_v21.10",
            "unhappy_code",
            "unhappy_v10.11",
            "unhappy_v20_11",
            "unhappy_no_tasks",
            "unhappy_tasks",
        ],
    )
    def test_init_sets_slurm_settings(self, code, version, tasks, should_error, monkeypatch, config):
        # When running under slurm and init() called
        monkeypatch.setenv("SLURM_JOB_ID", "123456")
        if tasks is not None:
            monkeypatch.setenv("SLURM_TASKS_PER_NODE", tasks)
        else:
            monkeypatch.delenv("SLURM_TASKS_PER_NODE", raising=False)
        monkeypatch.setenv("SCM_SRUN_OPTIONS", "")

        mock_result = MagicMock()
        mock_stdout = MagicMock()
        mock_stdout.decode = lambda: f"x {version}"
        mock_result.returncode = code
        mock_result.stdout = mock_stdout

        with patch("scm.plams.core.functions.subprocess.run", return_value=mock_result):
            init()

            # Then config should have defaults set
            assert_config_as_expected(config)

            # And slurm settings set when happy, otherwise None
            if not should_error:
                assert config.slurm == Settings({"slurm_version": version.split("."), "tasks_per_node": [int(tasks)]})
                assert os.environ["SCM_SRUN_OPTIONS"].startswith("-m block:block:block,NoPack --use-min-nodes")
            else:
                assert config.slurm is None

    def test_finish_awaits_plams_threads(self, config, mock_jobmanager):
        # When set flag on settings object on a different thread
        config.init = True
        config.default_jobmanager = mock_jobmanager

        def set_flag():
            time.sleep(0.1)
            config.flag = True

        thread = threading.Thread(target=set_flag, name="plamsthread")
        thread.start()

        # Then thread awaited and flag therefore set
        finish()
        assert not thread.is_alive()
        assert config.flag

    def test_finish_calls_clean_on_job_managers(self, config):
        # When finish called with job manager arguments
        config.init = True
        jobmanagers = [MagicMock() for _ in range(3)]
        config.default_jobmanager = jobmanagers[0]
        finish(jobmanagers[1:])

        # Then clean called on both default and passed job managers
        for mock_jobmanager in jobmanagers:
            assert mock_jobmanager._clean.call_count == 1

    @pytest.mark.parametrize("erase_workdir", [True, False])
    def test_finish_respects_default_job_manager_erase_workdir_flag(self, erase_workdir, config, mock_jobmanager):
        with patch("shutil.rmtree") as mock_rmtree:
            # When erase_workdir configured
            config.init = True
            config.default_jobmanager = mock_jobmanager
            config.erase_workdir = erase_workdir
            finish()
            # Then rmtree only called if flag enabled
            assert mock_rmtree.call_count == (1 if erase_workdir else 0)

    def test_finish_updates_init_flag_on_config(self, config, mock_jobmanager):
        # When call finish
        config.init = True
        config.default_jobmanager = mock_jobmanager
        finish()

        # Then init flag always set to False
        assert not config.init

    def test_finish_idempotent(self, config, mock_jobmanager):
        # When call finish twice
        config.init = True
        config.default_jobmanager = mock_jobmanager
        finish()
        finish()

        # Then second call a no-op
        assert config.default_jobmanager._clean.call_count == 1

    def test_finish_cleans_again_when_config_flag_reset(self, config, mock_jobmanager):
        # When call finish twice
        config.init = True
        config.default_jobmanager = mock_jobmanager
        finish()
        config.init = True
        finish()

        # Then second call effective
        assert config.default_jobmanager._clean.call_count == 2

    def test_init_then_finish_as_expected(self, config, mock_jobmanager):
        # When call init then finish
        init()
        finish()

        # Then config defaults initialised and job manager created and cleaned
        assert_config_as_expected(config, init=False)
        assert config.default_jobmanager._clean.call_count == 1
        assert mock_jobmanager.call_args_list[0].args == (
            {"counter_len": 3, "hashing": "input", "remove_empty_directories": True},
            None,
            None,
            False,
        )

    def test_init_then_finish_successive_calls_as_expected(self, config, mock_jobmanager):
        # When call init then finish successively
        init()
        config.preview = True
        config.log.stdout = 3
        config.job.runscript.stdout_redirect = True
        config.foo = "bar"
        finish()
        init()
        finish()

        # Then on the second call the default settings should be reset, but the custom ones remain
        assert_config_as_expected(config, init=False)
        assert config.foo == "bar"
        assert config.default_jobmanager._clean.call_count == 2
        assert mock_jobmanager.call_args_list[0].args == (
            {"counter_len": 3, "hashing": "input", "remove_empty_directories": True},
            None,
            None,
            False,
        )

    def test_init_then_finish_in_loop_as_expected(self, config, mock_jobmanager):
        # When call init then finish in loop
        for i in range(3):
            init(folder=f"folder{i}")
            finish()

            # Then config defaults initialised and job manager created and cleaned
            assert_config_as_expected(config, init=False)
            assert config.default_jobmanager._clean.call_count == i + 1
            assert mock_jobmanager.call_args_list[i].args == (
                {"counter_len": 3, "hashing": "input", "remove_empty_directories": True},
                None,
                f"folder{i}",
                False,
            )
