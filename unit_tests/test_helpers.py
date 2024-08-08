import builtins
from scm.plams.core.settings import (
    SafeRunSettings,
    LogSettings,
    RunScriptSettings,
    JobSettings,
    JobManagerSettings,
    ConfigSettings,
)

real_import = builtins.__import__


def get_mock_import_function(import_error_name):
    """
    Gets a mock version of the import function, which will raise an ImportError for the specified module,
    but import all other modules as usual.
    """

    def import_with_error(name, globals=None, locals=None, fromlist=(), level=0):
        if name == import_error_name:
            raise ImportError(f"Mock import error for '{import_error_name}'")
        return real_import(name, globals=globals, locals=locals, fromlist=fromlist, level=level)

    return import_with_error


def assert_config_as_expected(
    config, init=True, preview=False, stdout_redirect=False, stdout=3, verify_derived_types=True
):
    """
    Assert the passed config object has the expected values.
    These are the default values, unless overridden via the args.
    """
    assert config.init == init
    assert config.preview == preview
    assert config.sleepstep == 5
    assert config.ignore_failure
    assert config.daemon_threads
    assert not config.erase_workdir
    assert config.jobmanager.counter_len == 3
    assert config.jobmanager.hashing == "input"
    assert config.jobmanager.remove_empty_directories
    assert config.job.pickle
    assert config.job.pickle_protocol == -1
    assert config.job.keep == "all"
    assert config.job.save == "all"
    assert config.job.runscript.shebang == "#!/bin/sh"
    assert config.job.runscript.stdout_redirect == stdout_redirect
    assert config.job.link_files
    assert config.log.file == 5
    assert config.log.stdout == stdout
    assert config.log.time
    assert config.log.date
    assert config.saferun.repeat == 10
    assert config.saferun.delay == 1

    if verify_derived_types:
        assert isinstance(config.saferun, SafeRunSettings)
        assert isinstance(config.log, LogSettings)
        assert isinstance(config.job, JobSettings)
        assert isinstance(config.job.runscript, RunScriptSettings)
        assert isinstance(config.jobmanager, JobManagerSettings)
        assert isinstance(config, ConfigSettings)
