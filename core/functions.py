import glob
import os
import shutil
import sys
import threading
import time
import types

from os.path import join as opj
from os.path import isfile, isdir, expandvars, dirname

from .errors import PlamsError
from .settings import Settings

__all__ = ['init', 'finish', 'log', 'load', 'load_all', 'add_to_class', 'add_to_instance', 'config', 'read_molecules']

config = Settings()

#===========================================================================


def init(path=None, folder=None):
    """Initialize PLAMS environment. Create global ``config`` and the default |JobManager|.

    An empty |Settings| instance is created and populated with default settings by executing ``plams_defaults``. The following locations are used to search for the defaults file, in order of precedence:

    *   If ``$PLAMSDEFAULTS`` variable is in your environment and it points to a file, this file is used (executed as a Python script).
    *   If ``$ADFHOME`` variable is in your environment and ``$ADFHOME/scripting/scm/plams/plams_defaults`` exists, it is used.
    *   Otherwise, the path ``../plams_defaults`` relative to the current file (``functions.py``) is checked. If defaults file is not found there, an exception is raised.

    Then a |JobManager| instance is created as ``config.default_jobmanager`` using *path* and *folder* to determine the main working folder. Settings for this instance are taken from ``config.jobmanager``. If *path* is not supplied, the current directory is used. If *folder* is not supplied, ``plams_workdir`` is used.

    .. warning::
      This function **must** be called before any other PLAMS command can be executed. Trying to do anything without it results in a crash. See also |master-script|.
    """

    if 'PLAMSDEFAULTS' in os.environ and isfile(expandvars('$PLAMSDEFAULTS')):
        defaults = expandvars('$PLAMSDEFAULTS')
    elif 'ADFHOME' in os.environ and isfile(opj(expandvars('$ADFHOME'), 'scripting', 'scm', 'plams', 'plams_defaults')):
        defaults = opj(expandvars('$ADFHOME'), 'scripting', 'scm', 'plams', 'plams_defaults')
    else:
        defaults = opj(dirname(dirname(__file__)), 'plams_defaults')
        if not isfile(defaults):
            raise PlamsError('plams_defaults not found, please set PLAMSDEFAULTS or ADFHOME in your environment')
    exec(compile(open(defaults).read(), defaults, 'exec'))

    from .jobmanager import JobManager
    config.default_jobmanager = JobManager(config.jobmanager, path, folder)

    log('Running PLAMS located in {}'.format(dirname(dirname(__file__))) ,5)
    log('Using Python {}.{}.{} located in {}'.format(*sys.version_info[:3], sys.executable), 5)
    log('PLAMS defaults were loaded from {}'.format(defaults) ,5)

    log('PLAMS environment initialized', 5)
    log('PLAMS working folder: {}'.format(config.default_jobmanager.workdir), 1)

    try:
        import dill
    except ImportError:
        log('WARNING: importing dill package failed. Falling back to the default pickle module. Expect problems with pickling', 1)


#===========================================================================


def finish(otherJM=None):
    """Wait for all threads to finish and clean the environment.

    This function must be called at the end of your script for |cleaning| to take place. See |master-script| for details.

    If you used some other job managers than just the default one, they need to be passed as *otherJM* list.
    """
    for thread in threading.enumerate():
        if thread.name == 'plamsthread':
            thread.join()

    config.default_jobmanager._clean()
    if otherJM:
        for jm in otherJM:
            jm._clean()
    log('PLAMS environment cleaned up successfully', 5)
    log('PLAMS run finished. Goodbye', 3)

    if config.erase_workdir is True:
        shutil.rmtree(config.default_jobmanager.workdir)


#===========================================================================


def load(filename):
    """Load previously saved job from ``.dill`` file. This is just a shortcut for |load_job| method of the default |JobManager| ``config.default_jobmanager``."""
    return config.default_jobmanager.load_job(filename)


#===========================================================================


def load_all(path, jobmanager=None):
    """Load all jobs from *path*.

    This function works as multiple executions of |load_job|. It searches for ``.dill`` files inside the directory given by *path*, yet not directly in it, but one level deeper. In other words, all files matching ``path/*/*.dill`` are used. That way a path to the main working folder of a previously run script can be used to import all the jobs run by that script.

    In case of partially failed |MultiJob| instances (some children jobs finished successfully, but not all) the function will search for ``.dill`` files in children folders. That means, if ``path/[multijobname]/`` contains some subfolders (for children jobs) but does not contail a ``.dill`` file (the |MultiJob| was not fully successful), it will look into these subfolders. This behavior is recursive up to any folder tree depth.

    The purpose of this function is to provide a quick way of restarting a script. Loading all successful jobs from the previous run prevents double work and allows the new execution of the script to proceed directly to the place where the previous execution failed.

    Jobs are loaded using default job manager stored in ``config.default_jobmanager``. If you wish to use a different one you can pass it as *jobmanager* argument of this function.

    Returned value is a dictionary containing all loaded jobs as values and absolute paths to ``.dill`` files as keys.
    """
    jm = jobmanager or config.default_jobmanager
    loaded_jobs = {}
    for foldername in filter(lambda x: isdir(opj(path,x)), os.listdir(path)):
        maybedill = opj(path,foldername,foldername+'.dill')
        if isfile(maybedill):
            job = jm.load_job(maybedill)
            if job:
                loaded_jobs[os.path.abspath(maybedill)] = job
        else:
            loaded_jobs.update(load_all(path=opj(path,foldername), jobmanager=jm))
    return loaded_jobs


#===========================================================================


def read_molecules(folder, formats=None):
    """Read all molecules from *folder*.

    Read all the files present in *folder* with extensions compatible with :meth:`Molecule.read<scm.plams.core.basemol.Molecule.read>`. Returned value is a dictionary with keys being molecule names (filename without extension) and values being |Molecule| instances.

    The optional argument *formats* can be used to narrow down the search to files with specified extensions::

        molecules = read_molecules('mymols', formats=['xyz', 'pdb'])

    """
    from .basemol import Molecule
    extensions = formats or list(Molecule._readformat.keys())
    is_valid = lambda x: isfile(opj(folder,x)) and any([x.endswith('.'+ext) for ext in extensions])
    filenames = filter(is_valid, os.listdir(folder))
    ret = {}
    for f in filenames:
        m = Molecule(opj(folder,f))
        ret[m.properties.name] = m
    return ret


#===========================================================================


_stdlock = threading.Lock()
_filelock = threading.Lock()

def log(message, level=0):
    """Log *message* with verbosity *level*.

    Logs are printed independently to the text logfile (a file called ``logfile`` in the main working folder) and to the standard output. If *level* is equal or lower than verbosity (defined by ``config.log.file`` or ``config.log.stdout``) the message is printed. Date and/or time can be added based on ``config.log.date`` and ``config.log.time``. All logging activity is thread safe.
    """
    if 'log' in config:
        if level <= config.log.file or level <= config.log.stdout:
            message = str(message)
            prefix = ''
            if config.log.date:
                prefix += '%d.%m|'
            if config.log.time:
                prefix += '%H:%M:%S'
            if prefix:
                prefix = '[' + prefix.rstrip('|') + '] '
                message = time.strftime(prefix) + message
            if level <= config.log.stdout:
                with _stdlock:
                    print(message)
            if level <= config.log.file and 'default_jobmanager' in config:
                with _filelock, open(config.default_jobmanager.logfile, 'a') as f:
                    f.write(message + '\n')


#===========================================================================


def add_to_class(classname):
    """Add decorated function as a method to the whole class *classname*.

    The decorated function should follow a method-like syntax, with the first argument ``self`` that references the class instance.
    Example usage::

        @add_to_class(ADFResults)
        def get_energy(self):
            return self.readkf('Energy', 'Bond Energy')

    After executing the above code all instances of ``ADFResults`` in the current script (even the ones created beforehand) are enriched with ``get_energy`` method that can be invoked by::

        someadfresults.get_energy()

    The added method is accessible also from subclasses of *classname* so ``@add_to_class(Results)`` in the above example will work too.

    If *classname* is |Results| or any of its subclasses, the added method will be wrapped with the thread safety guard (see |parallel|).
    """
    from .results import _restrict, _MetaResults
    def decorator(func):
        if isinstance(classname, _MetaResults):
            func = _restrict(func)
        setattr(classname, func.__name__, func)
    return decorator


#===========================================================================


def add_to_instance(instance):
    """Add decorated function as a method to one particular *instance*.

    The decorated function should follow a method-like syntax, with the first argument ``self`` that references the class instance.
    Example usage::

        results = myjob.run()

        @add_to_instance(results)
        def get_energy(self):
            return self.readkf('Energy', 'Bond Energy')

        results.get_energy()

    The added method is accessible only for that one particular instance and it overrides any methods with the same name defined on a class level (in original class' source) or added with :func:`add_to_class` decorator.

    If *instance* is an instance of |Results| or any of its subclasses, the added method will be wrapped with the thread safety guard (see |parallel|).
    """
    from .results import _restrict, Results
    def decorator(func):
        if isinstance(instance, Results):
            func = _restrict(func)
        func = types.MethodType(func, instance)
        setattr(instance, func.__func__.__name__, func)
    return decorator
