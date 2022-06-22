from .scmjob import SCMJob, SCMResults
from ...core.settings import Settings
import os

__all__ = ['ParAMSJob', 'ParAMSResults']


class ParAMSResults(SCMResults):
    #_kfext = '.t61'
    #_rename_map = {'TAPE61':'$JN'+_kfext}

    def get_molecule(self, *args, **kwargs):
        raise PlamsError('ParAMSResults do not support get_molecule() method.')


class ParAMSJob(SCMJob):
    """A class for running ParAMS jobs.

    There are some convenient ways of setting up the input file:

    To set training set, validation set, job collection, engine collection, and parameter interface, use the special attributes

    * ``training_set`` and ``validation_set`` : these map to the first and second ``DataSet`` in the job settings. Note that when assigning to this variable you can use a str which will set the Path in the settings. Usage: ``job.training_set = 'path/to/training_set.yaml'``, OR ``job.training_set.Path = os.path.abspath('path/to/training_set.yaml')``, OR ``job.training_set = {'Path': os.path.abspath('training_set.yaml'), 'LoggerFrequency': {...}}``. In the former case the path is automatically converted to an absolute path.

    * ``job_collection``, ``engine_collection``, ``parameter_interface`` maps directly to the yaml file. Usage: ``job.job_collection = 'path/to/training_set.yaml'``. Here the path is automatically converted to an absolute path.

    Note: If you do not specify paths to the above files they will be read with
    the default names from the current working directory if they exist
    (converted to absolute paths in the input file).

    To add exit conditions, stoppers, or optimizers, you can use the
    ``add_exit_condition``, ``add_stopper``, and ``add_optimizer`` functions.
    They are more convenient to use than accessing the settings directly
    because of the repeated keys in the input file.

    
    """
    _result_type = ParAMSResults
    _command = 'params'
    _json_definitions = 'glompo'
    _subblock_end = 'End'
    _top = ['task', 'parameterinterface', 'dataset', 'jobcollection', 'enginecollection']

    def _ensure_training_set_in_settings(self):
        if 'DataSet' not in self.settings.input: #DataSet is already set if using .from_inputfile()
            self.settings.input.DataSet = [
                Settings({'Name': 'training_set'}),
            ]
    def _ensure_validation_set_in_settings(self):
        self._ensure_training_set_in_settings()
        if len(self.settings.input.DataSet) < 2:
            self.settings.input.DataSet.append(Settings({'Name': 'validation_set'}))

    def __init__(self, **kwargs):
        SCMJob.__init__(self,**kwargs)
        self.settings.ignore_molecule = True
        self._ensure_training_set_in_settings()

    @property
    def training_set(self):
        self._ensure_training_set_in_settings()
        return self.settings.input.DataSet[0]

    @training_set.setter
    def training_set(self, value):
        self._ensure_training_set_in_settings()
        if isinstance(value, str):
            self.settings.input.DataSet[0].path = self._abspath(value)
        elif isinstance(value, dict):
            self.settings.input.DataSet[0] = value
        else:
            raise TypeError(f"When assigning training set can only accept str (for the path) or a Settings instance. Tried {value} of type {type(value)}")

    @property
    def validation_set(self):
        self._ensure_validation_set_in_settings()
        return self.settings.input.DataSet[1]

    @validation_set.setter
    def validation_set(self, value):
        self._ensure_validation_set_in_settings()
        if isinstance(value, str):
            self.settings.input.DataSet[1].path = self._abspath(value)
        elif isinstance(value, dict):
            self.settings.input.DataSet[1] = value
        else:
            raise TypeError(f"When assigning validation set can only accept str (for the path) or a Settings instance. Tried {value} of type {type(value)}")


    @staticmethod
    def _abspath(value):
        if isinstance(value, str) and os.path.exists(value):
            return os.path.abspath(value).replace('\\', '/')
        return value.replace('\\', '/')

    @property
    def parameter_interface(self):
        return self.settings.input.ParameterInterface

    @parameter_interface.setter
    def parameter_interface(self, value):
        self.settings.input.ParameterInterface = self._abspath(value)

    @property
    def job_collection(self):
        return self.settings.input.JobCollection

    @job_collection.setter
    def job_collection(self, value):
        self.settings.input.JobCollection = self._abspath(value)

    @property
    def engine_collection(self):
        return self.settings.input.EngineCollection

    @engine_collection.setter
    def engine_collection(self, value):
        self.settings.input.EngineCollection = self._abspath(value)

    def _add_repeated_block(self, block, key, value):
        """
            Adds a repeated block to self.settings.input.

            Example: block='ExitCondition', key='MaxOptimizersConverged', value=1

            will give the following in the input file

            ```
            ExitCondition
                Type MaxOptimizersConverged
                MaxOptimizersConverged 1
            End
            ```
        """
        my_value = Settings(value) if isinstance(value, dict) else value
        if block in self.settings.input and not isinstance(self.settings.input[block], list):
            self.settings.input[block] = [self.settings.input[block]]
        if block not in self.settings.input:
            self.settings.input[block] = []
        s = Settings()
        s['Type'] = key
        s[key] = my_value
        self.settings.input[block].append(s)

    def add_exit_condition(self, key, value):
        self._add_repeated_block('ExitCondition', key, value)

    def add_stopper(self, key, value):
        self._add_repeated_block('Stopper', key, value)

    def add_optimizer(self, key, value=None):
        if value is None:
            value = Settings()
        self._add_repeated_block('Optimizer', key, value)

    def get_runscript(self):
        return super().get_runscript().replace('<', '-i ')

    def get_input(self, use_defaults=True, validate=True):
        """
            use_defaults : bool
                Try to locate some yaml files in the current working directory if they are not explicitly specified in the settings

            validate : bool
                Validate the written input file, calls ``validate_input``
        """

        original_settings = self.settings.copy()
        try:
            if use_defaults:
                if 'path' not in self.training_set and os.path.exists('training_set.yaml'):
                    self.training_set.path = os.path.abspath('training_set.yaml')
                if os.path.exists('validation_set.yaml'):
                    if (len(self.settings.input.DataSet) > 1 and 'path' not in self.validation_set) or len(self.settings.input.DataSet) <= 1:
                        self.validation_set.path = os.path.abspath('validation_set.yaml')
                if 'parameterinterface' not in self.settings.input and os.path.exists('parameter_interface.yaml'):
                    self.parameter_interface = 'parameter_interface.yaml'
                if 'parameterinterface' not in self.settings.input and os.path.exists('parameters.yaml'):
                    self.parameter_interface = 'parameters.yaml'
                if 'jobcollection' not in self.settings.input and os.path.exists('job_collection.yaml'):
                    self.job_collection = 'job_collection.yaml'
                if 'enginecollection' not in self.settings.input and os.path.exists('engine_collection.yaml'):
                    self.engine_collection = 'engine_collection.yaml'
                if 'enginecollection' not in self.settings.input and os.path.exists('job_collection_engines.yaml'):
                    self.engine_collection = 'job_collection_engines.yaml'
                if 'exitcondition' not in self.settings.input:
                    self.add_exit_condition('MaxOptimizersConverged', 1)

            inp = super().get_input()
        except Exception as e:
            self.settings = original_settings.copy()
            raise e
        self.settings = original_settings.copy()
        if validate:
            self.validate_input(inp)
        return inp

    def validate_input(self, inp):
        try:
            from scm.input_parser import InputParser
        except ImportError:  # Try to load the parser from $AMSHOME/scripting
            with UpdateSysPath():
                from scm.input_parser import InputParser
        with InputParser() as parser:
            s = parser.to_settings(self._json_definitions, inp)

    def collect(self):
        pass

    @classmethod
    def from_inputfile(cls, path):
        """
            Initializes a ParAMSJob with settings taken from an input file (e.g. params.in)

            Paths are replaced with absolute paths to enable the running of new jobs
        """

        job = super().from_inputfile(path)
        sett_dir = os.path.abspath(os.path.dirname(path))
        if 'DataSet' in job.settings.input:
            # I hope that the Path is required when specifying a data set and doesn't have default value
            if not os.path.isabs(job.training_set.Path):
                job.training_set.Path = os.path.join(sett_dir, job.training_set.Path).replace('\\', '/')
            if len(job.settings.input.DataSet) >= 2:
                if not os.path.isabs(job.validation_set.Path):
                    job.validation_set.Path = os.path.join(sett_dir, job.validation_set.Path).replace('\\', '/')
        if 'EngineCollection' in job.settings.input and not os.path.isabs(job.engine_collection):
            job.engine_collection = os.path.join(sett_dir, job.engine_collection).replace('\\', '/')
        if 'JobCollection' in job.settings.input and not os.path.isabs(job.job_collection):
            job.job_collection = os.path.join(sett_dir, job.job_collection).replace('\\', '/')
        if 'ParameterInterface' in job.settings.input and not os.path.isabs(job.parameter_interface):
            job.parameter_interface = os.path.join(sett_dir, job.parameter_interface).replace('\\', '/')

        return job

    @classmethod
    def load_external(cls, path, name=None):
        assert os.path.isdir(path), f"Path {path} does not exist or is not a directory"
        path = os.path.abspath(path)
        sett_dir = os.path.join(path, 'settings_and_initial_data')
        assert os.path.isdir(sett_dir), f"Directory {sett_dir} does not exist. settings_and_initial_data needs to be a subdirectory of {path}"
        params_in = os.path.join(sett_dir,  'params.in')
        if os.path.exists(params_in):
            job = ParAMSJob.from_inputfile(params_in)
        else:
            job = ParAMSJob()

        if name:
            job.name = name

        return job







