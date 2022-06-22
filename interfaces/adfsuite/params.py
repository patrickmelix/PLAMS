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
    """
    _result_type = ParAMSResults
    _command = 'params'
    _json_definitions = 'glompo'
    _subblock_end = 'End'
    _top = ['task', 'parameterinterface', 'dataset', 'jobcollection', 'enginecollection']

    def __init__(self, **kwargs):
        SCMJob.__init__(self,**kwargs)
        self.settings.ignore_molecule = True
        self.settings.input.DataSet = [
            Settings({'Name': 'training_set'}),
        ]

    @property
    def training_set(self):
        return self.settings.input.DataSet[0]

    @property
    def validation_set(self):
        if len(self.settings.input.DataSet) < 2:
            self.settings.input.DataSet.append(Settings({'Name': 'validation_set'}))
        return self.settings.input.DataSet[1]

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

    def add_optimizer(self, key, value):
        self._add_repeated_block('Optimizer', key, value)

    def get_runscript(self):
        return super().get_runscript().replace('<', '-i ')

    def get_input(self, validate=True):
        inp = super().get_input()
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







