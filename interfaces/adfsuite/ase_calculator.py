"""Implementation of the AMSPipeCalculator class.

"""
import numpy as np

from .amsworker import AMSWorker
from .ams import AMSJob
from ...core.settings import Settings
from ..molecule.ase import fromASE, toASE

from scm.plams import log
__all__ = ['AMSCalculator', 'BasePropertyExtractor'] 

try:
    from ase.calculators.calculator import Calculator, all_changes
    from ase.units import Hartree, Bohr
except ImportError:
    #empty interface if ase does not exist:
    __all__ = []
    class Calculator:
        def __init__(self, *args, **kwargs):
            raise NotImplementedError('AMSCalculator can not be used without ASE')
    all_changes = []
    Hartree = Bohr = 1


class BasePropertyExtractor:
    def __call__(self, ams_results, atoms):
        return self.extract(ams_results, atoms)

    def extract(self, ams_result, atoms):
        raise NotImplementedError

    def set_settings(self, settings):
        pass

    def check_settings(self, settings):
        s = settings.copy()
        self.set_settings(s)
        return s == settings


class EnergyExtractor(BasePropertyExtractor):
    name = 'energy'
    def extract(self, ams_results, atoms):
        return ams_results.get_energy() * Hartree


class ForceExtractor(BasePropertyExtractor):
    name = 'forces'
    def extract(self, ams_results, atoms):
        return -ams_results.get_gradients() * Hartree / Bohr

    def set_settings(self, settings):
        settings.input.ams.Properties.Gradients = "Yes"


class StressExtractor(BasePropertyExtractor):
    name = 'stress'
    def extract(self, ams_results, atoms):
            D = sum(atoms.get_pbc())
            if isinstance(atoms.get_pbc(), bool): D *= 3
            st = ams_results.get_stresstensor() * Hartree / Bohr**D
            xx, yy, zz, yz, xz, xy = (0., 0., 0., 0., 0., 0.)
            if D >= 1:
                xx = st[0][0]
                if D >= 2:
                    yy = st[1][1]
                    xy = st[0][1]
                    if D >= 3:
                        zz = st[2][2]
                        yz = st[1][2]
                        xz = st[0][2]
            return np.array([xx, yy, zz, yz, xz, xy])

    def set_settings(self, settings):
        settings.input.ams.Properties.StressTensor = "Yes"


class AMSCalculator(Calculator):
    """
    ASE Calculator which runs the AMS Driver. 

    Parameters:

    settings  : Settings
                A Settings object representing the input for an AMSJob or AMSWorker.
    name      : str, optional
                Name of the rundir of calculations done by this calculator. A counter 
                is appended to the name for every calculation.
    amsworker : bool , optional
                If True, use the AMSWorker to set up an interactive session. Otherwise
                use AMSJob to set up an io session.
    restart   : bool , optional
                Allow the engine to restart based on previous calculations.
    molecule  : Molecule , optional
                A Molecule object for which the calculation has to be performed. If 
                settings.input.ams.system is defined it overrides the molecule argument.
                If AMSCalculator.calculate(atoms = atoms) is called with an atoms argument 
                it overrides any earlier definition of the system and remembers it.
    extractors: List[BasePropertyExtractor] , optional
                Define extractors for additional properties.
    

    Examples:
    """
    #counters are a dict as a class variable. This is to support deepcopying/multiple instances with the same name
    _counter = {}

    def __new__(cls, settings = None, name = '', amsworker = False, restart = True, molecule = None, extractors = []):
        """Dispatch object creation to AMSPipeCalculator or AMSJobCalculator depending on |amsworker|"""
        if cls == AMSCalculator:
            if amsworker:
                obj = object.__new__(AMSPipeCalculator)
            else:
                obj = object.__new__(AMSJobCalculator)
        else:
            obj = object.__new__(cls)
        return obj

    def __init__(self, settings = None, name='', amsworker = False, restart = True, molecule = None, extractors = []):
        settings = settings.copy()
        self.settings = settings.copy()
        self.amsworker = amsworker
        self.name = name
        self.restart = restart
        self.molecule = molecule
        extractors = settings.pop('Extractors', [])
        self.extractors = [EnergyExtractor(), ForceExtractor(), StressExtractor()]
        self.extractors += [e for e in extractors if not e in self.extractors]
        
        if 'system' in self.settings.input.ams:
            mol_dict = AMSJob.settings_to_mol(settings)
            atoms = toASE(mol_dict['']) if mol_dict else None
        elif molecule:
            atoms = toASE(molecule)
        else:
            atoms = None

        super().__init__()
        self.atoms = atoms

        self.prev_ams_results = None
        self.results = dict()

    @property
    def counter(self):
        #this is needed for deepcopy/pickling etc
        if not self.name in self._counter:
            self.set_counter()
        self._counter[self.name] += 1
        return self._counter[self.name]

    def set_counter(self, value = 0):
        self._counter[self.name] = value

    #def __getnewargs__(self):
    #    return self.settings, self.name, self.amsworker, self.restart, self.molecule, self.extractors

    @property
    def implemented_properties(self):
        """Returns the list of properties that this calculator has implemented"""
        return [ extractor.name for extractor in self.extractors ]

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """Calculate the requested properties. If atoms is not set, it will reuse the last known Atoms object."""
        log("I was asked to compute "+str( properties)+ " for the following system", 0)
        log(str(atoms), 0)
        if atoms is not None:
            #no need to redo the calculation, we already have everything.
            if self.atoms == atoms and all([p in self.results for p in properties]):
                return
            self.atoms = atoms

        if self.atoms is None:
            raise ValueError("No atoms object was set.")

        if len(system_changes) == 0:
            return

        molecule = fromASE(self.atoms)
        ams_results = self._get_ams_results(molecule, properties)
        if not ams_results.ok():
            self.results = dict()
            return

        self.results_from_ams_results(ams_results, self._get_job_settings(properties))
        self.prev_ams_results = ams_results

    def results_from_ams_results(self, ams_results, job_settings):
        """Populates the self.results dictionary by having extractors act on an AMSResults object."""
        for extractor in self.extractors:
            if extractor.check_settings(job_settings):
                self.results[extractor.name] = extractor.extract(ams_results, self.atoms)

    def _get_ams_results(self, molecule, properties):
        raise NotImplementedError("Subclasses of AMSCalculator should implement this")

    def _get_job_settings(self, properties):
        """Returns a Settings object which ensures that an AMS calculation is run from which all requested
        properties can be extracted"""
        settings = self.settings.copy()
        for extractor in self.extractors:
            if extractor.name in properties:
                extractor.set_settings(settings)
        return settings


class AMSPipeCalculator(AMSCalculator):
    """This class should be instantiated through AMSCalculator with settings.Calculator.Pipe defined"""
    def __init__(self, settings = None, name = '', amsworker = False, restart = True, molecule = None, extractors = []):
        super().__init__(settings, name, amsworker, restart, molecule, extractors)

        worker_settings = self.settings.copy()
        del worker_settings.input.ams.Task

        self.worker = AMSWorker(worker_settings, use_restart_cache = self.restart)

    def _get_ams_results(self, molecule, properties):
        job_settings = self._get_job_settings(properties)
        
        #AMSWorker expects no engine definition at this point.
        s = Settings()
        s.input.ams = job_settings.input.ams
        if 'amsworker' in job_settings:
            s.amsworker = job_settings.amsworker
        s.amsworker.prev_results = self.prev_ams_results
        job_settings = s
        #run the worker from _solve_from_settings
        return self.worker._solve_from_settings(name     = self.name+str(self.counter),
                                                molecule = molecule,
                                                settings = job_settings
                                               )


class AMSJobCalculator(AMSCalculator):
    """This class should be instantiated through AMSCalculator with settings.Calculator.AMSJob defined"""
    def _get_ams_results(self, molecule, properties):
        settings = self.settings.copy()
        job_settings = self._get_job_settings(properties)
        if self.restart and self.prev_ams_results:
            job_settings.input.ams.EngineRestart = self.prev_ams_results.rkfpath(file='engine')

        return AMSJob( name     = self.name+str(self.counter),
                       molecule = molecule, 
                       settings = job_settings
                       ).run()

    
