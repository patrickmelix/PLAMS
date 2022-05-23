"""Implementation of the AMSPipeCalculator class.

"""
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree, Bohr

from .amsworker import AMSWorker
from .ams import AMSJob
from ...core.settings import Settings
from ..molecule.ase import fromASE, toASE

__all__ = ['AMSCalculator', 'BasePropertyExtractor'] 


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
        settings.input.ams.Properties.Gradients = True


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
        settings.input.ams.Properties.StressTensor = True


class AMSCalculator(Calculator):
    """
    ASE Calculator which runs the AMS Driver. 

    Parameters:

    settings : Settings
               A Settings object representing the input for an AMSJob or AMSWorker.
    molecule : Molecule , optional
               A Molecule object for which the calculation has to be performed. If 
               settings.input.ams.system is defined it overrides the molecule argument.
               If AMSCalculator.calculate(atoms = atoms) is called with an atoms argument 
               it overrides any earlier definition of the system and remembers it.

    Additional settings for the calculator can be passed in settings.Calculator:

    settings.Calculator.Pipe: use the AMSWorker to set up an interactive session.
    settings.Calculator.Pipe.worker: keyword arguments for the AMSWorker.

    settings.Calculator.AMSJob: use the AMSJob to set up an io session.

    settings.Calculator.Extractors: define extractors for additional properties.
    

    Examples:
    """


    def __new__(cls, settings = None, molecule = None):
        """Dispatch object creation to AMSPipeCalculator or AMSJobCalculator depending on |settings|"""
        if cls == AMSCalculator:
            if 'Pipe' in settings.Calculator:
                obj = object.__new__(AMSPipeCalculator)
            elif 'AMSJob' in settings.Calculator:
                obj = object.__new__(AMSJobCalculator)
            else:
                raise NotImplementedError("Only settings.Calculator.Pipe and .AMSJob are implemented.")
        else:
            obj = object.__new__(cls)
        return obj

    def __init__(self, settings, molecule = None):
        settings = settings.copy()
        self.settings = settings.copy()
        self.molecule = molecule
        self.calculator_settings = settings.Calculator.copy()
        del settings.Calculator

        if 'system' in self.settings.input.ams:
            mol_dict = AMSJob.settings_to_mol(settings)
            atoms = toASE(mol_dict['']) if mol_dict else None
        elif molecule:
            atoms = toASE(molecule)
        else:
            atoms = None

        super().__init__()
        self.atoms = atoms
        self.ams_settings = settings
        
        self.counter = 0
        self.prev_ams_results = None
        self.results = dict()

        self.extractors = [EnergyExtractor(), ForceExtractor(), StressExtractor()]
        self.extractors += self.calculator_settings.get('extractors', [])

    def __getnewargs__(self):
        return self.settings, self.molecule

    @property
    def implemented_properties(self):
        """Returns the list of properties that this calculator has implemented"""
        return [ extractor.name for extractor in self.extractors ]

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """Calculate the requested properties. If atoms is not set, it will reuse the last known Atoms object."""
        if atoms is not None:
            self.atoms = atoms

        if self.atoms is None:
            raise ValueError("No atoms object was set.")

        if len(system_changes) == 0:
            return

        molecule = fromASE(self.atoms)
        self.counter += 1
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
        settings = self.ams_settings.copy()
        for extractor in self.extractors:
            if extractor.name in properties:
                extractor.set_settings(settings)
        return settings


class AMSPipeCalculator(AMSCalculator):
    """This class should be instantiated through AMSCalculator with settings.Calculator.Pipe defined"""
    def __init__(self, settings, molecule = None):
        super().__init__(settings)

        worker_settings = self.ams_settings.copy()
        del worker_settings.input.ams.Task

        worker_kwargs = settings.Calculator.Pipe.Worker.as_dict()
        self.worker = AMSWorker(worker_settings, **worker_kwargs)

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
        return self.worker._solve_from_settings(name     = str(self.counter),
                                                molecule = molecule,
                                                settings = job_settings
                                               )


class AMSJobCalculator(AMSCalculator):
    """This class should be instantiated through AMSCalculator with settings.Calculator.AMSJob defined"""
    def _get_ams_results(self, molecule, properties):
        settings = self.ams_settings.copy()
        job_settings = self._get_job_settings(properties)
        return AMSJob( name     = str(self.counter),
                       molecule = molecule, 
                       settings = job_settings
                       ).run()

    
