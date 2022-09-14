import os
import shutil
from scm.plams import Settings, JobError, AMSJob, CRSJob, Molecule, AMSResults, CRSResults, KFFile
from collections import OrderedDict
from ..core.functions import add_to_instance
from ..core.basejob import MultiJob
from ..core.results import Results
from ..core.settings import Settings
from ..mol.molecule import Molecule
from ..mol.atom import Atom
from ..interfaces.adfsuite.ams import AMSJob
from ..tools.units import Units
from ..interfaces.adfsuite.quickjobs import model_to_settings
import numpy as np

__all__ = ['ADFCOSMORSCompoundJob', 'ADFCOSMORSCompoundResults']

class ADFCOSMORSCompoundResults(Results):
    """Results class for ADFCOSMORSCompoundJob
    """
    def coskfpath(self):
        """
            Returns the path to the resulting .coskf
        """

        return os.path.join(self.job.path, self.job.name+".coskf")

    def get_main_molecule(self):
        """
            Returns the optimized molecule
        """

        return self.job.children['solv'].results.get_main_molecule()

    def get_input_molecule(self):
        """
            Returns the input molecule
        """
        for job in self.job.children.values():
            return job.results.get_input_molecule()


class ADFCOSMORSCompoundJob(MultiJob):
    """A class for performing the equivalent of Task COSMO-RS Compound in the AMS GUI
    """

    _result_type = ADFCOSMORSCompoundResults

    def __init__(self, molecule:Molecule, preoptimization=None, singlepoint=False, settings=None, **kwargs):
        """

            Class for running the equivalent of "COSMO-RS Compound" in the AMS
            GUI. Note that these are ADF calculations, not COSMO-RS
            calculations!

            Initialize two or three jobs:

            (optional): Preoptimization with force field or semi-empirical method
            1. Gasphase optimization (BP86, DZP)
            2. Gasphase optimization (BP86, TZP, BeckeGrid Quality Good)
            3. Take optimized structure and run singlepoint with implicit solvation

            Access the result .coskf file with ``job.results.coskfpath()``.
            Note: this file will be called jobname.coskf, where jobname is the
            name of the ADFCOSMORSCompoundJob.

            molecule  : Molecule
                The initial structure

            preoptimization : str or None
                If None, do not preoptimize with a fast engine (then initial optimization is done with ADF). Otherwise, can be one of 'UFF', 'GAFF', 'GFNFF', 'GFN1-xTB', 'ANI-2x'. Note that you need valid licenses for ForceField or DFTB or MLPotential to use these preoptimizers.

            singlepoint : bool
                Only run a singlepoint with solvation to generate the .coskf file on the given Molecule. (no geometry optimization)

            settings : settings
                settings.runscript.nproc, settings.input.adf.custom_options. If 'adf' is in settings.input it should be provided without the solvation block.

            Usage:

            .. code-block:: python
                
                mol = from_smiles('O')
                job = ADFCOSMORSCompoundJob(name='water', molecule=mol, preoptimization='UFF')
                job.run()
                print(job.results.coskfpath())

        """
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)
        self.input_molecule = molecule
        self.settings = settings or Settings()

        if not singlepoint:
            if preoptimization:
                preoptimization_s = Settings()
                preoptimization_s.runscript.nproc = 1
                preoptimization_s.input.ams.Task = 'GeometryOptimization'
                preoptimization_s += model_to_settings(preoptimization)
                preoptimization_job = AMSJob(settings=preoptimization_s, name='preoptimization', molecule=molecule)
                self.children['preoptimization'] = preoptimization_job

            gas_s = Settings()
            gas_s.input.ams.Task = 'GeometryOptimization'
            gas_s += self.adf_settings(solvation=False, settings=self.settings)
            gas_job = AMSJob(settings=gas_s, name='gas')

            if preoptimization:
                @add_to_instance(gas_job)
                def prerun(self):
                    self.molecule = self.parent.children['preoptimization'].results.get_main_molecule()
            else:
                gas_job.molecule = molecule

            self.children['gas'] = gas_job

        solv_s = Settings()
        solv_s.input.ams.Task = 'SinglePoint'
        solv_s += self.adf_settings(solvation=True, settings=self.settings)
        solv_job = AMSJob(settings=solv_s, name='solv')

        if singlepoint:
            @add_to_instance(solv_job)
            def prerun(self):
                self.molecule = self.parent.input_molecule
        else:
            @add_to_instance(solv_job)
            def prerun(self):
                gas_job.results.wait()
                self.settings.input.ams.EngineRestart = "../gas/adf.rkf" 
                self.settings.input.ams.LoadSystem.File = "../gas/ams.rkf"
                #self.settings.input.ams.EngineRestart = self.parent.children['gas'].results.rkfpath(file='adf') # this doesn't work with PLAMS restart since the file will refer to the .res directory (so the job is rerun needlessly)
                #self.settings.input.ams.LoadSystem.File = self.parent.children['gas'].results.rkfpath(file='ams')
                # cannot copy to gasphase-ams.rkf etc. because that conflicts with PLAMS restarts
                #shutil.copyfile(gas_job.results.rkfpath(file='ams'), os.path.join(self.path, 'gasphase-ams.rkf'))
                #shutil.copyfile(gas_job.results.rkfpath(file='adf'), os.path.join(self.path, 'gasphase-adf.rkf'))

        @add_to_instance(solv_job)
        def postrun(self):
            self.parent.convert_to_coskf(self.results.rkfpath(file='adf'), os.path.join(self.parent.path, self.parent.name+'.coskf'))

        self.children['solv'] = solv_job

    @staticmethod
    def solvation_settings() -> Settings:
        sett = Settings()
        sett.input.adf.solvation = {
            'surf': 'Delley',
            'solv': 'name=CRS cav0=0.0 cav1=0.0',
            'charged': 'method=Conj corr',
            'c-mat': 'Exact',
            'scf': 'Var All',
            'radii': {
                'H': 1.30,
                'C': 2.00,
                'N': 1.83,
                'O': 1.72,
                'F': 1.72,
                'Si': 2.48,
                'P': 2.13,
                'S': 2.16,
                'Cl': 2.05,
                'Br': 2.16,
                'I': 2.32
            }
        }
        return sett


    @staticmethod
    def adf_settings(solvation:bool, settings=None) -> Settings:
        """
        Returns ADF settings with or without solvation

        If solvation == True, then also include the solvation block.
        """
        

        s = Settings()
        if settings:
            s = settings.copy()
        if 'adf' not in s.input:
            s.input.adf.Basis.Type = 'TZP'
            s.input.adf.Basis.Core = 'Small'
            s.input.adf.XC.GGA = 'BP86'
            s.input.adf.Symmetry = 'NOSYM'
            s.input.adf.BeckeGrid.Quality = 'Good'
        if solvation:
            s += ADFCOSMORSCompoundJob.solvation_settings()

        return s

    @staticmethod
    def convert_to_coskf(rkf: str, coskf: str):
        """ rkf: absolute path to adf.rkf, coskf: path to write out the resulting .coskf file  """
        f = KFFile(rkf)
        cosmo = f.read_section("COSMO")
        coskf = KFFile(coskf)
        for k,v in cosmo.items():
            coskf.write("COSMO",k,v)



