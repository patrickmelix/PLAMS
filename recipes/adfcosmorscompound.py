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
        try: 
            return self.job.children['preoptimization'].results.get_input_molecule()
        except KeyError:
            return self.job.children['dft_preoptimization'].results.get_input_molecule()


class ADFCOSMORSCompoundJob(MultiJob):
    """A class for performing the equivalent of Task COSMO-RS Compound in the AMS GUI
    """

    _result_type = ADFCOSMORSCompoundResults

    def __init__(self, molecule:Molecule, preoptimization=None, **kwargs):
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

            Usage:

            .. code-block:: python
                
                mol = from_smiles('O')
                job = ADFCOSMORSCompoundJob(name='water', molecule=mol, preoptimization='UFF')
                job.run()
                print(job.results.coskfpath())

        """
        MultiJob.__init__(self, children=OrderedDict(), **kwargs)

        if preoptimization:
            preoptimization_s = Settings()
            preoptimization_s.runscript.nproc = 1
            preoptimization_s.input.ams.Task = 'GeometryOptimization'
            preoptimization_s += model_to_settings(preoptimization)
            preoptimization_job = AMSJob(settings=preoptimization_s, name='preoptimization', molecule=molecule)
            self.children['preoptimization'] = preoptimization_job

        dft_preoptimization_s = Settings()
        dft_preoptimization_s.input.ams.Task = 'GeometryOptimization'
        dft_preoptimization_s.input.adf.XC.GGA = 'BP86'
        dft_preoptimization_s.input.adf.Basis.Type = 'DZP'
        dft_preoptimization_job = AMSJob(settings=dft_preoptimization_s, name='dft_preoptimization', molecule=molecule)
        if preoptimization:
            @add_to_instance(dft_preoptimization_job)
            def prerun(self):
                self.molecule = self.parent.children['preoptimization'].results.get_main_molecule()
        self.children['dft_preoptimization'] = dft_preoptimization_job

        gas_s = Settings()
        gas_s.input.ams.Task = 'GeometryOptimization'
        gas_s += self.adf_settings(solvation=False)
        gas_job = AMSJob(settings=gas_s, name='gas')
        @add_to_instance(gas_job)
        def prerun(self):
            self.molecule = self.parent.children['dft_preoptimization'].results.get_main_molecule()
        self.children['gas'] = gas_job

        solv_s = Settings()
        solv_s.input.ams.Task = 'SinglePoint'
        solv_s += self.adf_settings(solvation=True)
        solv_job = AMSJob(settings=solv_s, name='solv')

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
    def adf_settings(solvation:bool) -> Settings:
        """
        Returns ADF settings with or without solvation

        If solvation == True, then also include the solvation block.
        """
        
        s = """
        Engine ADF
            Basis
                Type TZP
                Core Small
            End
            XC
                GGA BP86
            End
            Symmetry NOSYM
        """
        if solvation:
            s += """

            Solvation
                Radii
                    Br       2.16
                    C        2.0
                    Cl       2.05
                    F        1.72
                    H        1.3
                    I        2.32
                    N        1.83
                    O        1.72
                    P        2.13
                    S        2.16
                    Si       2.48
                End
                Surf Delley
                Solv name=CRS cav0=0.0 cav1=0.0
                Charged method=CONJ corr
                C-Mat EXACT
                SCF VAR ALL
                CSMRSP
            End

            """
        s += """

            BeckeGrid
                Quality Good
            End
        EndEngine
        """

        s = AMSJob.from_input(s).settings

        return s

    @staticmethod
    def convert_to_coskf(rkf: str, coskf: str):
        """ rkf: absolute path to adf.rkf, coskf: path to write out the resulting .coskf file  """
        f = KFFile(rkf)
        cosmo = f.read_section("COSMO")
        coskf = KFFile(coskf)
        for k,v in cosmo.items():
            coskf.write("COSMO",k,v)



