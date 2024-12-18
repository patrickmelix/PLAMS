import re
import dill as pickle
import pytest
from unittest.mock import MagicMock, patch
from collections import namedtuple
from io import StringIO
import os
import time
import threading

from scm.plams.interfaces.adfsuite.ams import AMSJob, AMSResults
from scm.plams.core.settings import Settings
from scm.plams.mol.molecule import Atom, Molecule
from scm.plams.unit_tests.test_helpers import skip_if_no_scm_pisa, skip_if_no_scm_libbase


class TestAMSJob:
    """
    Test suite for AMSJob without using PISA / CS for input.
    Sets up a geometry optimization of water.
    """

    @pytest.fixture
    def job_input(self):
        JobInput = namedtuple("JobInput", "molecule settings input")
        molecule = self.get_input_molecule()
        settings = self.get_input_settings()
        expected_input = self.get_expected_input()
        return JobInput(molecule, settings, expected_input)

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        molecule = Molecule()
        molecule.add_atom(Atom(symbol="O", coords=(0, 0, 0)))
        molecule.add_atom(Atom(symbol="H", coords=(1, 0, 0)))
        molecule.add_atom(Atom(symbol="H", coords=(0, 1, 0)))
        return molecule

    @staticmethod
    def get_input_settings():
        """
        Get instance of the Settings class passed to the AMSJob
        """
        settings = Settings()
        settings.input.ams.Task = "GeometryOptimization"
        settings.input.ams.Properties.NormalModes = "Yes"
        settings.input.DFTB.Model = "GFN1-xTB"
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """Properties
  NormalModes Yes
End

System
  Atoms
              O       0.0000000000       0.0000000000       0.0000000000
              H       1.0000000000       0.0000000000       0.0000000000
              H       0.0000000000       1.0000000000       0.0000000000
  End
End

Task GeometryOptimization

Engine DFTB
  Model GFN1-xTB
EndEngine

"""

    def test_init_deep_copies_molecule(self, job_input):
        # Given job with molecule
        job = AMSJob(molecule=job_input.molecule)

        # When get molecule from job
        # Then job molecule is a deep copy
        if job_input.molecule is not None:
            assert job.molecule is not job_input.molecule
            assert job.molecule.atoms is not job_input.molecule.atoms
        else:
            assert job.molecule is None

    def test_pickle_dumps_and_loads_job_successfully(self, job_input):
        # Given job with molecule and settings
        job = AMSJob(molecule=job_input.molecule, settings=job_input.settings)

        # When round trip via pickling
        pickle_bytes = pickle.dumps(job)
        job2 = pickle.loads(pickle_bytes)

        # Then job is still of correct type
        assert isinstance(job2, AMSJob)

    def test_get_input_generates_expected_input_string(self, job_input):
        # Given job with molecule and settings
        job = AMSJob(molecule=job_input.molecule, settings=job_input.settings)

        # When get the job input for the input file
        # Then the input matches the expected input
        assert job.get_input() == job_input.input

    def test_get_runscript_generates_expected_string(self, job_input):
        # Given job
        job = AMSJob(molecule=job_input.molecule, settings=job_input.settings)

        # When get the runscript
        runscript = job.get_runscript()

        # Then standard runscript returned
        assert (
            runscript
            == """unset AMS_SWITCH_LOGFILE_AND_STDOUT
unset SCM_LOGFILE
AMS_JOBNAME="plamsjob" AMS_RESULTSDIR=. $AMSBIN/ams --input="plamsjob.in" < /dev/null

"""
        )

    def test_get_runscript_with_runscript_settings_generates_expected_string(self, job_input, config):
        # Given job with additional runscript settings
        job_input.settings.runscript.preamble_lines = ["# Start"]
        job_input.settings.runscript.postamble_lines = ["# End"]
        job_input.settings.runscript.nproc = 8
        job_input.settings.runscript.stdout_redirect = True
        job = AMSJob(molecule=job_input.molecule, settings=job_input.settings)

        # When get the runscript
        with patch("scm.plams.interfaces.adfsuite.ams.config", config):
            config.slurm = None  # Remove any specific settings when running under slurm
            runscript = job.get_runscript()

        # Then runscript with additional lines returned
        assert (
            runscript
            == """unset AMS_SWITCH_LOGFILE_AND_STDOUT
unset SCM_LOGFILE
# Start
AMS_JOBNAME="plamsjob" AMS_RESULTSDIR=. $AMSBIN/ams -n 8 --input="plamsjob.in" < /dev/null >"plamsjob.out"

# End

"""
        )

    @pytest.mark.parametrize(
        "status,expected",
        [
            ["NORMAL TERMINATION", True],
            ["NORMAL TERMINATION with warnings", True],
            ["NORMAL TERMINATION with errors", False],
            ["Input error", False],
            [None, False],
        ],
    )
    def test_check_returns_true_for_normal_termination_with_no_errors_otherwise_false(self, status, expected):
        # Given job with results of certain status
        job = AMSJob()
        job.results = MagicMock(spec=AMSResults)
        job.results.readrkf.return_value = status

        # When check the job
        # Then job check is ok only for normal termination with no errors
        assert job.check() == expected

    @pytest.mark.parametrize(
        "status,expected", [["NORMAL TERMINATION", None], ["NORMAL TERMINATION with errors", "something bad"]]
    )
    def test_get_errormsg_returns_message_from_logfile_on_error_otherwise_none(self, status, expected):
        # Given job with results of certain status
        job = AMSJob()
        job.results = MagicMock(spec=AMSResults)
        job.results.readrkf.return_value = status
        job.results.grep_file.return_value = ["ERROR: something bad"]

        # When get the error message
        # Then the message is none only for normal termination with no errors
        assert job.get_errormsg() == expected


class TestAMSJobWithPisa(TestAMSJob):
    """
    Test suite for AMSJob using PISA for settings input.
    Sets up a geometry optimization of water.
    """

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        skip_if_no_scm_pisa()
        from scm.input_classes.drivers import AMS
        from scm.input_classes.engines import DFTB

        settings = Settings()
        driver = AMS()
        driver.Task = "GeometryOptimization"
        driver.Properties.NormalModes = "True"
        driver.Engine = DFTB()
        driver.Engine.Model = "GFN1-xTB"
        settings.input = driver
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """Properties
  NormalModes True
End
Task GeometryOptimization

Engine DFTB
  Model GFN1-xTB
EndEngine

System
  Atoms
              O       0.0000000000       0.0000000000       0.0000000000
              H       1.0000000000       0.0000000000       0.0000000000
              H       0.0000000000       1.0000000000       0.0000000000
  End
End
"""


class TestAMSJobWithChemicalSystem(TestAMSJob):
    """
    Test suite for AMSJob using ChemicalSystem for molecule input
    """

    @staticmethod
    def get_input_molecule():
        """
        Instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

        molecule = ChemicalSystem()
        molecule.add_atom("O", coords=[0, 0, 0])
        molecule.add_atom("H", coords=[1, 0, 0])
        molecule.add_atom("H", coords=[0, 1, 0])
        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """Properties
  NormalModes Yes
End

System
  Atoms
     O    0.0000000000000000  0.0000000000000000  0.0000000000000000
     H    1.0000000000000000  0.0000000000000000  0.0000000000000000
     H    0.0000000000000000  1.0000000000000000  0.0000000000000000
  End
End

Task GeometryOptimization

Engine DFTB
  Model GFN1-xTB
EndEngine

"""


class TestAMSJobWithMultipleMolecules(TestAMSJob):
    """
    Test suite for AMSJob using multiple molecules.
    Sets up a NEB calculation for the isomerisation of HCN.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        main_molecule = Molecule()
        main_molecule.add_atom(Atom(symbol="C", coords=(0, 0, 0)))
        main_molecule.add_atom(Atom(symbol="N", coords=(1.18, 0, 0)))
        main_molecule.add_atom(Atom(symbol="H", coords=(2.196, 0, 0)))
        final_molecule = main_molecule.copy()
        final_molecule.atoms[1].x = 1.163
        final_molecule.atoms[2].x = -1.078

        molecule = {"": main_molecule, "final": final_molecule}

        return molecule

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        settings = Settings()
        settings.input.ams.Task = "NEB"
        settings.input.ams.NEB.Images = 9
        settings.input.ams.NEB.Iterations = 100
        settings.input.DFTB.Model = "DFTB3"
        settings.input.DFTB.ResourcesDir = "DFTB.org/3ob-3-1"
        settings.input.DFTB.DispersionCorrection = "D3-BJ"
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """NEB
  Images 9
  Iterations 100
End

System
  Atoms
              C       0.0000000000       0.0000000000       0.0000000000
              N       1.1800000000       0.0000000000       0.0000000000
              H       2.1960000000       0.0000000000       0.0000000000
  End
End
System final
  Atoms
              C       0.0000000000       0.0000000000       0.0000000000
              N       1.1630000000       0.0000000000       0.0000000000
              H      -1.0780000000       0.0000000000       0.0000000000
  End
End

Task NEB

Engine DFTB
  DispersionCorrection D3-BJ
  Model DFTB3
  ResourcesDir DFTB.org/3ob-3-1
EndEngine

"""

    def test_init_deep_copies_molecule(self, job_input):
        # Given job with molecule
        job = AMSJob(molecule=job_input.molecule)

        # When get molecule from job
        # Then job molecule is a deep copy
        assert job.molecule is not job_input.molecule
        for name, mol in job.molecule.items():
            assert mol is not job_input.molecule[name]
            assert mol.atoms is not job_input.molecule[name].atoms


class TestAMSJobWithMultipleMoleculesAndPisa(TestAMSJobWithMultipleMolecules):
    """
    Test suite for AMSJob using multiple molecules and PISA.
    Sets up a NEB calculation for the isomerisation of HCN.
    """

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        from scm.input_classes.drivers import AMS
        from scm.input_classes.engines import DFTB

        settings = Settings()
        driver = AMS()
        driver.Task = "NEB"
        driver.NEB.Images = 9
        driver.NEB.Iterations = 100
        driver.Engine = DFTB()
        driver.Engine.Model = "DFTB3"
        driver.Engine.ResourcesDir = "DFTB.org/3ob-3-1"
        driver.Engine.DispersionCorrection = "D3-BJ"
        settings.input = driver
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """NEB
  Images 9
  Iterations 100
End
Task NEB

Engine DFTB
  DispersionCorrection D3-BJ
  Model DFTB3
  ResourcesDir DFTB.org/3ob-3-1
EndEngine

System
  Atoms
              C       0.0000000000       0.0000000000       0.0000000000
              N       1.1800000000       0.0000000000       0.0000000000
              H       2.1960000000       0.0000000000       0.0000000000
  End
End

System final
  Atoms
              C       0.0000000000       0.0000000000       0.0000000000
              N       1.1630000000       0.0000000000       0.0000000000
              H      -1.0780000000       0.0000000000       0.0000000000
  End
End
"""


class TestAMSJobWithMultipleChemicalSystems(TestAMSJobWithMultipleMolecules):
    """
    Test suite for AMSJob using multiple Chemical Systems.
    Sets up a NEB calculation for the isomerisation of HCN.
    """

    @staticmethod
    def get_input_molecule():
        """
        Instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem

        main_molecule = ChemicalSystem()
        main_molecule.add_atom("C", coords=(0, 0, 0))
        main_molecule.add_atom("N", coords=(1, 0, 0))
        main_molecule.add_atom("H", coords=(2, 0, 0))
        final_molecule = main_molecule.copy()
        final_molecule.atoms[2].coords[0] = -1
        molecule = {"": main_molecule, "final": final_molecule}

        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """NEB
  Images 9
  Iterations 100
End

System
  Atoms
     C    0.0000000000000000  0.0000000000000000  0.0000000000000000
     N    1.0000000000000000  0.0000000000000000  0.0000000000000000
     H    2.0000000000000000  0.0000000000000000  0.0000000000000000
  End
End
System final
  Atoms
     C    0.0000000000000000  0.0000000000000000  0.0000000000000000
     N    1.0000000000000000  0.0000000000000000  0.0000000000000000
     H   -1.0000000000000000  0.0000000000000000  0.0000000000000000
  End
End

Task NEB

Engine DFTB
  DispersionCorrection D3-BJ
  Model DFTB3
  ResourcesDir DFTB.org/3ob-3-1
EndEngine

"""


class TestAMSJobWithChainOfMolecules(TestAMSJob):
    """
    Test suite for AMSJob with a chain of water molecules.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        mol = TestAMSJob.get_input_molecule()
        mol.lattice = [[3, 0, 0]]
        return mol.supercell(4)

    @staticmethod
    def get_input_settings():
        """
        Instance of the Settings class passed to the AMSJob
        """
        settings = Settings()
        settings.input.Mopac.SCF.ConvergenceThreshold = 1.0e-8
        settings.input.Mopac.model = "pm6"
        settings.input.AMS.Task = "SinglePoint"
        settings.input.AMS.Properties.Gradients = "Yes"
        settings.input.AMS.NumericalDifferentiation.Parallel.nCoresPerGroup = 1
        settings.input.AMS.NumericalDifferentiation.NuclearStepSize = 0.0001
        settings.input.AMS.EngineDebugging.IgnoreGradientsRequest = "No"
        settings.input.AMS.System.ElectrostaticEmbedding.ElectricField = "0.0 0.0 0.0"
        settings.input.AMS.Task = "SinglePoint"
        settings.input.AMS.Properties.Gradients = "Yes"
        settings.input.AMS.NumericalDifferentiation.Parallel.nCoresPerGroup = 1
        settings.input.AMS.NumericalDifferentiation.NuclearStepSize = 0.0001
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """EngineDebugging
  IgnoreGradientsRequest No
End

NumericalDifferentiation
  NuclearStepSize 0.0001
  Parallel
    nCoresPerGroup 1
  End
End

Properties
  Gradients Yes
End

System
  Atoms
              O       0.0000000000       0.0000000000       0.0000000000
              H       1.0000000000       0.0000000000       0.0000000000
              H       0.0000000000       1.0000000000       0.0000000000
              O       3.0000000000       0.0000000000       0.0000000000
              H       4.0000000000       0.0000000000       0.0000000000
              H       3.0000000000       1.0000000000       0.0000000000
              O       6.0000000000       0.0000000000       0.0000000000
              H       7.0000000000       0.0000000000       0.0000000000
              H       6.0000000000       1.0000000000       0.0000000000
              O       9.0000000000       0.0000000000       0.0000000000
              H      10.0000000000       0.0000000000       0.0000000000
              H       9.0000000000       1.0000000000       0.0000000000
  End
  ElectrostaticEmbedding
    ElectricField 0.0 0.0 0.0
  End
  Lattice
        12.0000000000     0.0000000000     0.0000000000
  End
End

Task SinglePoint

Engine Mopac
  SCF
    ConvergenceThreshold 1e-08
  End
  model pm6
EndEngine

"""


class TestAMSJobWithSystemBlockSettings(TestAMSJob):
    """
    Test suite for AMSJob with system block overrides/settings in the settings object.
    Sets up MD of Lennard-Jones system.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        molecule = Molecule()
        molecule.add_atom(Atom(symbol="Ar", coords=(0, 0, 0)))
        molecule.add_atom(Atom(symbol="Ar", coords=(1.605, 0.9266471820493496, 2.605)))
        molecule.lattice = [[3.21, 0.0, 0.0], [1.605, 2.779941546148048, 0.0], [0.0, 0.0, 5.21]]
        molecule.properties.charge = 42  # value to be overridden
        return molecule

    @staticmethod
    def get_input_settings():
        """
        Get instance of the Settings class passed to the AMSJob
        """
        settings = Settings()
        settings.input.ams.Task = "MolecularDynamics"
        settings.input.ams.MolecularDynamics.nSteps = 200
        settings.input.ams.MolecularDynamics.TimeStep = 5.0
        settings.input.ams.MolecularDynamics.Thermostat.Type = "NHC"
        settings.input.ams.MolecularDynamics.Thermostat.Temperature = 298.15
        settings.input.ams.MolecularDynamics.Thermostat.Tau = 100
        settings.input.ams.MolecularDynamics.Trajectory.SamplingFreq = 20
        settings.input.ams.MolecularDynamics.InitialVelocities.Type = "random"
        settings.input.ams.MolecularDynamics.InitialVelocities.Temperature = 200
        settings.input.ams.System.SuperCell = "4 4 4"
        settings.input.ams.System.PerturbCoordinates = 0.1
        settings.input.ams.System.Charge = 0

        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
MolecularDynamics
  InitialVelocities
    Temperature 200
    Type random
  End
  Thermostat
    Tau 100
    Temperature 298.15
    Type NHC
  End
  TimeStep 5.0
  Trajectory
    SamplingFreq 20
  End
  nSteps 200
End

System
  Atoms
             Ar       0.0000000000       0.0000000000       0.0000000000
             Ar       1.6050000000       0.9266471820       2.6050000000
  End
  Charge 0
  Lattice
         3.2100000000     0.0000000000     0.0000000000
         1.6050000000     2.7799415461     0.0000000000
         0.0000000000     0.0000000000     5.2100000000
  End
  PerturbCoordinates 0.1
  SuperCell 4 4 4
End

Task MolecularDynamics

"""


class TestAMSJobWithSystemBlockSettingsAndPisa(TestAMSJobWithSystemBlockSettings):
    """
    Test suite for AMSJob with system block overrides/settings in the settings object and PISA.
    Sets up MD of Lennard-Jones system.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        return None

    @staticmethod
    def get_input_settings():
        """
        Get instance of the Settings class passed to the AMSJob
        """
        from scm.input_classes.drivers import AMS

        settings = Settings()
        driver = AMS()
        driver.Task = "MolecularDynamics"
        driver.MolecularDynamics.NSteps = 200
        driver.MolecularDynamics.TimeStep = 5.0
        driver.MolecularDynamics.Thermostat.Type = "NHC"
        driver.MolecularDynamics.Thermostat.Temperature = [298.15]
        driver.MolecularDynamics.Thermostat.Tau = 100
        driver.MolecularDynamics.Trajectory.SamplingFreq = 20
        driver.MolecularDynamics.InitialVelocities.Type = "random"
        driver.MolecularDynamics.InitialVelocities.Temperature = 200
        driver.System.Atoms = [
            "Ar 0.0000000000       0.0000000000       0.0000000000",
            "Ar 1.6050000000       0.9266471820       2.6050000000",
        ]
        driver.System.Lattice = [
            "3.2100000000     0.0000000000     0.0000000000",
            "1.6050000000     2.7799415461     0.0000000000",
            "0.0000000000     0.0000000000     5.2100000000",
        ]
        driver.System.SuperCell = [4, 4, 4]
        driver.System.PerturbCoordinates = 0.1
        driver.System.Charge = 0
        settings.input = driver

        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
MolecularDynamics
  InitialVelocities
    Temperature 200.0
    Type Random
  End
  NSteps 200
  Thermostat
    Tau 100.0
    Temperature 298.15
    Type NHC
  End
  TimeStep 5.0
  Trajectory
    SamplingFreq 20
  End
End
System
  Atoms
    Ar 0.0000000000       0.0000000000       0.0000000000
    Ar 1.6050000000       0.9266471820       2.6050000000
  End
  Charge 0.0
  Lattice
    3.2100000000     0.0000000000     0.0000000000
    1.6050000000     2.7799415461     0.0000000000
    0.0000000000     0.0000000000     5.2100000000
  End
  PerturbCoordinates 0.1
  SuperCell 4 4 4
End
Task MolecularDynamics
"""


class TestAMSJobWithSystemBlockSettingsAndChemicalSystem(TestAMSJobWithSystemBlockSettings):
    """
    Test suite for AMSJob with system block overrides/settings in the settings object and a chemical system.
    Sets up MD of Lennard-Jones system.
    """

    @staticmethod
    def get_input_molecule():
        """
        Get instance of the Molecule class passed to the AMSJob
        """
        skip_if_no_scm_libbase()
        from scm.libbase import UnifiedChemicalSystem as ChemicalSystem, UnifiedLattice as Lattice

        molecule = ChemicalSystem()
        molecule.add_atom("Ar", coords=(0, 0, 0))
        molecule.add_atom("Ar", coords=(1.605, 0.9266471820493496, 2.605))
        molecule.lattice = Lattice([[3.21, 0.0, 0.0], [1.605, 2.779941546148048, 0.0], [0.0, 0.0, 5.21]])
        molecule.charge = 42  # value to be overridden
        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """\
MolecularDynamics
  InitialVelocities
    Temperature 200
    Type random
  End
  Thermostat
    Tau 100
    Temperature 298.15
    Type NHC
  End
  TimeStep 5.0
  Trajectory
    SamplingFreq 20
  End
  nSteps 200
End

System
  Atoms
     Ar   0.0000000000000000  0.0000000000000000  0.0000000000000000
     Ar   1.6050000000000000  0.9266471820493496  2.6050000000000000
  End
  Charge 0
  Lattice
     3.2100000000000000   0.0000000000000000   0.0000000000000000
     1.6050000000000000   2.7799415461480477   0.0000000000000000
     0.0000000000000000   0.0000000000000000   5.2100000000000000
  End
  PerturbCoordinates 0.1
  SuperCell 4 4 4
End

Task MolecularDynamics

"""


class TestAMSJobRun:

    def test_run_with_watch_forwards_ams_logs_to_stdout(self, config):
        # Patch the config and the stdout
        from scm.plams.core.logging import get_logger

        with patch("scm.plams.interfaces.adfsuite.ams.config", config), patch(
            "sys.stdout", new_callable=StringIO
        ) as mock_stdout:
            logger = get_logger("test_run_with_watch_forwards_ams_logs_to_stdout")
            with patch("scm.plams.core.functions._logger", logger):
                config.log.date = False
                config.log.time = False

                # Given a dummy job
                job = AMSJob()

                # Which writes logs to logfile periodically on background thread
                logfile = os.path.join(config.default_jobmanager.workdir, job.name, "ams.log")

                def write_logs():
                    for i in range(10):
                        time.sleep(0.01)
                        try:
                            with open(logfile, "a") as f:
                                f.write(f"<Oct21-2024> <09:34:54>  line {i}\n")
                                f.flush()
                        except FileNotFoundError:
                            pass

                background_thread = threading.Thread(target=write_logs)

                def get_runscript() -> str:
                    background_thread.start()
                    return "sleep 1"

                job.get_runscript = get_runscript

                # When run job and watching output
                job.run(watch=True)

                # Then ams logs are also forwarded to the standard output
                stdout = mock_stdout.getvalue().replace("\r\n", "\n").split("\n")
                postrun_lines = [l for l in stdout if re.fullmatch("plamsjob: line \d", l)]
                status_lines = [l for l in stdout if re.fullmatch("JOB plamsjob (STARTED|RUNNING|FINISHED|FAILED)", l)]

                assert len(postrun_lines) == 10
                assert len(status_lines) == 4

                logger.close()

    def test_get_errormsg_populates_correctly_for_different_scenarios(self):
        from scm.plams.core.errors import FileError

        # Invalid license
        job = AMSJob()
        job._error_msg = None
        results = MagicMock(spec=AMSResults)
        results.readrkf.side_effect = FileError()
        results.grep_file.side_effect = FileError()
        results.get_output_chunk.return_value = [
            "LICENSE INVALID",
            "---------------",
            "Your license does not include module AMS version 2024.206 on this machine.",
            "Module AMS",
            "Version 2024.206",
            "Machine: Foo",
            "License file: ./license.txt",
        ]
        job.results = results

        assert (
            job.get_errormsg()
            == """LICENSE INVALID
---------------
Your license does not include module AMS version 2024.206 on this machine.
Module AMS
Version 2024.206
Machine: Foo
License file: ./license.txt"""
        )

        # Invalid input
        job._error_msg = None
        results.grep_file.side_effect = None
        results.grep_file.return_value = [
            '<Dec05-2024> <12:03:49>  ERROR: Input error: value "foo" found in line 13 for multiple choice key "Task" is not an allowed choice'
        ]
        assert (
            job.get_errormsg()
            == 'Input error: value "foo" found in line 13 for multiple choice key "Task" is not an allowed choice'
        )

        # Error in calculation
        job._error_msg = None
        results.readrkf.return_value = "NORMAL TERMINATION with errors"
        results.readrkf.side_effect = None
        results.grep_file.return_value = [
            "<Dec05-2024> <13:44:55>  ERROR: Geometry optimization failed! (Did not converge.)"
        ]
        assert job.get_errormsg() == "Geometry optimization failed! (Did not converge.)"

        # Error in prerun
        job._error_msg = "RuntimeError: something went wrong"
        assert job.get_errormsg() == "RuntimeError: something went wrong"
