import dill as pickle
import pytest
from unittest.mock import MagicMock
from collections import namedtuple

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
        settings.input.ams.Properties.NormalModes = "True"
        settings.input.DFTB.Model = "GFN1-xTB"
        return settings

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """Properties
  NormalModes True
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

    def test_get_runscript_with_runscript_settings_generates_expected_string(self, job_input):
        # Given job with additional runscript settings
        job_input.settings.runscript.preamble_lines = ["# Start"]
        job_input.settings.runscript.postamble_lines = ["# End"]
        job_input.settings.runscript.nproc = 8
        job_input.settings.runscript.stdout_redirect = True
        job = AMSJob(molecule=job_input.molecule, settings=job_input.settings)

        # When get the runscript
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
        molecule.add_atom("O", coords=[0, 0, 0], unit="A")
        molecule.add_atom("H", coords=[1, 0, 0], unit="A")
        molecule.add_atom("H", coords=[0, 1, 0], unit="A")
        return molecule

    @staticmethod
    def get_expected_input():
        """
        Get expected input file
        """
        return """Properties
  NormalModes True
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

    def test_pickle_dumps_and_loads_job_successfully(self, job_input):
        pytest.skip("Cannot pickle ChemicalSystem")
