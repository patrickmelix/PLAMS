import pytest

from scm.plams.interfaces.molecule.ase import toASE, fromASE
from scm.plams.unit_tests.test_helpers import get_mock_find_spec
from scm.plams.core.errors import MissingOptionalPackageError


class TestASE:
    """
    Test suite for conversions to/from PLAMS molecule and ASE.
    """

    def roundtrip_and_assert(self, molecules, from_mol, to_mol):
        for name, orig_mol in molecules.items():
            print(f"Testing roundtrip for molecule '{name}'")
            converted_mol = from_mol(orig_mol)
            final_mol = to_mol(converted_mol)
            assert final_mol.label(4) == orig_mol.label(4)

    def test_to_ase_from_ase_roundtrip(self, plams_mols):
        for name, orig_mol in plams_mols.items():
            print(f"Testing roundtrip for molecule '{name}'")
            converted_mol = toASE(orig_mol)
            final_mol = fromASE(converted_mol)
            assert final_mol.label(4) == orig_mol.label(4)

    def test_to_ase_requires_ase_package(self, plams_mols):
        with get_mock_find_spec("scm.plams.core.functions", "ase"):
            with pytest.raises(MissingOptionalPackageError):
                toASE(plams_mols["water"])
