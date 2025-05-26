import pytest
from dataclasses import dataclass
from typing import Optional, List

from scm.plams.mol.molecule import Molecule
from scm.plams.interfaces.molecule.rdkit import from_smiles
from scm.plams.interfaces.molecule.packmol import PackMolStructure


class TestPackmolStructure:

    @dataclass
    class TestCase:
        """
        Dataclass holding inputs for PackMolStructure and the expected values for attributes after initialization.
        """

        molecule: Molecule
        n_molecules: Optional[int] = None
        n_atoms: Optional[int] = None
        box_bounds: Optional[List[float]] = None
        density: Optional[float] = None
        fixed: bool = False
        sphere: bool = False
        expected_n_molecules: Optional[int] = None
        expected_box_bounds: Optional[List[float]] = None
        expected_fixed: Optional[bool] = False
        expected_sphere: Optional[bool] = False
        expected_error: Optional[str] = None
        expected_input: Optional[str] = None

    water = from_smiles("O")
    water_in_box = water.copy()
    water_in_box.lattice = [[10, 0, 0], [0, 10, 0], [0, 0, 10]]

    test_cases = [
        # Happy
        TestCase(
            water,
            fixed=True,
            expected_n_molecules=1,
            expected_fixed=True,
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        TestCase(
            water,
            fixed=True,
            n_molecules=1,
            expected_n_molecules=1,
            expected_fixed=True,
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        TestCase(
            water,
            fixed=True,
            box_bounds=[0, 0, 0, 5, 5, 5],
            expected_n_molecules=1,
            expected_fixed=True,
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        TestCase(
            water_in_box,
            fixed=True,
            expected_n_molecules=1,
            expected_fixed=True,
            expected_box_bounds=[0, 0, 0, 10, 10, 10],
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        TestCase(
            water,
            fixed=True,
            sphere=True,
            expected_n_molecules=1,
            expected_fixed=True,
            expected_input="""\
structure pm
  number 1
  fixed 0. 0. 0. 0. 0. 0.
  avoid_overlap yes
end structure
""",
        ),
        TestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            density=0.9,
            expected_n_molecules=30,
            expected_box_bounds=[0, 0, 0, 10, 10, 10],
            expected_input="""\
structure pm
  number 30
  inside box 0.05 0.05 0.05 9.95 9.95 9.95
end structure
""",
        ),
        TestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            n_molecules=20,
            expected_n_molecules=20,
            expected_box_bounds=[0, 0, 0, 10, 10, 10],
            expected_input="""\
structure pm
  number 20
  inside box 0.05 0.05 0.05 9.95 9.95 9.95
end structure
""",
        ),
        TestCase(
            water,
            n_molecules=25,
            density=1.0,
            expected_n_molecules=25,
            expected_box_bounds=[0.0, 0.0, 0.0, 9.077035146666434, 9.077035146666434, 9.077035146666434],
            expected_input="""\
structure pm
  number 25
  inside box 0.05 0.05 0.05 9.027035146666433 9.027035146666433 9.027035146666433
end structure
""",
        ),
        TestCase(
            water,
            n_atoms=75,
            density=1.0,
            expected_n_molecules=25,
            expected_box_bounds=[0.0, 0.0, 0.0, 9.077035146666434, 9.077035146666434, 9.077035146666434],
            expected_input="""\
structure pm
  number 25
  inside box 0.05 0.05 0.05 9.027035146666433 9.027035146666433 9.027035146666433
end structure
""",
        ),
        TestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            n_atoms=60,
            sphere=True,
            expected_n_molecules=20,
            expected_box_bounds=[0, 0, 0, 10, 10, 10],
            expected_sphere=True,
            expected_input="""\
structure pm
  number 20
  inside sphere 0. 0. 0. 6.203504908994001
end structure
""",
        ),
        # Unhappy
        TestCase(water, fixed=True, n_molecules=2, expected_error="n_molecules must be 1"),
        TestCase(water, fixed=True, density=42, expected_error="density cannot be set"),
        TestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            density=0.9,
            n_molecules=10,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        TestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            density=0.9,
            n_atoms=30,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        TestCase(
            water,
            n_molecules=30,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        TestCase(
            water,
            box_bounds=[0, 0, 0, 10, 10, 10],
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        TestCase(
            water,
            density=1,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
        TestCase(
            water,
            n_atoms=10,
            expected_error="exactly two of box_bounds, density and n_molecules/n_atoms must be set",
        ),
    ]

    def get_structure_from_test_case(self, test_case):
        return PackMolStructure(
            molecule=test_case.molecule,
            n_molecules=test_case.n_molecules,
            n_atoms=test_case.n_atoms,
            box_bounds=test_case.box_bounds,
            density=test_case.density,
            fixed=test_case.fixed,
            sphere=test_case.sphere,
        )

    @pytest.mark.parametrize("test_case", test_cases)
    def test_structure_as_expected_given_inputs(self, test_case):

        if test_case.expected_error:
            with pytest.raises(ValueError, match=test_case.expected_error):
                self.get_structure_from_test_case(test_case)
        else:
            structure = self.get_structure_from_test_case(test_case)
            assert structure.n_molecules == test_case.expected_n_molecules
            assert structure.box_bounds == test_case.expected_box_bounds
            assert structure.fixed == test_case.expected_fixed
            assert structure.sphere == test_case.expected_sphere

    @pytest.mark.parametrize("test_case", test_cases)
    def test_get_input_block(self, test_case):
        if test_case.expected_error:
            with pytest.raises(ValueError, match=test_case.expected_error):
                self.get_structure_from_test_case(test_case)
        else:
            structure = self.get_structure_from_test_case(test_case)
            input_block = structure.get_input_block("pm", 0.1)

            assert input_block == test_case.expected_input
