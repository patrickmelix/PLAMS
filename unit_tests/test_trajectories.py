import pytest
from pathlib import Path

from scm.plams.mol.molecule import Molecule
from scm.plams.trajectories.rkfhistoryfile import RKFHistoryFile


@pytest.fixture
def conformers_rkf(rkf_folder):
    return str(Path(rkf_folder) / "conformers" / "conformers.rkf")


class TestRKFHistoryFile:

    def test_frames(self, conformers_rkf):
        history_file = RKFHistoryFile(conformers_rkf)

        num_frames = history_file.get_length()
        assert num_frames == 12

        input_mol = history_file.get_plamsmol()
        assert input_mol.get_formula() == "C2H7NO"

        for i in range(num_frames):
            mol = Molecule()
            crds, _ = history_file.read_frame(1, mol)
            assert mol.label(3) == input_mol.label(3)
            assert mol.label(4) != input_mol.label(4)
