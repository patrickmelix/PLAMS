import pytest
from unittest.mock import patch
from pathlib import Path


@pytest.fixture(autouse=True)
def config():
    """
    Instead of re-using the same global config object, patch with a fresh config settings instance.
    """
    from scm.plams.core.functions import config
    from scm.plams.core.settings import ConfigSettings

    # Set the workdir created by the "real" config (triggered by the loading of the module) to be deleted
    config.erase_workdir = True

    config_settings = ConfigSettings()
    config_settings.init = False

    with patch("scm.plams.core.functions.config", config_settings):
        yield config_settings


@pytest.fixture
def xyz_folder():
    """
    Returns the path to the XYZ folder
    """
    p = Path(__file__).parent.absolute() / "xyz"
    assert p.exists()
    return p


@pytest.fixture
def pdb_folder():
    """
    Returns the path to the PDB folder
    """
    p = Path(__file__).parent.absolute() / "pdb"
    assert p.exists()
    return p


@pytest.fixture
def rkf_folder():
    """
    Returns the path to the RKF folder
    """
    p = Path(__file__).parent.absolute() / "rkf"
    assert p.exists()
    return p
