__all__ = [
    "PlamsError",
    "FileError",
    "ResultsError",
    "JobError",
    "PTError",
    "UnitsError",
    "MoleculeError",
    "TrajectoryError",
    "MissingOptionalPackageError",
]


class PlamsError(Exception):
    """General PLAMS error."""


class FileError(PlamsError):
    """File or filesystem related error."""


class ResultsError(PlamsError):
    """|Results| related error."""


class JobError(PlamsError):
    """|Job| related error."""


class PTError(PlamsError):
    """:class:`Periodic table<scm.plams.utils.PeriodicTable>` error."""


class UnitsError(PlamsError):
    """:class:`Units converter<scm.plams.utils.Units>` error."""


class MoleculeError(PlamsError):
    """|Molecule| related error."""


class TrajectoryError(PlamsError):
    """:class:`Trajectory<scm.plams.trajectories.TrajectoryFile>` error."""


class MissingOptionalPackageError(PlamsError):
    """Missing optional package related error."""

    def __init__(self, package_name: str):
        super().__init__(
            f"The optional package '{package_name}' is required for this PLAMS functionality, but is not available. Please install and try again."
        )
