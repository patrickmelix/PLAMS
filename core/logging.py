import logging
import sys
from typing import Optional

__all__ = ["Logger"]


class Logger:
    """
    Wrapper around default logger, which handles logging to stdout and a text logfile depending on the options in the
    global logging config.
    """

    def __init__(self, name: str):
        """
        Set up logger with given name.
        """
        self._logger = logging.Logger(name)
        self._stdout_handler: Optional[logging.StreamHandler] = None
        self._file_handler: Optional[logging.FileHandler] = None
        self._formatter: Optional[logging.Formatter] = None

    def configure_stdout(self, level: int) -> None:
        """
        Configure the logging to stdout.

        Note that for backwards compatibility, level is between 0-7, with 0 indicating no logging and 7 the most verbose logging.
        Note that this is only a PLAMS logging level, it will be mapped to a level between INFO and WARNING for the
        python logger.
        """
        if self._stdout_handler is None:
            self._stdout_handler = logging.StreamHandler(sys.stdout)
            self._stdout_handler.setFormatter(self._formatter)
            self._logger.addHandler(self._stdout_handler)

        if level != self._stdout_handler.level:
            self._stdout_handler.setLevel(28 - level)

    def configure_logfile(self, path: Optional[str], level: int = 0) -> None:
        """
        Configure the logging to a logfile.

        Path is the file path for the logfile. If set to None this will remove file logging.

        For backwards compatibility, level is between 0-7, with 0 indicating no logging and 7 the most verbose logging.
        Note that this is only a PLAMS logging level, it will be mapped to a level between INFO and WARNING for the
        python logger.
        """
        # Remove and close existing file handler if present and required
        if self._file_handler is not None and (path is None or path != self._file_handler.baseFilename):
            self._logger.removeHandler(self._file_handler)
            self._file_handler.close()

        # Add new file handler if required
        if path is not None and (self._file_handler is None or path != self._file_handler.baseFilename):
            self._file_handler = logging.FileHandler(path)
            self._file_handler.setFormatter(self._formatter)
            self._logger.addHandler(self._file_handler)

        if self._file_handler is not None:
            self._file_handler.setLevel(28 - level)

    def configure_formatter(self, include_date: bool, include_time: bool) -> None:
        """
        Configure the formatting of the logging both to stdout and to the logfile.
        Choose whether to begin the logline with a date and/or time stamp.
        """
        datefmt = None
        if include_date and include_time:
            datefmt = "[%d.%m|%H:%M:%S]"
        elif include_date:
            datefmt = "[%d.%m]"
        elif include_time:
            datefmt = "[%H:%M:%S]"

        if self._formatter is None or datefmt != self._formatter.datefmt:
            fmt = "%(asctime)s %(message)s" if datefmt is not None else None
            self._formatter = logging.Formatter(fmt, datefmt=datefmt)
            if self._stdout_handler is not None:
                self._stdout_handler.setFormatter(self._formatter)
            if self._file_handler is not None:
                self._file_handler.setFormatter(self._formatter)

    def log(self, message: str, level: int) -> None:
        """
        Log a message with the given level of verbosity.
        """
        self._logger.log(28 - level, message)
