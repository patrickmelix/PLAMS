import logging
import sys
from typing import Optional
import threading


from scm.plams.core.errors import FileError

__all__ = ["Logger", "LoggerManager"]


class LoggerManager:
    """
    Manages PLAMS logger instances.
    The manager should not be instantiated directly, but loggers accessed through the ``get_logger`` method.
    """

    _loggers = {}

    def __new__(cls, *args, **kwargs):
        raise TypeError("LoggerManager cannot be directly instantiated.")

    @classmethod
    def get_logger(cls, name: str) -> "Logger":
        """
        Get a logger with the specified name.
        If there is an existing logger with the same name this is returned, otherwise the logger is created.
        """
        if name not in cls._loggers:
            cls._loggers[name] = Logger(name)
        return cls._loggers[name]


class Logger:
    """
    Wrapper around default logger, which handles logging to stdout and a text logfile depending on the options in the
    global logging config.
    """

    def __init__(self, name: str):
        """
        Get a logger with given name.
        """
        self._logger = logging.Logger(name)
        self._stdout_handler: Optional[logging.StreamHandler] = None
        self._file_handler: Optional[logging.FileHandler] = None
        self._formatter: Optional[logging.Formatter] = None
        self._lock = threading.Lock()

    def configure_stdout(self, level: int) -> None:
        """
        Configure the logging to stdout.

        Note that for backwards compatibility, level is between 0-7, with 0 indicating no logging and 7 the most verbose logging.
        Note that this is only a PLAMS logging level, it will be mapped to a level between INFO and WARNING for the
        python logger.
        """
        with self._lock:
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
        A path must be unique across all loggers, otherwise a |FileError| is raised.

        For backwards compatibility, level is between 0-7, with 0 indicating no logging and 7 the most verbose logging.
        Note that this is only a PLAMS logging level, it will be mapped to a level between INFO and WARNING for the
        python logger.
        """
        with self._lock:
            # Remove and close existing file handler if present and required
            if self._file_handler is not None and (path is None or path != self._file_handler.baseFilename):
                self._logger.removeHandler(self._file_handler)
                self._file_handler.flush()
                self._file_handler.close()

            # Add new file handler if required
            if path is not None and (self._file_handler is None or path != self._file_handler.baseFilename):

                # Check logfile is not already in use by another logger
                for name, logger in LoggerManager._loggers.items():
                    if logger._file_handler is not None and logger._file_handler.baseFilename == path:
                        raise FileError(f"Logger '{name}' already exists with path '{path}'")

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

        with self._lock:
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
        # Shouldn't really have to take a lock here as logging itself is thread-safe
        # but in the PLAMS log function the logfile can be reconfigured on a per-call basis
        # which could easily lead to dropped logs if this were multi-threaded.
        with self._lock:
            self._logger.log(28 - level, message)
