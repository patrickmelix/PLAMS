import csv
import logging
import re
import sys
import threading
import time
from abc import ABC
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Literal, Mapping, Optional, Union

from scm.plams.core.errors import FileError

__all__ = ["get_logger", "Logger"]


def get_logger(name: str, format: Optional[Literal["csv"]] = "str") -> "Logger":
    """
    Get a logger with the specified name.
    If there is an existing logger with the same name this is returned, otherwise a new logger is created.
    """
    return LogManager.get_logger(name, format=format)


class LogManager:
    """
    Manages PLAMS logger instances.
    The manager should not be instantiated directly, but loggers accessed through the ``get_logger`` method.
    """

    _loggers: Dict[str, "Logger"] = {}

    def __new__(cls, *args, **kwargs):
        raise TypeError("LoggerManager cannot be directly instantiated.")

    @classmethod
    def get_logger(cls, name: str, format: Optional[Literal["csv", "str"]] = "str") -> "Logger":
        """
        Get a logger with the specified name.
        If there is an existing logger with the same name this is returned, otherwise a new logger is created.
        """
        if name not in cls._loggers:
            if format == "csv":
                logger = CSVLogger(name)
            else:
                logger = StrLogger(name)
            cls._loggers[name] = logger
        return cls._loggers[name]


class Logger(ABC):
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

    def configure(
        self,
        stdout_level: int = 0,
        logfile_level: int = 0,
        logfile_path: Optional[str] = None,
        include_date: bool = False,
        include_time: bool = False,
    ) -> None:
        """
        Configure logging to stdout and the logfile, and its formatting.

        For backwards compatibility, the logging level is between 0-7, with 0 indicating no logging and 7 the most verbose logging.
        Note that this is only a PLAMS logging level, it will be mapped to a level between INFO and WARNING for the
        standard python logger.

        The logfile path must be unique across all loggers, otherwise a |FileError| is raised. If set to None this will disable file logging.

        :param stdout_level: value between 0-7, with 0 indicating no logging and 7 the most verbose logging to stdout
        :param logfile_level: value between 0-7, with 0 indicating no logging and 7 the most verbose logging to logfile
        :param logfile_path: path for the logfile, if set to None this will remove file logging
        :param include_date: whether to include date stamp at the start of a log line
        :param include_time: whether to include time stamp at the start of a log line
        """
        with self._lock:

            # Initialise the stdout handler once
            if self._stdout_handler is None:
                self._stdout_handler = logging.StreamHandler(sys.stdout)
                self._stdout_handler.setFormatter(self._formatter)
                self._logger.addHandler(self._stdout_handler)

            # Update the stdout handler level if required
            if stdout_level != self._stdout_handler.level:
                self._stdout_handler.setLevel(28 - stdout_level)

            # Remove and close existing file handler if present and required
            if self._file_handler is not None and (
                logfile_path is None or logfile_path != self._file_handler.baseFilename
            ):
                self._logger.removeHandler(self._file_handler)
                self._file_handler.flush()
                self._file_handler.close()

            # Add new file handler if required
            if logfile_path is not None and (
                self._file_handler is None or logfile_path != self._file_handler.baseFilename
            ):
                # Check logfile is not already in use by another logger
                # This will cause permission errors on Windows and generally is not a good idea
                for name, logger in LogManager._loggers.items():
                    if logger._file_handler is not None and logger._file_handler.baseFilename == logfile_path:
                        raise FileError(f"Logger '{name}' already exists with logfile path '{logfile_path}'")

                self._file_handler = logging.FileHandler(logfile_path)
                self._file_handler.setFormatter(self._formatter)
                self._logger.addHandler(self._file_handler)

            # Update the logfile handler level if required
            if self._file_handler is not None:
                self._file_handler.setLevel(28 - logfile_level)

            # Configure the formatter if required
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

    def log(self, message: Any, level: int) -> None:
        """
        Log a message with the given level of verbosity.
        """
        pass


class StrLogger(Logger):
    """
    Wrapper around default logger, which handles logging to stdout and a text logfile depending on the options in the
    global logging config.
    """

    def __init__(self, name: str):
        """
        Get a logger with given name.
        """
        super().__init__(name)

    def configure(self, stdout_level=0, logfile_level=0, logfile_path=None, include_date=False, include_time=False):
        return super().configure(stdout_level, logfile_level, logfile_path, include_date, include_time)

    def log(self, message: str, level: int) -> None:
        """
        Log a message with the given level of verbosity.
        """
        # Shouldn't really have to take a lock here as logging itself is thread-safe
        # but in the PLAMS log function the logfile can be configured on a per-call basis
        # which could easily lead to dropped logs if this were multi-threaded.
        with self._lock:
            self._logger.log(28 - level, message)


class CSVFormatter(logging.Formatter):
    def __init__(
        self,
        fmt: Optional[str] = None,
        datefmt: Optional[str] = None,
        style: Union[Literal["%"], Literal["{"], Literal["$"]] = "%",
        validate: bool = True,
        log_level: bool = True,
        log_logger_name: bool = False,
        log_time: bool = True,
        log_lineno: bool = False,
        **csv_args,
    ) -> None:
        """
        Initialize the CSVFormatter.

        :param fmt: Format string for time and other placeholders.
        :param datefmt: Date format string.
        :param style: Formatting style, one of '%', '{', or '$'.
        :param validate: Validate the format string.
        :param default_headers: List of default headers for the CSV. If None, headers are inferred dynamically.
        """
        super().__init__(fmt, datefmt, style, validate)

        self.log_level = log_level
        self.log_logger_name = log_logger_name
        self.log_time = log_time
        self.log_lineno = log_lineno

        self.headers = None
        self.written_headers = False

    def format(self, record: logging.LogRecord) -> str:
        """
        Format the log record into a CSV row.

        :param record: The log record to format.
        :return: A CSV-formatted string.
        """
        # Extract core log fields
        log_record = {}
        if self.log_level:
            log_record["level"] = record.levelname
        if self.log_logger_name:
            log_record["logger_name"] = record.name
        if self.log_time:
            log_record["asctime"] = self.formatTime(record, self.datefmt)
        if self.log_logger_name:
            log_record["filename"] = record.filename
        if self.log_logger_name:
            log_record["lineno"] = record.lineno

        if isinstance(record.msg, dict):
            log_record.update(record.msg)
        else:
            log_record["message"] = record.getMessage()

        if self.headers is None:
            self.headers = list(log_record.keys())

        row = StringIO()
        csv_writer = csv.DictWriter(row, fieldnames=self.headers, lineterminator="\r\n")
        if not self.written_headers:
            csv_writer.writeheader()
            self.written_headers = True

        csv_writer.writerow(log_record)
        return row.getvalue().strip()


class CSVLogger(Logger):
    """
    Logger that logs dictionary messages in CSV format with optional date and time stamps.
    """

    def __init__(self, name: str):
        """
        Get a logger with given name.
        """
        super().__init__(name)

    def configure(
        self,
        stdout_level=0,
        logfile_level=0,
        logfile_path=None,
        include_date=False,
        include_time=False,
        log_level=False,
        log_logger_name=False,
        log_lineno=False,
        csv_formatter_cls=CSVFormatter,
        **csv_args,
    ):
        super().configure(stdout_level, logfile_level, logfile_path, include_date, include_time)

        datefmt = None
        log_time = True
        if include_date and include_time:
            datefmt = "%d.%m %H:%M:%S"
        elif include_date:
            datefmt = "%d.%m"
        elif include_time:
            datefmt = "%H:%M:%S"
        else:
            log_time = False

        if self._formatter is None or datefmt != self._formatter.datefmt:
            self._formatter: Optional[logging.Formatter] = csv_formatter_cls(
                fmt=None,
                datefmt=datefmt,
                log_time=log_time,
                log_level=log_level,
                log_logger_name=log_logger_name,
                log_lineno=log_lineno,
                **csv_args,
            )
            if self._stdout_handler is not None:
                self._stdout_handler.setFormatter(self._formatter)
            if self._file_handler is not None:
                # this is needed to write down the headers as swell
                self._file_handler.setFormatter(
                    csv_formatter_cls(
                        fmt=None,
                        datefmt=datefmt,
                        log_time=log_time,
                        log_level=log_level,
                        log_logger_name=log_logger_name,
                        log_lineno=log_lineno,
                        **csv_args,
                    )
                )

    def log(self, message: Mapping[str, Any], level: int = 2) -> None:
        """
        Log a dictionary message in CSV format with verbosity level.

        Logs are written to the CSV logfile and optionally printed to stdout.

        :param message: The message to log as a dictionary.
        :type message: dict
        :param level: Verbosity level (1=important, 3=normal, 5=verbose, 7=debug), defaults to 7.
        :type level: int
        """
        with self._lock:
            self._logger.log(28 - level, message)
