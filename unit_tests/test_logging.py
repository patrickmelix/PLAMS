import re
import uuid
from unittest.mock import patch
from io import StringIO

import pytest

from scm.plams.core.errors import FileError
from scm.plams.core.logging import LoggerManager
from scm.plams.unit_tests.test_helpers import temp_file_path


class TestLoggerManager:

    def test_get_logger_returns_existing_or_creates_new(self):
        name1 = str(uuid.uuid4())
        name2 = str(uuid.uuid4())
        logger1 = LoggerManager.get_logger(name1)
        logger2 = LoggerManager.get_logger(name2)
        logger3 = LoggerManager.get_logger(name1)

        assert logger1 == logger3 != logger2


class TestLogger:

    def test_no_logging_to_stdout_by_default(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            logger = LoggerManager.get_logger(str(uuid.uuid4()))
            logger.log("hello", 1)

            assert mock_stdout.getvalue() == ""

    def test_configure_stdout_writes_to_stdout_up_to_and_including_level(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            logger = LoggerManager.get_logger(str(uuid.uuid4()))
            logger.configure_stdout(3)
            for i in range(10):
                logger.log(f"log line {i}", i)

            assert (
                mock_stdout.getvalue()
                == """log line 0
log line 1
log line 2
log line 3
"""
            )

    def test_configure_logfile_writes_to_file_up_to_and_including_level(self):
        with temp_file_path(suffix=".log") as temp_log_file:
            logger = LoggerManager.get_logger(str(uuid.uuid4()))
            logger.configure_logfile(temp_log_file, 3)
            for i in range(10):
                logger.log(f"log line {i}", i)

            with open(temp_log_file) as tf:
                assert (
                    tf.read()
                    == """log line 0
log line 1
log line 2
log line 3
"""
                )
            logger.configure_logfile(None)

    def test_multiple_loggers_cannot_write_to_same_file(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file:
                logger1 = LoggerManager.get_logger(str(uuid.uuid4()))
                logger2 = LoggerManager.get_logger(str(uuid.uuid4()))
                logger1.configure_stdout(2)
                logger2.configure_stdout(3)
                logger1.configure_logfile(temp_log_file, 2)
                with pytest.raises(FileError):
                    logger2.configure_logfile(temp_log_file, 3)

    def test_multiple_loggers_can_write_to_stdout_but_different_files(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file1, temp_file_path(suffix=".log") as temp_log_file2:
                logger1 = LoggerManager.get_logger(str(uuid.uuid4()))
                logger2 = LoggerManager.get_logger(str(uuid.uuid4()))
                logger1.configure_stdout(2)
                logger2.configure_stdout(3)
                logger1.configure_logfile(temp_log_file1, 1)
                logger2.configure_logfile(temp_log_file2, 2)

                for i in range(5):
                    logger1.log(f"From 1, level {i}", i)
                    logger2.log(f"From 2, level {i}", i)

                assert (
                    mock_stdout.getvalue()
                    == """From 1, level 0
From 2, level 0
From 1, level 1
From 2, level 1
From 1, level 2
From 2, level 2
From 2, level 3
"""
                )

                with open(temp_log_file1) as tf1:
                    assert (
                        tf1.read()
                        == """From 1, level 0
From 1, level 1
"""
                    )

                with open(temp_log_file2) as tf2:
                    assert (
                        tf2.read()
                        == """From 2, level 0
From 2, level 1
From 2, level 2
"""
                    )
            logger1.configure_logfile(None)
            logger2.configure_logfile(None)

    def test_same_logger_can_switch_write_files(self):
        with temp_file_path(suffix=".log") as temp_log_file1, temp_file_path(suffix=".log") as temp_log_file2:
            logger = LoggerManager.get_logger(str(uuid.uuid4()))
            logger.configure_logfile(temp_log_file1, 2)

            for i in range(5):
                logger.log(f"To 1, level {i}", i)

            logger.configure_logfile(temp_log_file2, 1)

            for i in range(5):
                logger.log(f"To 2, level {i}", i)

            logger.configure_logfile(None, 1)

            for i in range(5):
                logger.log(f"To None, level {i}", i)

            logger.configure_logfile(temp_log_file1, 2)
            for i in range(5):
                logger.log(f"To 1 again, level {i}", i)

            with open(temp_log_file1) as tf1:
                assert (
                    tf1.read()
                    == """To 1, level 0
To 1, level 1
To 1, level 2
To 1 again, level 0
To 1 again, level 1
To 1 again, level 2
"""
                )

            with open(temp_log_file2) as tf2:

                assert (
                    tf2.read()
                    == """To 2, level 0
To 2, level 1
"""
                )
            logger.configure_logfile(None)

    def test_configure_formatter_prefixes_date_and_or_time_for_stdout_and_file(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file:
                logger = LoggerManager.get_logger(str(uuid.uuid4()))
                logger.configure_stdout(4)
                logger.configure_logfile(temp_log_file, 4)

                dts = [[tf1, tf2] for tf1 in [True, False] for tf2 in [True, False]]
                for date, time in dts:
                    logger.configure_formatter(date, time)
                    for i in range(3, 8):
                        logger.log(f"d={date}, t={time}, level={i}", i)

                logger.configure_formatter(True, False)

                with open(temp_log_file) as tf:
                    for i, (l1, l2) in enumerate(
                        zip(
                            [l for l in mock_stdout.getvalue().replace("\r\n", "\n").split("\n") if l],
                            [l for l in tf.read().replace("\r\n", "\n").split("\n") if l],
                        )
                    ):
                        date, time = dts[i // 2]
                        pattern = f"d={date}, t={time}, level={i % 2 + 3}"
                        if date and time:
                            pattern = "\[\d{2}\.\d{2}\|\d{2}:\d{2}:\d{2}\] " + pattern
                        elif date:
                            pattern = "\[\d{2}\.\d{2}\] " + pattern
                        elif time:
                            pattern = "\[\d{2}:\d{2}:\d{2}\] " + pattern
                        assert re.fullmatch(pattern, l1) is not None and re.fullmatch(pattern, l2) is not None

                logger.configure_logfile(None)

    def test_configure_order_invariant(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with temp_file_path(suffix=".log") as temp_log_file1, temp_file_path(suffix=".log") as temp_log_file2:
                logger = LoggerManager.get_logger(str(uuid.uuid4()))

                logger.configure_stdout(10)
                logger.configure_formatter(True, True)
                logger.configure_logfile(temp_log_file2)
                logger.configure_formatter(True, False)
                logger.configure_logfile(None, 10)
                logger.configure_stdout(2)
                logger.configure_formatter(False, False)
                logger.configure_logfile(temp_log_file1, 2)

                logger.log("A log line", 1)

                with open(temp_log_file1) as tf1:
                    assert (
                        mock_stdout.getvalue()
                        == tf1.read()
                        == """A log line
"""
                    )

                logger.configure_logfile(None)
