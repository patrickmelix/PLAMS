import re
from unittest.mock import patch
from io import StringIO
import tempfile

from scm.plams.core.logging import Logger


class TestLogger:

    def test_no_logging_to_stdout_by_default(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            logger = Logger("test")
            logger.log("hello", 1)

            assert mock_stdout.getvalue() == ""

    def test_configure_stdout_writes_to_stdout_up_to_and_including_level(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            logger = Logger("test")
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
        with tempfile.NamedTemporaryFile() as temp_log_file:
            logger = Logger("test")
            logger.configure_logfile(temp_log_file.name, 3)
            for i in range(10):
                logger.log(f"log line {i}", i)

            assert (
                temp_log_file.read().decode()
                == """log line 0
log line 1
log line 2
log line 3
"""
            )

    def test_multiple_loggers_write_to_stdout_but_different_files(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with tempfile.NamedTemporaryFile() as temp_log_file1, tempfile.NamedTemporaryFile() as temp_log_file2:
                logger1 = Logger("test1")
                logger2 = Logger("test2")
                logger1.configure_stdout(2)
                logger2.configure_stdout(3)
                logger1.configure_logfile(temp_log_file1.name, 1)
                logger2.configure_logfile(temp_log_file2.name, 2)

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

                assert (
                    temp_log_file1.read().decode()
                    == """From 1, level 0
From 1, level 1
"""
                )
                assert (
                    temp_log_file2.read().decode()
                    == """From 2, level 0
From 2, level 1
From 2, level 2
"""
                )

    def test_same_logger_can_switch_write_files(self):
        with tempfile.NamedTemporaryFile() as temp_log_file1, tempfile.NamedTemporaryFile() as temp_log_file2:
            logger = Logger("test")
            logger.configure_logfile(temp_log_file1.name, 2)

            for i in range(5):
                logger.log(f"To 1, level {i}", i)

            logger.configure_logfile(temp_log_file2.name, 1)

            for i in range(5):
                logger.log(f"To 2, level {i}", i)

            logger.configure_logfile(None, 1)

            for i in range(5):
                logger.log(f"To None, level {i}", i)

            assert (
                temp_log_file1.read().decode()
                == """To 1, level 0
To 1, level 1
To 1, level 2
"""
            )
            assert (
                temp_log_file2.read().decode()
                == """To 2, level 0
To 2, level 1
"""
            )

    def test_configure_formatter_prefixes_date_and_or_time_for_stdout_and_file(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with tempfile.NamedTemporaryFile() as temp_log_file:
                logger = Logger("test")
                logger.configure_stdout(4)
                logger.configure_logfile(temp_log_file.name, 4)

                dts = [[tf1, tf2] for tf1 in [True, False] for tf2 in [True, False]]
                for date, time in dts:
                    logger.configure_formatter(date, time)
                    for i in range(3, 8):
                        logger.log(f"d={date}, t={time}, level={i}", i)

                logger.configure_formatter(True, False)

                for i, (l1, l2) in enumerate(
                    zip(
                        [l for l in mock_stdout.getvalue().replace("\r\n", "\n").split("\n") if l],
                        [l for l in temp_log_file.read().decode().replace("\r\n", "\n").split("\n") if l],
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

    def test_configure_order_invariant(self):
        with patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            with tempfile.NamedTemporaryFile() as temp_log_file1, tempfile.NamedTemporaryFile() as temp_log_file2:
                logger = Logger("test")

                logger.configure_stdout(10)
                logger.configure_formatter(True, True)
                logger.configure_logfile(temp_log_file2.name)
                logger.configure_formatter(True, False)
                logger.configure_logfile(None, 10)
                logger.configure_stdout(2)
                logger.configure_formatter(False, False)
                logger.configure_logfile(temp_log_file1.name, 2)

                logger.log("A log line", 1)

                assert (
                    mock_stdout.getvalue()
                    == temp_log_file1.read().decode()
                    == """A log line
"""
                )
