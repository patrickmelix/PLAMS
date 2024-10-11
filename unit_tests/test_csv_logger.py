import os
from scm.plams.core.csv_logger import LoggerCSV


def test_csv_logger():
    csv_logger = LoggerCSV()
    message = {"aaa": 5, "bbb": 33}
    csv_logger.log(message=message)
    os.remove(csv_logger.log_path)
