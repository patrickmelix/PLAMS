#!/usr/bin/env amspython
# coding: utf-8

# # Logging

# ## Logging in python

import logging

logger = logging.getLogger(__name__)
logging.basicConfig(filename=None, level=logging.DEBUG)
logger.debug("This message should go to the log file")
logger.info("So should this")
logger.warning("And this, too")
logger.error("And non-ASCII stuff, too, like Øresund and Malmö")


logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR


# ## Logging in plams

# ### Default logger

import scm.plams as plams

plams.log("asda", level=0)


# ### Extra Loggers

from scm.plams.core.errors import FileError
from scm.plams.core.logging import get_logger, LogManager
from scm.plams.unit_tests.test_helpers import temp_file_path

log_i = get_logger("my_special_workflow")
log_i.configure(stdout_level=7, logfile_level=0, logfile_path=None, include_date=True, include_time=True)
log_i.log("I am not important", 10)
log_i.log("I am very important", 0)


# ### CVSLogger

from scm.plams.core.errors import FileError
from scm.plams.core.logging import get_logger, LogManager
from scm.plams.unit_tests.test_helpers import temp_file_path
from pathlib import Path


logfile_path = Path("prova.log")
logfile_path.unlink(missing_ok=True)


log_csv = get_logger("csv_format_log", format="csv")
log_csv.configure(stdout_level=5, logfile_level=5, logfile_path=logfile_path, include_date=True, include_time=True)


log_csv.log(message={"p": 25, "v": "asdadf", "SSS": 1000}, level=3)


log_csv.log(message={"p": 25, "v": "asdadf", "SSS": 1000}, level=3)


# ### Logging jobs in csv

from scm.plams.core.job_csv_formatter import JobCSVFormatter
from scm.plams.core.logging import get_logger
from pathlib import Path
import scm.plams as plams

logfile_path = Path("prova.log")
logfile_path.unlink(missing_ok=True)
log_csv_job = get_logger("job_csv_logger", format="csv")
log_csv_job.configure(
    stdout_level=0,
    logfile_level=3,
    logfile_path=logfile_path,
    include_date=True,
    include_time=True,
    csv_formatter_cls=JobCSVFormatter,
    counter_len=3,
)


s = plams.Settings()
s.input.ams.Task = "SinglePoint"
s.input.dftb
job = plams.AMSJob(molecule=plams.from_smiles("C"), settings=s)
log_csv_job.log(job, level=0)


job = plams.AMSJob(molecule=plams.from_smiles("C"), settings=s)
log_csv_job.log(job, level=0)


# ### JobManager integrated job logger

import scm.plams as plams

plams.init()
s = plams.Settings()
s.input.ams.Task = "SinglePoint"
s.input.dftb
job = plams.AMSJob(molecule=plams.from_smiles("C"), settings=s)
job.run()
plams.finish()


plams.config.default_jobmanager.logger_csv.log(job, level=6)

