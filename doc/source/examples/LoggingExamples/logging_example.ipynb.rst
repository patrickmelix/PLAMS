Worked Example
--------------

Logging
=======

Logging in python
-----------------

.. code:: ipython3

    import logging
    
    logger = logging.getLogger(__name__)
    logging.basicConfig(filename=None, level=logging.DEBUG)
    logger.debug("This message should go to the log file")
    logger.info("So should this")
    logger.warning("And this, too")
    logger.error("And non-ASCII stuff, too, like Øresund and Malmö")


.. parsed-literal::

    DEBUG:__main__:This message should go to the log file
    INFO:__main__:So should this
    WARNING:__main__:And this, too
    ERROR:__main__:And non-ASCII stuff, too, like Øresund and Malmö


.. code:: ipython3

    logging.DEBUG, logging.INFO, logging.WARNING, logging.ERROR




.. parsed-literal::

    (10, 20, 30, 40)



Logging in plams
----------------

Default logger
~~~~~~~~~~~~~~

.. code:: ipython3

    import scm.plams as plams
    
    plams.log("asda", level=0)


.. parsed-literal::

    [22.11|16:30:31] asda


Extra Loggers
~~~~~~~~~~~~~

.. code:: ipython3

    from scm.plams.core.errors import FileError
    from scm.plams.core.logging import get_logger, LogManager
    from scm.plams.unit_tests.test_helpers import temp_file_path
    
    log_i = get_logger("my_special_workflow")
    log_i.configure(stdout_level=7, logfile_level=0, logfile_path=None, include_date=True, include_time=True)
    log_i.log("I am not important", 10)
    log_i.log("I am very important", 0)


.. parsed-literal::

    [22.11|16:30:33] I am very important


CVSLogger
~~~~~~~~~

.. code:: ipython3

    from scm.plams.core.errors import FileError
    from scm.plams.core.logging import get_logger, LogManager
    from scm.plams.unit_tests.test_helpers import temp_file_path
    from pathlib import Path

.. code:: ipython3

    logfile_path = Path("prova.log")
    logfile_path.unlink(missing_ok=True)

.. code:: ipython3

    log_csv = get_logger("csv_format_log", format="csv")
    log_csv.configure(stdout_level=5, logfile_level=5, logfile_path=logfile_path, include_date=True, include_time=True)

.. code:: ipython3

    log_csv.log(message={"p": 25, "v": "asdadf", "SSS": 1000}, level=3)


.. parsed-literal::

    asctime,p,v,SSS
    2024-11-22 16:31:58,25,asdadf,1000


.. code:: ipython3

    log_csv.log(message={"p": 25, "v": "asdadf", "SSS": 1000}, level=3)


.. parsed-literal::

    2024-11-22 16:31:59,25,asdadf,1000


Logging jobs in csv
~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

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

.. code:: ipython3

    s = plams.Settings()
    s.input.ams.Task = "SinglePoint"
    s.input.dftb
    job = plams.AMSJob(molecule=plams.from_smiles("C"), settings=s)
    log_csv_job.log(job, level=0)


.. parsed-literal::

    2024-11-22 16:33:54,plamsjob,created,,,,,,,


.. code:: ipython3

    job = plams.AMSJob(molecule=plams.from_smiles("C"), settings=s)
    log_csv_job.log(job, level=0)


.. parsed-literal::

    2024-11-22 16:33:55,plamsjob,created,,,,,,,


JobManager integrated job logger
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    import scm.plams as plams
    
    plams.init()
    s = plams.Settings()
    s.input.ams.Task = "SinglePoint"
    s.input.dftb
    job = plams.AMSJob(molecule=plams.from_smiles("C"), settings=s)
    job.run()
    plams.finish()


.. parsed-literal::

    PLAMS working folder: /path/plams/plams_workdir.002
    [22.11|14:08:53] JOB plamsjob STARTED
    [22.11|14:08:53] JOB plamsjob RUNNING
    [22.11|14:08:55] JOB plamsjob FINISHED
    [22.11|14:08:55] JOB plamsjob SUCCESSFUL
    PLAMS run finished. Goodbye


.. code:: ipython3

    plams.config.default_jobmanager.logger_csv.log(job, level=6)
