Worked Example
--------------

Logging in PLAMS
~~~~~~~~~~~~~~~~

PLAMS has built-in logging which aims to simplify tracking the progress
and status of jobs. This consists of progress logging to stdout and a
logfile, and writing job summaries to CSV files. Each of these is
explained below.

Progress Logger
~~~~~~~~~~~~~~~

PLAMS writes job progress to stdout and a plain text logfile. The
location of this logfile is determined by the working directory of the
default job manager.

Users can also write logs to the same locations using the ``log``
function. This takes a ``level`` argument. By convention in PLAMS, the
level should be between 0-7, with 0 being the most and 7 the least
important logging.

The level of logging that is written to stdout and the logfile can be
changed through the ``config.LogSettings``.

.. code:: ipython3

    from scm.plams import Settings, AMSJob, from_smiles, log, config
    
    counter = 0
    
    
    def get_test_job():
        global counter
        s = Settings()
        s.input.ams.Task = "SinglePoint"
        s.input.dftb
        counter += 1
        return AMSJob(name=f"test{counter}", molecule=from_smiles("C"), settings=s)

.. code:: ipython3

    config.log.stdout = 3
    config.log.file = 5
    config.jobmanager.hashing = None  # Force PLAMS to re-run identical test jobs

.. code:: ipython3

    job = get_test_job()
    job.run()
    log("Test job finished", 5)


.. parsed-literal::

    [28.11|09:43:00] JOB test1 STARTED
    [28.11|09:43:00] JOB test1 RUNNING
    [28.11|09:43:01] JOB test1 FINISHED
    [28.11|09:43:01] JOB test1 SUCCESSFUL


.. code:: ipython3

    with open(config.default_jobmanager.logfile, "r") as f:
        print(f.read())


.. parsed-literal::

    [28.11|09:43:00] Starting test1.prerun()
    [28.11|09:43:00] test1.prerun() finished
    [28.11|09:43:00] JOB test1 RUNNING
    [28.11|09:43:00] Executing test1.run
    [28.11|09:43:01] Execution of test1.run finished with returncode 0
    [28.11|09:43:01] JOB test1 FINISHED
    [28.11|09:43:01] Starting test1.postrun()
    [28.11|09:43:01] test1.postrun() finished
    [28.11|09:43:01] JOB test1 SUCCESSFUL
    [28.11|09:43:01] Test job finished
    


Note that the logs from an AMS calculation can also be forwarded to the
progress logs using the ``watch = True`` flag.

.. code:: ipython3

    job = get_test_job()
    job.run(watch=True);


.. parsed-literal::

    [28.11|09:43:01] JOB test2 STARTED
    [28.11|09:43:01] JOB test2 RUNNING
    [28.11|09:43:01] test2: AMS 2024.206  RunTime: Nov28-2024 09:43:01  ShM Nodes: 1  Procs: 6
    [28.11|09:43:02] test2: DFTB: SCC cycle
    [28.11|09:43:02] test2: cyc=  1 err=1.1E+00 method=1 nvec= 1 mix=0.075 e=    0.0000
    [28.11|09:43:02] test2: cyc=  2 err=1.1E+00 method=1 nvec= 1 mix=0.154 e=    0.0000
    [28.11|09:43:02] test2: cyc=  3 err=8.9E-01 method=1 nvec= 2 mix=0.201 e=    0.0000
    [28.11|09:43:02] test2: cyc=  4 err=1.7E-02 method=1 nvec= 3 mix=0.207 e=    0.0000
    [28.11|09:43:02] test2: cyc=  5 err=6.8E-03 method=1 nvec= 4 mix=0.213 e=    0.0000
    [28.11|09:43:02] test2: cyc=  6 err=2.6E-03 method=1 nvec= 5 mix=0.219 e=    0.0000
    [28.11|09:43:02] test2: cyc=  7 err=7.2E-05 method=1 nvec= 6 mix=0.226 e=    0.0000
    [28.11|09:43:02] test2: cyc=  8 err=6.8E-05 method=1 nvec= 1 mix=0.233 e=    0.0000
    [28.11|09:43:02] test2: cyc=  9 err=4.2E-05 method=1 nvec= 2 mix=0.240 e=    0.0000
    [28.11|09:43:02] test2: cyc= 10 err=6.2E-07 method=1 nvec= 3 mix=0.247 e=    0.0000
    [28.11|09:43:02] test2: cyc= 11 err=5.8E-08 method=1 nvec= 3 mix=0.254 e=    0.0000
    [28.11|09:43:02] test2: cyc= 12 err=3.6E-08 method=1 nvec= 4 mix=0.262 e=    0.0000
    [28.11|09:43:02] test2: cyc= 13 err=9.0E-11 method=1 nvec= 4 mix=0.270 e=    0.0000
    [28.11|09:43:02] test2: SCC cycle converged!
    [28.11|09:43:02] test2: NORMAL TERMINATION
    [28.11|09:43:02] JOB test2 FINISHED
    [28.11|09:43:02] JOB test2 SUCCESSFUL


Job Summary Logger
~~~~~~~~~~~~~~~~~~

PLAMS also writes summaries of jobs to a CSV file, the location of which
by default is also determined by the job manager.

.. code:: ipython3

    from scm.plams import MultiJob
    
    
    multi_job = MultiJob(name="parent_job", children=[get_test_job() for _ in range(3)])
    multi_job.children[2].settings.input.ams.Task = "Not a task!"
    multi_job.run();


.. parsed-literal::

    [28.11|09:43:02] JOB parent_job STARTED
    [28.11|09:43:02] JOB parent_job RUNNING
    [28.11|09:43:02] JOB parent_job/test3 STARTED
    [28.11|09:43:02] JOB parent_job/test3 RUNNING
    [28.11|09:43:03] JOB parent_job/test3 FINISHED
    [28.11|09:43:03] JOB parent_job/test3 SUCCESSFUL
    [28.11|09:43:03] JOB parent_job/test4 STARTED
    [28.11|09:43:03] JOB parent_job/test4 RUNNING
    [28.11|09:43:04] JOB parent_job/test4 FINISHED
    [28.11|09:43:04] JOB parent_job/test4 SUCCESSFUL
    [28.11|09:43:04] JOB parent_job/test5 STARTED
    [28.11|09:43:04] JOB parent_job/test5 RUNNING
    [28.11|09:43:11] WARNING: Job test5 finished with nonzero return code
    [28.11|09:43:11] WARNING: Main KF file ams.rkf not present in /path/plams/examples/LoggingExamples/plams_workdir/parent_job/test5
    [28.11|09:43:11] JOB parent_job/test5 CRASHED
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] Obtaining results of test5 failed. Returned value is None
    [28.11|09:43:11] Obtaining results of test5 successful. However, no guarantee that they make sense
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] Obtaining results of test5 failed. Returned value is None
    [28.11|09:43:11] Obtaining results of test5 successful. However, no guarantee that they make sense
    [28.11|09:43:11] Could not read termination status from file None
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] Obtaining results of test5 failed. Returned value is None
    [28.11|09:43:11] Obtaining results of test5 successful. However, no guarantee that they make sense
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] Obtaining results of test5 failed. Returned value is None
    [28.11|09:43:11] Obtaining results of test5 successful. However, no guarantee that they make sense
    [28.11|09:43:11] Could not read termination status from file None
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] Obtaining results of test5 failed. Returned value is None
    [28.11|09:43:11] Obtaining results of test5 successful. However, no guarantee that they make sense
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] WARNING: Trying to obtain results of crashed or failed job test5
    [28.11|09:43:11] Obtaining results of test5 successful. However, no guarantee that they make sense
    [28.11|09:43:11] Obtaining results of test5 successful. However, no guarantee that they make sense
    [28.11|09:43:11] JOB parent_job FINISHED
    [28.11|09:43:11] JOB parent_job FAILED


These CSVs give overall information on the status of all jobs run by a
given job manager.

.. code:: ipython3

    try:
        import pandas as pd
    
        df = pd.read_csv(config.default_jobmanager.job_logger.logfile)
        print(df[["job_name", "job_status", "job_ok", "job_get_errormsg"]])
    except ImportError:
        pass


.. parsed-literal::

         job_name  job_status  job_ok  \
    0       test1  successful    True   
    1       test2  successful    True   
    2       test3  successful    True   
    3       test4  successful    True   
    4       test5     crashed   False   
    5  parent_job      failed   False   
    
                                        job_get_errormsg  
    0                                                NaN  
    1                                                NaN  
    2                                                NaN  
    3                                                NaN  
    4  Input error: value "Not a task!" found in line...  
    5                                                NaN  


Extra Loggers
~~~~~~~~~~~~~

PLAMS also provides access to text and csv loggers for use in your
workflows. These are obtained with the ``get_logger`` method. Loggers
need to be set up with the ``configure`` method, which determines which
file to write to and whether to write to stdout.

.. code:: ipython3

    from scm.plams import get_logger
    
    my_text_logger = get_logger("my_text_logger", "txt")
    my_csv_logger = get_logger("my_csv_logger", "csv")
    
    my_text_logger.configure(logfile_level=3, logfile_path="./my_text_logger.txt", include_date=True, include_time=True)
    my_csv_logger.configure(logfile_level=3, logfile_path="./my_csv_logger.csv")
    
    my_text_logger.log("I want to log this to a separate file!", 3)
    my_csv_logger.log({"molecule": "Acetylacetonate", "short_name": "acac"}, 3)
    my_csv_logger.log({"molecule": "Ethylenediaminetetraaceticacid acid", "short_name": "EDTA"}, 3)
    
    with open(my_text_logger.logfile) as f:
        print(f.read())
    
    with open(my_csv_logger.logfile) as f:
        print(f.read())


.. parsed-literal::

    [28.11|09:43:12] I want to log this to a separate file!
    
    molecule,short_name
    Acetylacetonate,acac
    Ethylenediaminetetraaceticacid acid,EDTA
    

