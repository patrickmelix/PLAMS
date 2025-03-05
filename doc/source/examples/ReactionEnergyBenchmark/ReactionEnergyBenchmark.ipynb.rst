Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

   from scm.plams import read_molecules, Settings, AMSJob, init, config, JobRunner

.. code:: ipython3

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/examples/plams_workdir.011

Summary Text File
~~~~~~~~~~~~~~~~~

Set up a results file, ``summary.txt``, which will hold a table with reaction energies.

.. code:: ipython3

   summary_fname = "summary.txt"
   with open(summary_fname, "w", buffering=1) as summary_file:
       summary_file.write(
           "#Method HC7_1 HC7_2 HC7_3 HC7_4 HC7_5 HC7_6 HC7_7 ISOL6_1 ISOL6_2 ISOL6_3 ISOL6_4 ISOL6_5 (kcal/mol)\n"
       )
       summary_file.write("Truhlar_ref 14.34 25.02 1.90 9.81 14.84 193.99 127.22 9.77 21.76 6.82 33.52 5.30\n")
       summary_file.write("Smith_wB97X_ref 28.77 41.09 1.75 6.26 9.30 238.83 157.65 9.32 20.80 1.03 26.43 0.40\n")

Benchmark Calculation Setup
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The molecules used in the benchmark calculations are loaded, and the settings for each engine created.

.. code:: ipython3

   mol_dict = read_molecules("HC7-11")
   mol_dict.update(read_molecules("ISOL6"))

.. code:: ipython3

   s_reaxff = Settings()
   s_reaxff.input.ReaxFF.ForceField = "CHON-2019.ff"

   s_ani2x = Settings()
   s_ani2x.input.MLPotential.Backend = "TorchANI"
   s_ani2x.input.MLPotential.Model = "ANI-2x"

   s_ani1ccx = Settings()
   s_ani1ccx.input.MLPotential.Backend = "TorchANI"
   s_ani1ccx.input.MLPotential.Model = "ANI-1ccx"

   s_dftb = Settings()
   s_dftb.input.DFTB.Model = "SCC-DFTB"
   s_dftb.input.DFTB.ResourcesDir = "DFTB.org/mio-1-1"

   s_band = Settings()
   s_band.input.BAND.XC.libxc = "wb97x"
   s_band.input.BAND.Basis.Type = "TZP"
   s_band.input.BAND.Basis.Core = "None"

   s_adf = Settings()
   s_adf.input.ADF.Basis.Type = "TZP"
   s_adf.input.ADF.Basis.Core = "None"
   s_adf.input.ADF.XC.libxc = "WB97X"

   s_mp2 = Settings()
   s_mp2.input.ADF.Basis.Type = "TZ2P"
   s_mp2.input.ADF.Basis.Core = "None"
   s_mp2.input.ADF.XC.MP2 = ""

The engines in ``s_engine_dict`` below will be used for the calculation - you can remove the engines that you would not like to run (for example, if you do not have the necessary license).

.. code:: ipython3

   s_engine_dict = {
       "ANI-1ccx": s_ani1ccx,
       "ANI-2x": s_ani2x,
       "DFTB": s_dftb,
       "ReaxFF": s_reaxff,
       "ADF": s_adf,
       "MP2": s_mp2,
       "BAND": s_band,
   }

.. code:: ipython3

   s_ams = Settings()
   s_ams.input.ams.Task = "SinglePoint"
   # s_ams.input.ams.Task = 'GeometryOptimization'
   # s_ams.input.ams.GeometryOptimization.CoordinateType = 'Cartesian'

   jobs = dict()

Running Benchmark Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Benchmark calculations are configured to run in parallel with one core per job.

.. code:: ipython3

   import multiprocessing

   maxjobs = multiprocessing.cpu_count()
   config.default_jobrunner = JobRunner(parallel=True, maxjobs=maxjobs)
   config.job.runscript.nproc = 1

   print(f"Running up to {maxjobs} jobs in parallel")

::

   Running up to 12 jobs in parallel

Note: calculations will take some time to run, around 1-2 hours on a modern laptop.

.. code:: ipython3

   with open(summary_fname, "a", buffering=1) as summary_file:
       for engine_name, s_engine in s_engine_dict.items():
           s = s_ams.copy() + s_engine.copy()
           jobs[engine_name] = dict()

           # call .run() for *all* jobs *before* accessing job.results.get_energy() for *any* job
           for mol_name, mol in mol_dict.items():
               jobs[engine_name][mol_name] = AMSJob(settings=s, molecule=mol, name=engine_name + "_" + mol_name)
               jobs[engine_name][mol_name].run()

       for engine_name in s_engine_dict:
           # for each engine, calculate reaction energies
           E = dict()
           for mol_name, job in jobs[engine_name].items():
               E[mol_name] = job.results.get_energy(unit="kcal/mol")
           deltaE_list = [
               ##### HC7/11 ######
               E["22"] - E["1"],
               E["31"] - E["1"],
               E["octane"] - E["2233tetramethylbutane"],
               5 * E["ethane"] - E["hexane"] - 4 * E["methane"],
               7 * E["ethane"] - E["octane"] - 6 * E["methane"],
               3 * E["ethylene"] + 2 * E["ethyne"] - E["adamantane"],
               3 * E["ethylene"] + 1 * E["ethyne"] - E["bicyclo222octane"],
               ###### start ISOL6 ######
               E["p_3"] - E["e_3"],
               E["p_9"] - E["e_9"],
               E["p_10"] - E["e_10"],
               E["p_13"] - E["e_13"],
               E["p_14"] - E["e_14"],
           ]

           out_str = engine_name
           for deltaE in deltaE_list:
               out_str += " {:.1f}".format(deltaE)
           out_str += "\n"

           print(out_str)
           summary_file.write(out_str)

::

   [05.03|09:22:09] JOB ANI-1ccx_ethylene STARTED
   [05.03|09:22:09] JOB ANI-1ccx_hexane STARTED
   [05.03|09:22:09] JOB ANI-1ccx_methane STARTED
   [05.03|09:22:09] JOB ANI-1ccx_2233tetramethylbutane STARTED
   [05.03|09:22:09] JOB ANI-1ccx_bicyclo222octane STARTED
   [05.03|09:22:09] JOB ANI-1ccx_ethylene RUNNING
   [05.03|09:22:09] JOB ANI-1ccx_31 STARTED
   [05.03|09:22:09] JOB ANI-1ccx_adamantane STARTED
   [05.03|09:22:09] JOB ANI-1ccx_hexane RUNNING
   [05.03|09:22:09] JOB ANI-1ccx_ethyne STARTED
   ... (PLAMS log lines truncated) ...
   [05.03|09:22:09] Waiting for job ANI-1ccx_ethylene to finish
   [05.03|09:22:18] Waiting for job ANI-1ccx_ethyne to finish
   [05.03|09:22:26] Waiting for job ANI-1ccx_22 to finish
   [05.03|09:22:26] Waiting for job ANI-1ccx_octane to finish
   [05.03|09:23:44] Waiting for job ANI-1ccx_e_14 to finish
   ANI-1ccx 15.8 30.7 -0.2 7.7 11.7 196.5 127.9 9.2 21.5 3.8 35.4 6.9

   [05.03|10:20:13] Waiting for job ANI-2x_22 to finish
   ANI-2x 42.4 48.1 -1.9 5.4 8.0 238.4 156.9 7.9 20.8 -1.5 26.7 0.5

   [05.03|10:23:23] Waiting for job DFTB_ethylene to finish
   DFTB 7.1 19.0 0.2 3.6 5.4 223.5 150.8 11.7 22.6 7.9 35.2 8.6

   ReaxFF 34.3 23.0 13.4 2.6 5.5 289.2 168.3 -9.2 24.2 70.7 30.9 89.3

   [05.03|10:23:24] Waiting for job ADF_e_10 to finish
   ADF 20.8 31.4 -1.9 6.7 9.9 214.5 139.7 9.5 20.5 4.7 31.0 4.8

   MP2 20.5 27.1 5.4 9.9 15.0 214.8 141.3 11.0 24.6 8.4 36.8 6.2

   [05.03|10:27:29] Waiting for job BAND_bicyclo222octane to finish
   BAND 19.2 30.2 -2.7 7.3 11.0 218.5 141.3 9.9 18.5 5.3 31.7 4.7

Results
~~~~~~~

Results can be loaded from the summary file and displayed using pandas. The pandas package is installable with ``amspackages``.

.. code:: ipython3

   with open(summary_fname) as summary_file:
       try:
           import pandas as pd
           from IPython.display import display, Markdown

           data = pd.read_csv(summary_file, delim_whitespace=True).iloc[:, :-1]
           display(Markdown(data.to_markdown()))

       except ImportError:
           data = summary_file.read()
           print(data)

== =============== ===== ===== ===== ===== ===== ====== ====== ======= ======= ======= ======= =======
\  #Method         HC7_1 HC7_2 HC7_3 HC7_4 HC7_5 HC7_6  HC7_7  ISOL6_1 ISOL6_2 ISOL6_3 ISOL6_4 ISOL6_5
== =============== ===== ===== ===== ===== ===== ====== ====== ======= ======= ======= ======= =======
0  Truhlar_ref     14.34 25.02 1.9   9.81  14.84 193.99 127.22 9.77    21.76   6.82    33.52   5.3
1  Smith_wB97X_ref 28.77 41.09 1.75  6.26  9.3   238.83 157.65 9.32    20.8    1.03    26.43   0.4
2  ANI-1ccx        15.8  30.7  -0.2  7.7   11.7  196.5  127.9  9.2     21.5    3.8     35.4    6.9
3  ANI-2x          42.4  48.1  -1.9  5.4   8     238.4  156.9  7.9     20.8    -1.5    26.7    0.5
4  DFTB            7.1   19    0.2   3.6   5.4   223.5  150.8  11.7    22.6    7.9     35.2    8.6
5  ReaxFF          34.3  23    13.4  2.6   5.5   289.2  168.3  -9.2    24.2    70.7    30.9    89.3
6  ADF             20.8  31.4  -1.9  6.7   9.9   214.5  139.7  9.5     20.5    4.7     31      4.8
7  MP2             20.5  27.1  5.4   9.9   15    214.8  141.3  11      24.6    8.4     36.8    6.2
8  BAND            19.2  30.2  -2.7  7.3   11    218.5  141.3  9.9     18.5    5.3     31.7    4.7
== =============== ===== ===== ===== ===== ===== ====== ====== ======= ======= ======= ======= =======
