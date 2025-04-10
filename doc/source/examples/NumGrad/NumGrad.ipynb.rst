Worked Example
--------------

Initialization
~~~~~~~~~~~~~~

.. code:: ipython3

   import numpy as np
   import multiprocessing
   from scm.plams.recipes.numgrad import NumGradJob
   from scm.plams import config, JobRunner, AMSJob, Settings, AMSResults, Units, init

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/NumGrad/plams_workdir.003

Parallelization
~~~~~~~~~~~~~~~

.. code:: ipython3

   config.default_jobrunner = JobRunner(parallel=True, maxjobs=multiprocessing.cpu_count())

Define molecule
~~~~~~~~~~~~~~~

.. code:: ipython3

   mol = AMSJob.from_input(
       """
   System
     Atoms
       C       0.000000000000       0.138569980000       0.355570700000
       O       0.000000000000       0.187935770000      -1.074466460000
       H       0.882876920000      -0.383123830000       0.697839450000
       H      -0.882876940000      -0.383123830000       0.697839450000
       H       0.000000000000       1.145042790000       0.750208830000
     End
   End
   """
   ).molecule[""]

Setup and run jobs
~~~~~~~~~~~~~~~~~~

.. code:: ipython3

   s = Settings()
   s.input.AMS.Task = "SinglePoint"
   s.input.ForceField.Type = "UFF"
   s.runscript.nproc = 1

   j = NumGradJob(name="numgrad", molecule=mol, settings=s, jobtype=AMSJob)
   r = j.run()
   r.wait()

   s_analytical = s.copy()
   s_analytical.input.AMS.Properties.Gradients = "Yes"
   analytical_job = AMSJob(name="analytical", molecule=mol, settings=s_analytical)
   analytical_job.run()

::

   [05.03|16:16:05] JOB numgrad STARTED
   [05.03|16:16:05] Waiting for job numgrad to finish
   [05.03|16:16:05] JOB numgrad RUNNING
   [05.03|16:16:05] JOB numgrad/1_x_-1 STARTED
   [05.03|16:16:05] JOB numgrad/1_x_1 STARTED
   [05.03|16:16:05] JOB numgrad/1_y_-1 STARTED
   [05.03|16:16:05] JOB numgrad/1_y_1 STARTED
   [05.03|16:16:05] JOB numgrad/1_z_-1 STARTED
   [05.03|16:16:05] JOB numgrad/1_z_1 STARTED
   [05.03|16:16:05] JOB numgrad/2_x_-1 STARTED
   [05.03|16:16:05] JOB numgrad/2_x_1 STARTED
   ... (PLAMS log lines truncated) ...




   <scm.plams.interfaces.adfsuite.ams.AMSResults at 0x77d85d1da100>



   [05.03|16:16:11] JOB analytical RUNNING
   [05.03|16:16:11] JOB analytical FINISHED
   [05.03|16:16:11] JOB analytical SUCCESSFUL

Print results
~~~~~~~~~~~~~

.. code:: ipython3

   numerical_gradients = np.array(r.get_gradient(AMSResults.get_energy)).reshape(-1, 3)
   analytical_gradients = np.array(analytical_job.results.get_gradients()).reshape(-1, 3) * Units.convert(
       1.0, "bohr^-1", "angstrom^-1"
   )

   print("Numerical Gradients (NumGradJob), hartree/angstrom:")
   print(numerical_gradients)
   print("Analytical Gradients, hartree/angstrom:")
   print(analytical_gradients)
   print("Error Gradients, hartree/angstrom:")
   print(numerical_gradients - analytical_gradients)

::

   Numerical Gradients (NumGradJob), hartree/angstrom:
   [[ 1.60868100e-08 -1.69669118e-02  4.89757712e-01]
    [-1.04939876e-09  1.59476070e-02 -4.59496019e-01]
    [-2.43299529e-02  1.43887838e-02 -9.57760819e-03]
    [ 2.43299384e-02  1.43887756e-02 -9.57760381e-03]
    [-5.15117463e-10 -2.77353210e-02 -1.11275734e-02]]
   Analytical Gradients, hartree/angstrom:
   [[ 1.60879471e-08 -1.70121112e-02  4.89788974e-01]
    [-1.04943735e-09  1.59456436e-02 -4.59495951e-01]
    [-2.43388668e-02  1.44021876e-02 -9.58077129e-03]
    [ 2.43388522e-02  1.44021793e-02 -9.58076691e-03]
    [-5.15150510e-10 -2.77378992e-02 -1.11314841e-02]]
   Error Gradients, hartree/angstrom:
   [[-1.13706585e-12  4.51994875e-05 -3.12614668e-05]
    [ 3.85875469e-14  1.96343429e-06 -6.77586056e-08]
    [ 8.91382542e-06 -1.34037122e-05  3.16309721e-06]
    [-8.91382498e-06 -1.34037118e-05  3.16309724e-06]
    [ 3.30475429e-14  2.57812910e-06  3.91074449e-06]]
