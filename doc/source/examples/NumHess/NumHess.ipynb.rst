Worked Example
--------------

Initialization
~~~~~~~~~~~~~~

.. code:: ipython3

   from scm.plams import *

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/NumHess/plams_workdir.002

Define molecule
~~~~~~~~~~~~~~~

Normally you would read the molecule from an xyz file

.. code:: ipython3

   #mol = Molecule('methanol.xyz')


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

Setup and run job
~~~~~~~~~~~~~~~~~

Here we use the fast DFTB engine, for ADF it is recommended to disable symmetry

.. code:: ipython3

   s = Settings()
   s.input.ams.task = 'SinglePoint'
   s.input.ams.Properties.Gradients = 'Yes'
   # s.input.adf.basis.type = 'DZP'
   # s.input.adf.symmetry = 'NOSYM'
   # s.input.adf.xc.gga = 'PW91'
   s.input.DFTB.Model = 'GFN1-xTB'
   s.runscript.nproc = 1

   j = NumHessJob(name='test', molecule=mol, settings=s, jobtype=AMSJob,
                  gradient = lambda x: x.get_gradients().reshape(-1))
   r = j.run(jobrunner=JobRunner(parallel=True, maxjobs=8))

::

   [05.03|19:46:22] JOB test STARTED
   [05.03|19:46:22] Renaming job test to test.003
   [05.03|19:46:22] JOB test.003 RUNNING
   [05.03|19:46:22] JOB test.003/1_0_-1 STARTED
   [05.03|19:46:22] JOB test.003/1_0_1 STARTED
   [05.03|19:46:22] JOB test.003/1_1_-1 STARTED
   [05.03|19:46:22] JOB test.003/1_1_1 STARTED
   [05.03|19:46:22] Job 1_0_-1 previously run as 1_0_-1, using old results
   [05.03|19:46:22] JOB test.003/1_2_-1 STARTED
   [05.03|19:46:22] Job 1_0_1 previously run as 1_0_1, using old results
   [05.03|19:46:22] JOB test.003/1_2_1 STARTED
   [05.03|19:46:22] JOB test.003/2_0_-1 STARTED
   ... (PLAMS log lines truncated) ...
   [05.03|19:46:22] Job 1_1_1 previously run as 1_1_1, using old results
   [05.03|19:46:22] Job 1_1_-1 previously run as 1_1_-1, using old results
   [05.03|19:46:22] Job 1_2_-1 previously run as 1_2_-1, using old results
   [05.03|19:46:22] Job 1_2_1 previously run as 1_2_1, using old results
   [05.03|19:46:22] Job 2_0_-1 previously run as 2_0_-1, using old results
   [05.03|19:46:22] Job 2_0_1 previously run as 2_0_1, using old results
   [05.03|19:46:22] Job 2_1_-1 previously run as 2_1_-1, using old results
   [05.03|19:46:22] Job 2_2_-1 previously run as 2_2_-1, using old results
   [05.03|19:46:22] Job 2_1_1 previously run as 2_1_1, using old results
   [05.03|19:46:22] Job 2_2_1 previously run as 2_2_1, using old results
   [05.03|19:46:22] Job 3_0_-1 previously run as 3_0_-1, using old results
   [05.03|19:46:22] Job 3_0_1 previously run as 3_0_1, using old results
   [05.03|19:46:22] Job 3_1_-1 previously run as 3_1_-1, using old results
   [05.03|19:46:22] Job 3_1_1 previously run as 3_1_1, using old results
   [05.03|19:46:22] Job 3_2_-1 previously run as 3_2_-1, using old results
   [05.03|19:46:22] Job 4_0_-1 previously run as 4_0_-1, using old results
   [05.03|19:46:22] Job 3_2_1 previously run as 3_2_1, using old results
   [05.03|19:46:22] Job 4_1_-1 previously run as 4_1_-1, using old results
   [05.03|19:46:22] Job 4_0_1 previously run as 4_0_1, using old results
   [05.03|19:46:22] Job 4_2_-1 previously run as 4_2_-1, using old results
   [05.03|19:46:22] Job 4_1_1 previously run as 4_1_1, using old results
   [05.03|19:46:22] Job 4_2_1 previously run as 4_2_1, using old results
   [05.03|19:46:22] Job 5_0_-1 previously run as 5_0_-1, using old results
   [05.03|19:46:22] Job 5_0_1 previously run as 5_0_1, using old results
   [05.03|19:46:22] Job 5_1_-1 previously run as 5_1_-1, using old results
   [05.03|19:46:22] Job 5_1_1 previously run as 5_1_1, using old results
   [05.03|19:46:22] Job 5_2_-1 previously run as 5_2_-1, using old results
   [05.03|19:46:22] Job 5_2_1 previously run as 5_2_1, using old results

Print results
~~~~~~~~~~~~~

.. code:: ipython3

   print(r.get_hessian(mass_weighted=True))

::

   [[ 2.43775998e-03 -1.24639788e-10  8.18469493e-11 -2.16941229e-05
      4.93482678e-12 -7.23564494e-12 -1.60687481e-02  1.15383294e-02
     -7.21385604e-03 -1.60687464e-02 -1.15383280e-02  7.21385523e-03
      3.47625118e-03  4.94064021e-11 -4.51157758e-11]
    [-1.24639788e-10  2.43617120e-03  1.02498504e-05  4.71685852e-12
     -2.19212829e-05  2.94502575e-05  1.14532498e-02 -3.28822545e-03
      4.41476047e-03 -1.14532483e-02 -3.28822477e-03  4.41475996e-03
     -8.08146241e-11 -2.20627126e-02 -9.41894822e-03]
    [ 8.18469493e-11  1.02498504e-05  2.14233489e-03  7.20976151e-13
      2.90464884e-05 -8.63468221e-04 -4.79807048e-03  3.01867574e-03
     -3.68101844e-03  4.79806952e-03  3.01867515e-03 -3.68101831e-03
     -4.67671617e-11 -6.61690630e-03 -4.46017443e-03]
    [-2.16941229e-05  4.71685852e-12  7.20976151e-13  1.49619166e-04
     -2.14462932e-12  8.63862652e-13 -3.85034452e-04 -5.25094481e-04
     -1.04652467e-03 -3.85034497e-04  5.25094463e-04  1.04652462e-03
     -1.35573946e-03 -4.80036012e-12  2.12546863e-11]
    [ 4.93482678e-12 -2.19212829e-05  2.90464884e-05 -2.14462932e-12
      1.50336368e-04 -3.08725828e-05 -5.36109419e-04 -1.06913048e-03
      6.27011757e-04  5.36109379e-04 -1.06913050e-03  6.27011688e-04
      1.39759271e-11  3.87441861e-06 -1.11023649e-03]
    [-7.23564494e-12  2.94502575e-05 -8.63468221e-04  8.63862652e-13
     -3.08725828e-05  1.04165349e-03 -7.28706914e-04  4.42296863e-04
     -2.04725214e-03  7.28706967e-04  4.42296858e-04 -2.04725210e-03
      1.96993795e-11 -7.46717127e-04 -2.14861643e-03]
    [-1.60687481e-02  1.14532498e-02 -4.79807048e-03 -3.85034452e-04
     -5.36109419e-04 -7.28706914e-04  2.30394724e-01 -1.34285111e-01
      7.67026070e-02 -2.16065901e-02  5.58885425e-03 -1.11180651e-02
     -1.12833086e-02  6.64825115e-04  3.15373400e-03]
    [ 1.15383294e-02 -3.28822545e-03  3.01867574e-03 -5.25094481e-04
     -1.06913048e-03  4.42296863e-04 -1.34285111e-01  8.15853116e-02
     -4.64977897e-02 -5.58885465e-03 -7.11299026e-03 -1.05553101e-02
      1.06516082e-02 -1.84750280e-02  1.40635542e-02]
    [-7.21385604e-03  4.41476047e-03 -3.68101844e-03 -1.04652467e-03
      6.27011757e-04 -2.04725214e-03  7.67026070e-02 -4.64977897e-02
      6.64040737e-02  1.11180651e-02 -1.05553100e-02  4.44753123e-03
      1.47911832e-02 -5.53928385e-03  5.49673781e-03]
    [-1.60687464e-02 -1.14532483e-02  4.79806952e-03 -3.85034497e-04
      5.36109379e-04  7.28706967e-04 -2.16065901e-02 -5.58885465e-03
      1.11180651e-02  2.30394705e-01  1.34285095e-01 -7.67025966e-02
     -1.12833095e-02 -6.64825208e-04 -3.15373371e-03]
    [-1.15383280e-02 -3.28822477e-03  3.01867515e-03  5.25094463e-04
     -1.06913050e-03  4.42296858e-04  5.58885425e-03 -7.11299026e-03
     -1.05553100e-02  1.34285095e-01  8.15853035e-02 -4.64977825e-02
     -1.06516072e-02 -1.84750276e-02  1.40635540e-02]
    [ 7.21385523e-03  4.41475996e-03 -3.68101831e-03  1.04652462e-03
      6.27011688e-04 -2.04725210e-03 -1.11180651e-02 -1.05553101e-02
      4.44753123e-03 -7.67025966e-02 -4.64977825e-02  6.64040713e-02
     -1.47911830e-02 -5.53928367e-03  5.49673789e-03]
    [ 3.47625118e-03 -8.08146241e-11 -4.67671617e-11 -1.35573946e-03
      1.39759271e-11  1.96993795e-11 -1.12833086e-02  1.06516082e-02
      1.47911832e-02 -1.12833095e-02 -1.06516072e-02 -1.47911830e-02
      2.46727295e-03 -2.61918768e-10  1.79709160e-11]
    [ 4.94064021e-11 -2.20627126e-02 -6.61690630e-03 -4.80036012e-12
      3.87441861e-06 -7.46717127e-04  6.64825115e-04 -1.84750280e-02
     -5.53928385e-03 -6.64825208e-04 -1.84750276e-02 -5.53928367e-03
     -2.61918768e-10  2.99742686e-01  1.01775207e-01]
    [-4.51157758e-11 -9.41894822e-03 -4.46017443e-03  2.12546863e-11
     -1.11023649e-03 -2.14861643e-03  3.15373400e-03  1.40635542e-02
      5.49673781e-03 -3.15373371e-03  1.40635540e-02  5.49673789e-03
      1.79709160e-11  1.01775207e-01  7.62490460e-02]]
