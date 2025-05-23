Worked Example
--------------

Generating coskf files from xyz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The example will first load all the molecules in the folder ``compounds_xyz`` and then optimize the gas geometry using ADF, and perform the ADF COSMO calculation for each compound. When the calculations are finished, we will find all the .coskf files in the ``test_coskfs_xyz`` directory.

.. code:: ipython3

   from scm.plams import from_smiles, read_molecules, init, JobRunner, config
   from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob
   import os

   # this line is not required in AMS2025+
   init();

::

   PLAMS working folder: /path/plams/examples/COSMORSCompound/plams_workdir

Enable the parallel calculation through ``JobRunner``. Here, we’ll assign one core to each job, and we can have up to eight jobs running all at once.

.. code:: ipython3

   config.default_jobrunner = JobRunner(parallel=True, maxjobs=8)  # Set the default jobrunner to be parallel
   config.default_jobmanager.settings.hashing = None  # Disable rerun prevention
   config.job.runscript.nproc = 1  # Number of cores for each job
   config.log.stdout = 1  # Suppress plams output

.. code:: ipython3

   molecules = read_molecules("./compounds_xyz")

   results = []
   for name, mol in molecules.items():
       job = ADFCOSMORSCompoundJob(molecule=mol, coskf_name=name, coskf_dir="test_coskfs_xyz")
       results.append(job.run())

::

   [19.03|15:14:40] JOB plamsjob STARTED
   [19.03|15:14:40] JOB plamsjob STARTED
   [19.03|15:14:40] JOB plamsjob/gas STARTED
   [19.03|15:14:40] JOB plamsjob/solv STARTED
   [19.03|15:14:40] JOB plamsjob/sigma STARTED

.. code:: ipython3

   for result in results:
       result.wait()

::

   [19.03|15:14:40] Waiting for job plamsjob to finish
   [19.03|15:14:40] JOB plamsjob.002/gas STARTED
   [19.03|15:14:40] JOB plamsjob.002/solv STARTED
   [19.03|15:14:40] Waiting for job gas to finish
   [19.03|15:14:40] JOB plamsjob.002/sigma STARTED
   [19.03|15:14:40] Waiting for job solv to finish
   [19.03|15:14:40] Waiting for job gas to finish
   [19.03|15:14:40] Waiting for job solv to finish
   [19.03|15:14:45] JOB plamsjob.002/gas SUCCESSFUL
   [19.03|15:14:49] JOB plamsjob.002/solv SUCCESSFUL
   [19.03|15:14:49] JOB plamsjob.002/sigma SUCCESSFUL
   [19.03|15:14:50] JOB plamsjob.002 SUCCESSFUL
   [19.03|15:14:57] JOB plamsjob/gas SUCCESSFUL
   [19.03|15:15:10] JOB plamsjob/solv SUCCESSFUL
   [19.03|15:15:10] JOB plamsjob/sigma SUCCESSFUL
   ... (PLAMS log lines truncated) ...

.. code:: ipython3

   print(f"coskf files generated: {', '.join([f for f in os.listdir('./test_coskfs_xyz')])}")

::

   coskf files generated: CO.coskf, H2O.coskf

Generating .coskf files from smiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, we will specify the smiles and name of a set of compounds and generate the initial geometry of each compound using ``from_smiles`` function. With the setting, ``nconfs=100`` and ``forcefield='uff'``, we will generate 100 conformers and find the one with the lowest energy using ‘uff’ forcefield. When the calculations are finished, we will find all the .coskf file in the ``test_coskfs_smiles`` directory.

.. code:: ipython3

   rd_smiles = ["O", "CO"]
   rd_names = ["H2O", "CO"]
   molecules = {}
   for name, smiles in zip(rd_names, rd_smiles):
       molecules[name] = from_smiles(smiles, nconfs=100, forcefield="uff")[0]  # lowest energy one in 100 conformers

Lastly, we give this information to the ``ADFCOSMORSCompoundJob`` class, including the name of the coskf files as well as the directory in which we’ll find them after the calculations complete. Using the setting, ``preoptimization='GFN1-xTB'`` and ``singlepoint=False``, it will utilize the DFTB for a quick pre-optimization. Subsequently, it will execute a gas phase optimization using ADF, followed by the solvation calculation.

.. code:: ipython3

   results = []
   for name, mol in molecules.items():
       job = ADFCOSMORSCompoundJob(
           molecule=mol,  # The initial structure
           coskf_name=name,  # a name to be used for coskf file
           coskf_dir="test_coskfs_smiles",  # a directory to put the .coskf files generated
           preoptimization="GFN1-xTB",  # perform preoptimize or not
           singlepoint=False,  # run a singlepoint in gasphase and solvation calculation without geometry optimization. Cannot be combined with `preoptimization`
           name=name,
       )  # an optional name for the calculation directory
       results.append(job.run())

::

   [19.03|15:15:16] JOB H2O STARTED
   [19.03|15:15:16] JOB CO STARTED
   [19.03|15:15:16] JOB H2O/preoptimization STARTED
   [19.03|15:15:16] JOB CO/preoptimization STARTED

.. code:: ipython3

   for result in results:
       result.wait()

::

   [19.03|15:15:16] Waiting for job H2O to finish
   [19.03|15:15:16] JOB CO/gas STARTED
   [19.03|15:15:16] JOB H2O/gas STARTED
   [19.03|15:15:16] JOB H2O/solv STARTED
   [19.03|15:15:16] JOB CO/solv STARTED
   [19.03|15:15:16] JOB H2O/sigma STARTED
   [19.03|15:15:16] JOB CO/sigma STARTED
   [19.03|15:15:16] Waiting for job gas to finish
   [19.03|15:15:16] Waiting for job preoptimization to finish
   [19.03|15:15:16] Waiting for job preoptimization to finish
   [19.03|15:15:16] Waiting for job gas to finish
   [19.03|15:15:16] Waiting for job solv to finish
   [19.03|15:15:16] Waiting for job solv to finish
   [19.03|15:15:16] JOB H2O/preoptimization SUCCESSFUL
   [19.03|15:15:16] JOB CO/preoptimization SUCCESSFUL
   [19.03|15:15:23] JOB H2O/gas SUCCESSFUL
   [19.03|15:15:27] JOB H2O/solv SUCCESSFUL
   ... (PLAMS log lines truncated) ...
   [19.03|15:15:31] Waiting for job CO to finish

.. code:: ipython3

   print(f"coskf files generated: {', '.join([f for f in os.listdir('./test_coskfs_smiles')])}")

::

   coskf files generated: CO.coskf, H2O.coskf
