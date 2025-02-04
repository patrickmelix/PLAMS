Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

<<<<<<< HEAD
   import multiprocessing
   from scm.plams import JobRunner, config, from_smiles, Settings, AMSJob, init

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/examples/BasisSetBenchmark/plams_workdir
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
    import multiprocessing
    from scm.plams import JobRunner, config, from_smiles, Settings, AMSJob
=======
    import sys
    import multiprocessing
    from scm.plams import JobRunner, config, from_smiles, Settings, AMSJob, init
    import numpy as np
    
    init();  # this line is not required in AMS2025+
>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

Set Up Job Runner
~~~~~~~~~~~~~~~~~

Set up job runner, running as many jobs as possible in parallel.

.. code:: ipython3

   config.default_jobrunner = JobRunner(parallel=True, maxjobs=multiprocessing.cpu_count())
   config.job.runscript.nproc = 1

Set Up Molecules
~~~~~~~~~~~~~~~~

Create the molecules we want to use in our benchmark from SMILES.

.. code:: ipython3

   # The molecules we want to use in our benchmark:
   mol_smiles = {"Methane": "C", "Ethane": "C-C", "Ethylene": "C=C", "Acetylene": "C#C"}
   molecules = {}
   for name, smiles in mol_smiles.items():
       # Compute 10 conformers, optimize with UFF and pick the lowest in energy.
       molecules[name] = from_smiles(smiles, nconfs=10, forcefield="uff")[0]
       print(name, molecules[name])

::

   Methane   Atoms: 
       1         C       0.000000       0.000000       0.000000
       2         H       0.538912       0.762358      -0.599295
       3         H       0.731244      -0.596616       0.583182
       4         H      -0.567129      -0.670302      -0.678108
       5         H      -0.703028       0.504560       0.694220
     Bonds: 
      (1)--1.0--(2)
      (1)--1.0--(3)
      (1)--1.0--(4)
      (1)--1.0--(5)

   Ethane   Atoms: 
       1         C      -0.757196      -0.040522       0.044605
       2         C       0.757196       0.040522      -0.044605
       3         H      -1.205222       0.185290      -0.945970
       4         H      -1.130281       0.694397       0.788688
       5         H      -1.061719      -1.061491       0.357407
       6         H       1.205222      -0.185290       0.945971
       7         H       1.130281      -0.694396      -0.788689
       8         H       1.061719       1.061491      -0.357406
     Bonds: 
      (1)--1.0--(2)
      (1)--1.0--(3)
      (1)--1.0--(4)
      (1)--1.0--(5)
      (2)--1.0--(6)
      (2)--1.0--(7)
      (2)--1.0--(8)

   Ethylene   Atoms: 
       1         C       0.664485       0.027988      -0.023685
       2         C      -0.664485      -0.027988       0.023685
       3         H       1.253433      -0.878614       0.070299
       4         H       1.167038       0.980564      -0.156575
       5         H      -1.253433       0.878614      -0.070299
       6         H      -1.167038      -0.980564       0.156575
     Bonds: 
      (1)--2.0--(2)
      (1)--1.0--(3)
      (1)--1.0--(4)
      (2)--1.0--(5)
      (2)--1.0--(6)

   Acetylene   Atoms: 
       1         C      -0.587409       0.175060      -0.002211
       2         C       0.587409      -0.094463       0.002211
       3         H      -1.618985       0.411721      -0.006095
       4         H       1.618985      -0.331124       0.006094
     Bonds: 
      (1)--3.0--(2)
      (1)--1.0--(3)
      (2)--1.0--(4)

Initialize Calculation Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Set up the settings which are common across jobs. The basis type is added later for each job.

.. code:: ipython3

   common_settings = Settings()
   common_settings.input.ams.Task = "SinglePoint"
   common_settings.input.ams.System.Symmetrize = "Yes"
   common_settings.input.adf.Basis.Core = "None"

.. code:: ipython3

   basis = ["QZ4P", "TZ2P", "TZP", "DZP", "DZ", "SZ"]
   reference_basis = "QZ4P"

Run Calculations
~~~~~~~~~~~~~~~~

.. code:: ipython3

<<<<<<< HEAD
   results = {}
   for bas in basis:
       for name, molecule in molecules.items():
           settings = common_settings.copy()
           settings.input.adf.Basis.Type = bas
           job = AMSJob(name=name + "_" + bas, molecule=molecule, settings=settings)
           results[(name, bas)] = job.run()
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
    results = {}
    for bas in basis:
        for name, molecule in molecules.items():
            settings = common_settings.copy()
            settings.input.adf.Basis.Type = bas
            job = AMSJob(name=name + "_" + bas, molecule=molecule, settings=settings)
            results[(name, bas)] = job.run()
=======
    results = {}
    jobs = []
    for bas in basis:
        for name, molecule in molecules.items():
            settings = common_settings.copy()
            settings.input.adf.Basis.Type = bas
            job = AMSJob(name=name + "_" + bas, molecule=molecule, settings=settings)
            jobs.append(job)
            results[(name, bas)] = job.run()
>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

::

<<<<<<< HEAD
   [10.02|15:01:11] JOB Methane_QZ4P STARTED
   [10.02|15:01:11] JOB Ethane_QZ4P STARTED
   [10.02|15:01:11] JOB Ethylene_QZ4P STARTED
   [10.02|15:01:11] JOB Methane_QZ4P RUNNING
   [10.02|15:01:11] JOB Acetylene_QZ4P STARTED
   [10.02|15:01:11] JOB Methane_TZ2P STARTED
   [10.02|15:01:11] JOB Ethane_TZ2P STARTED
   [10.02|15:01:11] JOB Ethylene_TZ2P STARTED
   [10.02|15:01:11] JOB Acetylene_TZ2P STARTED
   [10.02|15:01:11] JOB Methane_TZP STARTED
   ... (PLAMS log lines truncated) ...
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
.. parsed-literal::

    [13.08|15:09:01] JOB Methane_QZ4P STARTED
    [13.08|15:09:01] JOB Ethane_QZ4P STARTED
    [13.08|15:09:01] JOB Ethylene_QZ4P STARTED
    [13.08|15:09:01] JOB Acetylene_QZ4P STARTED
    [13.08|15:09:01] JOB Methane_TZ2P STARTED
    [13.08|15:09:01] JOB Methane_QZ4P RUNNING
    [13.08|15:09:01] JOB Ethane_TZ2P STARTED
    [13.08|15:09:01] JOB Ethane_QZ4P RUNNING
    [13.08|15:09:01] JOB Ethylene_TZ2P STARTED
    [13.08|15:09:01] JOB Ethylene_QZ4P RUNNING
    [13.08|15:09:01] JOB Acetylene_TZ2P STARTED
    [13.08|15:09:01] JOB Acetylene_QZ4P RUNNING
    [13.08|15:09:01] JOB Methane_TZ2P RUNNING
    [13.08|15:09:01] JOB Methane_TZP STARTED
    [13.08|15:09:01] JOB Ethane_TZ2P RUNNING
    [13.08|15:09:01] JOB Ethylene_TZ2P RUNNING
    [13.08|15:09:01] JOB Ethane_TZP STARTED
    [13.08|15:09:01] JOB Ethylene_TZP STARTED
    [13.08|15:09:01] JOB Acetylene_TZ2P RUNNING
    [13.08|15:09:01] JOB Acetylene_TZP STARTED
    [13.08|15:09:01] JOB Methane_DZP STARTED
    [13.08|15:09:01] JOB Methane_TZP RUNNING
    [13.08|15:09:01] JOB Ethane_DZP STARTED
    [13.08|15:09:01] JOB Ethane_TZP RUNNING
    [13.08|15:09:01] JOB Ethylene_DZP STARTED
    [13.08|15:09:01] JOB Ethylene_TZP RUNNING
    [13.08|15:09:01] JOB Acetylene_DZP STARTED
    [13.08|15:09:01] JOB Acetylene_TZP RUNNING
    [13.08|15:09:01] JOB Methane_DZ STARTED
    [13.08|15:09:01] JOB Methane_DZP RUNNING
    [13.08|15:09:01] JOB Ethane_DZ STARTED
    [13.08|15:09:01] JOB Ethane_DZP RUNNING
    [13.08|15:09:01] JOB Ethylene_DZ STARTED
    [13.08|15:09:01] JOB Ethylene_DZP RUNNING
    [13.08|15:09:01] JOB Acetylene_DZ STARTED
    [13.08|15:09:01] JOB Acetylene_DZP RUNNING
    [13.08|15:09:01] JOB Methane_SZ STARTED
    [13.08|15:09:01] JOB Methane_DZ RUNNING
    [13.08|15:09:01] JOB Ethane_SZ STARTED
    [13.08|15:09:01] JOB Ethane_DZ RUNNING
    [13.08|15:09:01] JOB Ethylene_SZ STARTED
    [13.08|15:09:01] JOB Ethylene_DZ RUNNING
    [13.08|15:09:01] JOB Acetylene_DZ RUNNING
    [13.08|15:09:01] JOB Acetylene_SZ STARTED
    [13.08|15:09:01] JOB Methane_SZ RUNNING
    [13.08|15:09:01] JOB Ethane_SZ RUNNING

=======
.. parsed-literal::

    [04.02|17:21:55] JOB Methane_QZ4P STARTED
    [04.02|17:21:55] JOB Ethane_QZ4P STARTED
    [04.02|17:21:55] JOB Ethylene_QZ4P STARTED
    [04.02|17:21:55] JOB Acetylene_QZ4P STARTED
    [04.02|17:21:55] JOB Methane_QZ4P RUNNING
    [04.02|17:21:55] JOB Methane_TZ2P STARTED
    [04.02|17:21:55] JOB Ethane_QZ4P RUNNING
    [04.02|17:21:55] JOB Ethane_TZ2P STARTED
    [04.02|17:21:55] JOB Ethylene_QZ4P RUNNING
    [04.02|17:21:55] JOB Ethylene_TZ2P STARTED
    [04.02|17:21:55] JOB Methane_TZ2P RUNNING
    [04.02|17:21:55] JOB Acetylene_TZ2P STARTED
    [04.02|17:21:55] JOB Acetylene_QZ4P RUNNING
    [04.02|17:21:55] JOB Acetylene_TZ2P RUNNING
    [04.02|17:21:55] JOB Methane_TZP STARTED
    [04.02|17:21:55] JOB Ethane_TZP STARTED
    [04.02|17:21:55] JOB Ethane_TZ2P RUNNING
    [04.02|17:21:55] JOB Ethylene_TZP STARTED
    [04.02|17:21:55] JOB Acetylene_TZP STARTED
    [04.02|17:21:55] JOB Methane_DZP STARTED
    [04.02|17:21:55] JOB Ethane_DZP STARTED
    [04.02|17:21:55] JOB Ethylene_TZ2P RUNNING
    [04.02|17:21:55] JOB Ethylene_DZP STARTED
    [04.02|17:21:55] JOB Acetylene_DZP STARTED
    [04.02|17:21:55] JOB Methane_DZ STARTED
    [04.02|17:21:55] JOB Ethane_DZ STARTED
    [04.02|17:21:55] JOB Ethylene_DZ STARTED
    [04.02|17:21:55] JOB Acetylene_DZ STARTED
    [04.02|17:21:55] JOB Ethylene_TZP RUNNING
    [04.02|17:21:55] JOB Methane_SZ STARTED
    [04.02|17:21:55] JOB Acetylene_TZP RUNNING
    [04.02|17:21:55] JOB Ethane_SZ STARTED
    [04.02|17:21:55] JOB Ethylene_SZ STARTED
    [04.02|17:21:55] JOB Methane_TZP RUNNING
    [04.02|17:21:55] JOB Acetylene_SZ STARTED
    [04.02|17:21:55] JOB Acetylene_DZP RUNNING
    [04.02|17:21:55] JOB Methane_DZP RUNNING
    [04.02|17:21:55] JOB Ethylene_DZP RUNNING
    [04.02|17:21:55] JOB Ethane_DZP RUNNING
    [04.02|17:21:55] JOB Ethylene_SZ RUNNING
    [04.02|17:21:55] JOB Ethane_TZP RUNNING
    [04.02|17:21:55] JOB Acetylene_DZ RUNNING
    [04.02|17:21:55] JOB Methane_DZ RUNNING
    [04.02|17:21:55] JOB Ethane_DZ RUNNING
    [04.02|17:21:55] JOB Methane_SZ RUNNING
    [04.02|17:21:55] JOB Ethylene_DZ RUNNING
    [04.02|17:21:55] JOB Acetylene_SZ RUNNING
    [04.02|17:21:55] JOB Ethane_SZ RUNNING
    [04.02|17:21:59] JOB Methane_TZ2P FINISHED
    [04.02|17:21:59] JOB Methane_TZ2P SUCCESSFUL
    [04.02|17:21:59] JOB Methane_QZ4P FINISHED
    [04.02|17:21:59] JOB Methane_QZ4P SUCCESSFUL
    [04.02|17:21:59] JOB Methane_TZP FINISHED
    [04.02|17:21:59] JOB Methane_TZP SUCCESSFUL
    [04.02|17:22:00] JOB Acetylene_TZP FINISHED
    [04.02|17:22:00] JOB Acetylene_TZP SUCCESSFUL
    [04.02|17:22:00] JOB Ethylene_TZ2P FINISHED
    [04.02|17:22:00] JOB Ethylene_TZ2P SUCCESSFUL
    [04.02|17:22:00] JOB Acetylene_TZ2P FINISHED
    [04.02|17:22:00] JOB Acetylene_TZ2P SUCCESSFUL
    [04.02|17:22:00] JOB Ethylene_TZP FINISHED
    [04.02|17:22:00] JOB Acetylene_DZP FINISHED
    [04.02|17:22:00] JOB Acetylene_DZP SUCCESSFUL
    [04.02|17:22:00] JOB Acetylene_QZ4P FINISHED
    [04.02|17:22:00] JOB Ethylene_TZP SUCCESSFUL
    [04.02|17:22:00] JOB Ethylene_QZ4P FINISHED
    [04.02|17:22:00] JOB Acetylene_QZ4P SUCCESSFUL
    [04.02|17:22:00] JOB Ethylene_QZ4P SUCCESSFUL
    [04.02|17:22:03] JOB Ethane_TZ2P FINISHED
    [04.02|17:22:03] JOB Ethane_TZ2P SUCCESSFUL
    [04.02|17:22:03] JOB Acetylene_SZ FINISHED
    [04.02|17:22:03] JOB Acetylene_SZ SUCCESSFUL
    [04.02|17:22:04] JOB Acetylene_DZ FINISHED
    [04.02|17:22:04] JOB Acetylene_DZ SUCCESSFUL
    [04.02|17:22:04] JOB Methane_DZP FINISHED
    [04.02|17:22:04] JOB Methane_DZP SUCCESSFUL
    [04.02|17:22:04] JOB Ethane_QZ4P FINISHED
    [04.02|17:22:04] JOB Ethane_QZ4P SUCCESSFUL
    [04.02|17:22:04] JOB Methane_DZ FINISHED
    [04.02|17:22:04] JOB Methane_SZ FINISHED
    [04.02|17:22:04] JOB Methane_DZ SUCCESSFUL
    [04.02|17:22:04] JOB Methane_SZ SUCCESSFUL
    [04.02|17:22:05] JOB Ethylene_DZ FINISHED
    [04.02|17:22:05] JOB Ethylene_DZ SUCCESSFUL
    [04.02|17:22:05] JOB Ethane_DZP FINISHED
    [04.02|17:22:05] JOB Ethane_SZ FINISHED
    [04.02|17:22:05] JOB Ethane_DZP SUCCESSFUL
    [04.02|17:22:05] JOB Ethane_SZ SUCCESSFUL
    [04.02|17:22:05] JOB Ethylene_DZP FINISHED
    [04.02|17:22:05] JOB Ethylene_DZP SUCCESSFUL
    [04.02|17:22:05] JOB Ethane_DZ FINISHED
    [04.02|17:22:05] JOB Ethane_DZ SUCCESSFUL
    [04.02|17:22:05] JOB Ethylene_SZ FINISHED
    [04.02|17:22:05] JOB Ethylene_SZ SUCCESSFUL
    [04.02|17:22:06] JOB Ethane_TZP FINISHED
    [04.02|17:22:06] JOB Ethane_TZP SUCCESSFUL

>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

Results
~~~~~~~

Extract the energy from each calculation. Calculate the average absolute error in bond energy per atom for each basis set.

.. code:: ipython3

<<<<<<< HEAD
   average_errors = {}
   for bas in basis:
       if bas != reference_basis:
           errors = []
           for name, molecule in molecules.items():
               reference_energy = results[(name, reference_basis)].get_energy(unit="kcal/mol")
               energy = results[(name, bas)].get_energy(unit="kcal/mol")
               errors.append(abs(energy - reference_energy) / len(molecule))
               print("Energy for {} using {} basis set: {} [kcal/mol]".format(name, bas, energy))
           average_errors[bas] = sum(errors) / len(errors)
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
    average_errors = {}
    for bas in basis:
        if bas != reference_basis:
            errors = []
            for name, molecule in molecules.items():
                reference_energy = results[(name, reference_basis)].get_energy(unit="kcal/mol")
                energy = results[(name, bas)].get_energy(unit="kcal/mol")
                errors.append(abs(energy - reference_energy) / len(molecule))
                print("Energy for {} using {} basis set: {} [kcal/mol]".format(name, bas, energy))
            average_errors[bas] = sum(errors) / len(errors)
=======
    try:
        # For AMS2025+ can use JobAnalysis class to perform results analysis
        from scm.plams import JobAnalysis
    
        ja = (
            JobAnalysis(jobs=jobs, std_fields=None)
            .add_formula_field()
            .add_smiles_field()
            .add_settings_field(("Input", "ADF", "Basis", "Type"), display_name="Basis")
            .add_field("NAtoms", lambda j: len(j.molecule))
            .add_field(
                "Energy", lambda j: j.results.get_energy(unit="kcal/mol"), display_name="Energy [kcal/mol]", fmt=".2f"
            )
            .sort_jobs(["NAtoms", "Energy"])
        )
    
        ref_ja = ja.copy().filter_jobs(lambda data: data["InputAdfBasisType"] == "QZ4P")
    
        ref_energies = {f: e for f, e in zip(ref_ja.Formula, ref_ja.Energy)}
    
        def get_average_error(job):
            return abs(job.results.get_energy(unit="kcal/mol") - ref_energies[job.molecule.get_formula()]) / len(
                job.molecule
            )
    
        ja.add_field("AvErr", get_average_error, display_name="Average Error [kcal/mol]", fmt=".2f")
    
        # Pretty-print if running in a notebook
        if "ipykernel" in sys.modules:
            ja.display_table()
        else:
            print(ja.to_table())
    
    except ImportError:
    
        average_errors = {}
        for bas in basis:
            if bas != reference_basis:
                errors = []
                for name, molecule in molecules.items():
                    reference_energy = results[(name, reference_basis)].get_energy(unit="kcal/mol")
                    energy = results[(name, bas)].get_energy(unit="kcal/mol")
                    errors.append(abs(energy - reference_energy) / len(molecule))
                    print("Energy for {} using {} basis set: {} [kcal/mol]".format(name, bas, energy))
                average_errors[bas] = sum(errors) / len(errors)
>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

::

<<<<<<< HEAD
   [10.02|15:01:11] JOB Acetylene_TZP RUNNING
   [10.02|15:01:11] Waiting for job Methane_QZ4P to finish
   [10.02|15:01:11] JOB Methane_DZ RUNNING
   [10.02|15:01:11] JOB Ethane_DZP RUNNING
   [10.02|15:01:11] JOB Ethylene_DZP RUNNING
   [10.02|15:01:11] JOB Acetylene_DZP RUNNING
   [10.02|15:01:11] JOB Methane_SZ RUNNING
   [10.02|15:01:11] JOB Ethane_DZ RUNNING
   [10.02|15:01:11] JOB Ethane_SZ RUNNING
   [10.02|15:01:11] JOB Ethylene_DZ RUNNING
   [10.02|15:01:11] JOB Acetylene_DZ RUNNING
   ... (PLAMS log lines truncated) ...
   Energy for Methane using TZ2P basis set: -572.110159165262 [kcal/mol]
   [10.02|15:01:15] Waiting for job Ethane_QZ4P to finish
   Energy for Ethane using TZ2P basis set: -971.8820186844459 [kcal/mol]
   Energy for Ethylene using TZ2P basis set: -769.4329031250381 [kcal/mol]
   Energy for Acetylene using TZ2P basis set: -555.6672902509043 [kcal/mol]
   Energy for Methane using TZP basis set: -571.0448969099549 [kcal/mol]
   Energy for Ethane using TZP basis set: -970.0758887573307 [kcal/mol]
   Energy for Ethylene using TZP basis set: -767.3275176578105 [kcal/mol]
   [10.02|15:01:18] Waiting for job Acetylene_TZP to finish
   Energy for Acetylene using TZP basis set: -552.9562856742521 [kcal/mol]
   Energy for Methane using DZP basis set: -569.1190156251034 [kcal/mol]
   [10.02|15:01:19] Waiting for job Ethane_DZP to finish
   Energy for Ethane using DZP basis set: -966.0916443143674 [kcal/mol]
   Energy for Ethylene using DZP basis set: -764.4132984010868 [kcal/mol]
   Energy for Acetylene using DZP basis set: -550.6461805496328 [kcal/mol]
   Energy for Methane using DZ basis set: -560.9344313073021 [kcal/mol]
   [10.02|15:01:19] Waiting for job Ethane_DZ to finish
   Energy for Ethane using DZ basis set: -951.1666971758781 [kcal/mol]
   Energy for Ethylene using DZ basis set: -750.1745108422972 [kcal/mol]
   Energy for Acetylene using DZ basis set: -537.1008020388887 [kcal/mol]
   Energy for Methane using SZ basis set: -723.550123154895 [kcal/mol]
   Energy for Ethane using SZ basis set: -1216.914233427825 [kcal/mol]
   Energy for Ethylene using SZ basis set: -934.6558200110123 [kcal/mol]
   Energy for Acetylene using SZ basis set: -647.5029836817757 [kcal/mol]
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
.. parsed-literal::

    [13.08|15:09:01] JOB Ethylene_SZ RUNNING
    [13.08|15:09:01] Waiting for job Methane_QZ4P to finish
    [13.08|15:09:01] JOB Acetylene_SZ RUNNING
    [13.08|15:09:04] JOB Methane_QZ4P FINISHED
    [13.08|15:09:04] JOB Methane_QZ4P SUCCESSFUL
    [13.08|15:09:04] Waiting for job Methane_TZ2P to finish
    [13.08|15:09:04] JOB Methane_TZ2P FINISHED
    [13.08|15:09:04] JOB Methane_TZ2P SUCCESSFUL
    Energy for Methane using TZ2P basis set: -572.1101591652508 [kcal/mol]
    [13.08|15:09:04] Waiting for job Ethane_QZ4P to finish
    [13.08|15:09:04] JOB Methane_TZP FINISHED
    [13.08|15:09:04] JOB Methane_TZP SUCCESSFUL
    [13.08|15:09:05] JOB Ethylene_QZ4P FINISHED
    [13.08|15:09:05] JOB Ethylene_QZ4P SUCCESSFUL
    [13.08|15:09:05] JOB Acetylene_QZ4P FINISHED
    [13.08|15:09:05] JOB Acetylene_QZ4P SUCCESSFUL
    [13.08|15:09:05] JOB Ethylene_TZ2P FINISHED
    [13.08|15:09:05] JOB Ethylene_TZ2P SUCCESSFUL
    [13.08|15:09:05] JOB Ethylene_TZP FINISHED
    [13.08|15:09:05] JOB Ethylene_TZP SUCCESSFUL
    [13.08|15:09:05] JOB Acetylene_TZP FINISHED
    [13.08|15:09:05] JOB Acetylene_TZP SUCCESSFUL
    [13.08|15:09:05] JOB Acetylene_TZ2P FINISHED
    [13.08|15:09:05] JOB Acetylene_TZ2P SUCCESSFUL
    [13.08|15:09:06] JOB Ethane_TZP FINISHED
    [13.08|15:09:06] JOB Ethane_TZP SUCCESSFUL
    [13.08|15:09:06] JOB Methane_DZP FINISHED
    [13.08|15:09:06] JOB Methane_DZP SUCCESSFUL
    [13.08|15:09:06] JOB Ethane_TZ2P FINISHED
    [13.08|15:09:06] JOB Ethane_TZ2P SUCCESSFUL
    [13.08|15:09:07] JOB Ethane_QZ4P FINISHED
    [13.08|15:09:07] JOB Ethane_QZ4P SUCCESSFUL
    Energy for Ethane using TZ2P basis set: -971.8820186845153 [kcal/mol]
    Energy for Ethylene using TZ2P basis set: -769.4329031250394 [kcal/mol]
    Energy for Acetylene using TZ2P basis set: -555.667290250868 [kcal/mol]
    Energy for Methane using TZP basis set: -571.0448969099632 [kcal/mol]
    Energy for Ethane using TZP basis set: -970.0758887574258 [kcal/mol]
    Energy for Ethylene using TZP basis set: -767.3275176577931 [kcal/mol]
    Energy for Acetylene using TZP basis set: -552.956285674204 [kcal/mol]
    Energy for Methane using DZP basis set: -569.1190156251367 [kcal/mol]
    [13.08|15:09:07] Waiting for job Ethane_DZP to finish
    [13.08|15:09:07] JOB Methane_DZ FINISHED
    [13.08|15:09:07] JOB Methane_DZ SUCCESSFUL
    [13.08|15:09:07] JOB Ethylene_DZP FINISHED
    [13.08|15:09:07] JOB Ethylene_DZP SUCCESSFUL
    [13.08|15:09:07] JOB Ethane_DZP FINISHED
    [13.08|15:09:08] JOB Ethane_DZP SUCCESSFUL
    Energy for Ethane using DZP basis set: -966.0916443143979 [kcal/mol]
    Energy for Ethylene using DZP basis set: -764.4132984011687 [kcal/mol]
    [13.08|15:09:08] Waiting for job Acetylene_DZP to finish
    [13.08|15:09:08] JOB Methane_SZ FINISHED
    [13.08|15:09:08] JOB Methane_SZ SUCCESSFUL
    [13.08|15:09:08] JOB Acetylene_DZP FINISHED
    [13.08|15:09:08] JOB Acetylene_DZP SUCCESSFUL
    Energy for Acetylene using DZP basis set: -550.6461805495554 [kcal/mol]
    Energy for Methane using DZ basis set: -560.9344313072968 [kcal/mol]
    [13.08|15:09:08] Waiting for job Ethane_DZ to finish
    [13.08|15:09:08] JOB Ethylene_DZ FINISHED
    [13.08|15:09:08] JOB Ethylene_DZ SUCCESSFUL
    [13.08|15:09:08] JOB Acetylene_DZ FINISHED
    [13.08|15:09:08] JOB Acetylene_DZ SUCCESSFUL
    [13.08|15:09:08] JOB Ethane_DZ FINISHED
    [13.08|15:09:08] JOB Ethane_DZ SUCCESSFUL
    Energy for Ethane using DZ basis set: -951.1666971758054 [kcal/mol]
    Energy for Ethylene using DZ basis set: -750.1745108423067 [kcal/mol]
    Energy for Acetylene using DZ basis set: -537.100802038877 [kcal/mol]
    Energy for Methane using SZ basis set: -723.5501231548906 [kcal/mol]
    [13.08|15:09:08] Waiting for job Ethane_SZ to finish
    [13.08|15:09:09] JOB Ethylene_SZ FINISHED
    [13.08|15:09:09] JOB Ethylene_SZ SUCCESSFUL
    [13.08|15:09:09] JOB Ethane_SZ FINISHED
    [13.08|15:09:09] JOB Ethane_SZ SUCCESSFUL
    Energy for Ethane using SZ basis set: -1216.91423342784 [kcal/mol]
    Energy for Ethylene using SZ basis set: -934.6558200110214 [kcal/mol]
    [13.08|15:09:09] Waiting for job Acetylene_SZ to finish
    [13.08|15:09:09] JOB Acetylene_SZ FINISHED
    [13.08|15:09:09] JOB Acetylene_SZ SUCCESSFUL
    Energy for Acetylene using SZ basis set: -647.50298368177 [kcal/mol]

=======
.. parsed-literal::

    [04.02|17:21:55] Waiting for job Methane_QZ4P to finish
    [04.02|17:21:59] Waiting for job Ethane_QZ4P to finish
    [04.02|17:22:04] Waiting for job Ethane_TZP to finish



.. raw:: html

    <div style="max-width: 100%; overflow-x: auto;">
    <table border="1" style="border-collapse: collapse; width: auto; ">
    <thead><tr><th>Formula<th>Smiles<th>Basis<th>NAtoms<th>Energy [kcal/mol]<th>Average Error [kcal/mol]</th></tr></thead>
    <tbody>
    <tr><td>C2H2   </td><td>C#C   </td><td>DZ   </td><td>4     </td><td>-537.10          </td><td>4.91                    </td></tr>
    <tr><td>C2H2   </td><td>C#C   </td><td>DZP  </td><td>4     </td><td>-550.65          </td><td>1.53                    </td></tr>
    <tr><td>C2H2   </td><td>C#C   </td><td>TZP  </td><td>4     </td><td>-552.96          </td><td>0.95                    </td></tr>
    <tr><td>C2H2   </td><td>C#C   </td><td>TZ2P </td><td>4     </td><td>-555.67          </td><td>0.27                    </td></tr>
    <tr><td>C2H2   </td><td>C#C   </td><td>QZ4P </td><td>4     </td><td>-556.76          </td><td>0.00                    </td></tr>
    <tr><td>C2H2   </td><td>C#C   </td><td>SZ   </td><td>4     </td><td>-647.50          </td><td>22.69                   </td></tr>
    <tr><td>CH4    </td><td>C     </td><td>DZ   </td><td>5     </td><td>-560.93          </td><td>2.34                    </td></tr>
    <tr><td>CH4    </td><td>C     </td><td>DZP  </td><td>5     </td><td>-569.12          </td><td>0.70                    </td></tr>
    <tr><td>CH4    </td><td>C     </td><td>TZP  </td><td>5     </td><td>-571.04          </td><td>0.32                    </td></tr>
    <tr><td>CH4    </td><td>C     </td><td>TZ2P </td><td>5     </td><td>-572.11          </td><td>0.10                    </td></tr>
    <tr><td>CH4    </td><td>C     </td><td>QZ4P </td><td>5     </td><td>-572.63          </td><td>0.00                    </td></tr>
    <tr><td>CH4    </td><td>C     </td><td>SZ   </td><td>5     </td><td>-723.55          </td><td>30.18                   </td></tr>
    <tr><td>C2H4   </td><td>C=C   </td><td>DZ   </td><td>6     </td><td>-750.17          </td><td>3.37                    </td></tr>
    <tr><td>C2H4   </td><td>C=C   </td><td>DZP  </td><td>6     </td><td>-764.41          </td><td>1.00                    </td></tr>
    <tr><td>C2H4   </td><td>C=C   </td><td>TZP  </td><td>6     </td><td>-767.33          </td><td>0.51                    </td></tr>
    <tr><td>C2H4   </td><td>C=C   </td><td>TZ2P </td><td>6     </td><td>-769.43          </td><td>0.16                    </td></tr>
    <tr><td>C2H4   </td><td>C=C   </td><td>QZ4P </td><td>6     </td><td>-770.41          </td><td>0.00                    </td></tr>
    <tr><td>C2H4   </td><td>C=C   </td><td>SZ   </td><td>6     </td><td>-934.66          </td><td>27.37                   </td></tr>
    <tr><td>C2H6   </td><td>CC    </td><td>SZ   </td><td>8     </td><td>-1216.91         </td><td>30.49                   </td></tr>
    <tr><td>C2H6   </td><td>CC    </td><td>DZ   </td><td>8     </td><td>-951.17          </td><td>2.73                    </td></tr>
    <tr><td>C2H6   </td><td>CC    </td><td>DZP  </td><td>8     </td><td>-966.09          </td><td>0.87                    </td></tr>
    <tr><td>C2H6   </td><td>CC    </td><td>TZP  </td><td>8     </td><td>-970.08          </td><td>0.37                    </td></tr>
    <tr><td>C2H6   </td><td>CC    </td><td>TZ2P </td><td>8     </td><td>-971.88          </td><td>0.14                    </td></tr>
    <tr><td>C2H6   </td><td>CC    </td><td>QZ4P </td><td>8     </td><td>-973.02          </td><td>0.00                    </td></tr>
    </tbody>
    </table>
    </div>

>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

.. code:: ipython3

<<<<<<< HEAD
   print("== Results ==")
   print("Average absolute error in bond energy per atom")
   for bas in basis:
       if bas != reference_basis:
           print("Error for basis set {:<4}: {:>10.3f} [kcal/mol]".format(bas, average_errors[bas]))
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
    print("== Results ==")
    print("Average absolute error in bond energy per atom")
    for bas in basis:
        if bas != reference_basis:
            print("Error for basis set {:<4}: {:>10.3f} [kcal/mol]".format(bas, average_errors[bas]))
=======
    print("== Results ==")
    print("Average absolute error in bond energy per atom")
    for bas in basis:
        if bas != reference_basis:
            if ja:
                av = np.average(ja.copy().filter_jobs(lambda data: data["InputAdfBasisType"] == bas).AvErr)
            else:
                av = average_errors[bas]
            print("Error for basis set {:<4}: {:>10.3f} [kcal/mol]".format(bas, av))
>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

::

<<<<<<< HEAD
   == Results ==
   Average absolute error in bond energy per atom
   Error for basis set TZ2P:      0.170 [kcal/mol]
   Error for basis set TZP :      0.537 [kcal/mol]
   Error for basis set DZP :      1.024 [kcal/mol]
   Error for basis set DZ  :      3.339 [kcal/mol]
   Error for basis set SZ  :     27.683 [kcal/mol]
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
.. parsed-literal::

    == Results ==
    Average absolute error in bond energy per atom
    Error for basis set TZ2P:      0.170 [kcal/mol]
    Error for basis set TZP :      0.537 [kcal/mol]
    Error for basis set DZP :      1.024 [kcal/mol]
    Error for basis set DZ  :      3.339 [kcal/mol]
    Error for basis set SZ  :     27.683 [kcal/mol]

=======
.. parsed-literal::

    == Results ==
    Average absolute error in bond energy per atom
    Error for basis set TZ2P:      0.170 [kcal/mol]
    Error for basis set TZP :      0.537 [kcal/mol]
    Error for basis set DZP :      1.024 [kcal/mol]
    Error for basis set DZ  :      3.339 [kcal/mol]
    Error for basis set SZ  :     27.683 [kcal/mol]


>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)
