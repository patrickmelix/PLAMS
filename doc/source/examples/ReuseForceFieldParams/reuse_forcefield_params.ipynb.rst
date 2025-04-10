Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

   import sys

   from scm.plams import AMSJob, Settings, init, from_smiles

.. code:: ipython3

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/examples/plams_workdir

Run/Load Job with ForceField Information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First run a reference calculation where charges are guessed (using DFTB by default):

.. code:: ipython3

   ref_job = AMSJob.from_input(
       """
   Task GeometryOptimization

   GeometryOptimization
      Convergence Step=1.0e-3
   End

   System
      Atoms
         C 0.0 0.0 0.0
         O 1.13 0.0 0.0
         C 0.0 0.0 2.1
         O 1.13 0.0 1.9
      End
   End

   Engine ForceField
      Verbosity Verbose
      GuessCharges True
   EndEngine
   """
   )

.. code:: ipython3

   ref_job.run();

::

   [18.03|13:53:01] JOB plamsjob STARTED
   [18.03|13:53:01] JOB plamsjob RUNNING
   [18.03|13:53:02] JOB plamsjob FINISHED
   [18.03|13:53:02] JOB plamsjob SUCCESSFUL

.. code:: ipython3

   # Alternatively, load a previously run calculation
   # ref_job = AMSJob.load_external("./plams_workdir/plamsjob/ams.rkf")

Reuse ForceField Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extract the charges and types from the job results and add them as properties on the molecule:

.. code:: ipython3

   charges, types, patch = ref_job.results.get_forcefield_params()

.. code:: ipython3

   mol = ref_job.molecule[""].copy()

   for i, at in enumerate(mol.atoms):
       at.properties.ForceField.Charge = charges[i]
       at.properties.ForceField.Type = types[i]

.. code:: ipython3

   sett = Settings()
   sett.input.AMS.Task = "SinglePoint"
   sett.input.ForceField.Type = "UFF"

.. code:: ipython3

   # Create a patch file if required
   if patch:
       with open("patch.dat", "w") as outfile:
           outfile.write(str(patch))
           outfile.close()
       # For example with:
       # sett.input.ForceField.GAFF.ForceFieldPatchFile = "patch.dat"

.. code:: ipython3

   job = AMSJob(molecule=mol, settings=sett)

.. code:: ipython3

   print(job.get_input())

::

   Task SinglePoint

   System
     Atoms
                 C       0.0000000000       0.0000000000       0.0000000000 ForceField.Charge=0.2881959744167275 ForceField.Type=C_1
                 O       1.1300000000       0.0000000000       0.0000000000 ForceField.Charge=-0.2676126103828702 ForceField.Type=O_2
                 C       0.0000000000       0.0000000000       2.1000000000 ForceField.Charge=0.2536150412119178 ForceField.Type=C_1
                 O       1.1300000000       0.0000000000       1.9000000000 ForceField.Charge=-0.27419840524497996 ForceField.Type=O_2
     End
   End

   Engine ForceField
     Type UFF
   EndEngine

.. code:: ipython3

   job.run();

::

   [18.03|13:53:02] JOB plamsjob STARTED
   [18.03|13:53:02] Renaming job plamsjob to plamsjob.002
   [18.03|13:53:02] JOB plamsjob.002 RUNNING
   [18.03|13:53:03] JOB plamsjob.002 FINISHED
   [18.03|13:53:03] JOB plamsjob.002 SUCCESSFUL
