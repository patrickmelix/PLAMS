Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

<<<<<<< HEAD
   import numpy as np
   from scm.plams import Settings, Molecule, Atom, AMSJob, init

   # this line is not required in AMS2025+
   init()

::

   PLAMS working folder: /path/plams/examples/He2DissociationCurve/plams_workdir
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
    import numpy as np
    from scm.plams import Settings, Molecule, Atom, AMSJob
=======
    import sys
    import numpy as np
    from scm.plams import Settings, Molecule, Atom, AMSJob, init
    
    init();  # this line is not required in AMS2025+


.. parsed-literal::

    PLAMS working folder: /path/plams/examples/He2DissociationCurve/plams_workdir.003

>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

Setup Dimer
~~~~~~~~~~~

Create Helium atoms and an array of interatomic distances at which to run calculation.

.. code:: ipython3

   # type of atoms
   atom1 = "He"
   atom2 = "He"

.. code:: ipython3

   # interatomic distance values
   dmin = 2.2
   dmax = 4.2
   step = 0.2

.. code:: ipython3

   # create a list with interatomic distances
   distances = np.arange(dmin, dmax, step)
   print(distances)

::

   [2.2 2.4 2.6 2.8 3.  3.2 3.4 3.6 3.8 4. ]

Calculation Settings
~~~~~~~~~~~~~~~~~~~~

The calculation settins are stored in a ``Settings`` object.

.. code:: ipython3

   # calculation parameters (single point, TZP/PBE+GrimmeD3)
   sett = Settings()
   sett.input.ams.task = "SinglePoint"
   sett.input.adf.basis.type = "TZP"
   sett.input.adf.xc.gga = "PBE"
   sett.input.adf.xc.dispersion = "Grimme3"

Create and Run Jobs
~~~~~~~~~~~~~~~~~~~

For each interatomic distance, create a Helium dimer molecule with the required geometry then the single point energy calculation job. Run the job and extract the energy.

.. code:: ipython3

<<<<<<< HEAD
   energies = []
   for d in distances:
       mol = Molecule()
       mol.add_atom(Atom(symbol=atom1, coords=(0.0, 0.0, 0.0)))
       mol.add_atom(Atom(symbol=atom2, coords=(d, 0.0, 0.0)))
       job = AMSJob(molecule=mol, settings=sett, name=f"dist_{d:.2f}")
       job.run()
       energies.append(job.results.get_energy(unit="kcal/mol"))
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
    energies = []
    for d in distances:
        mol = Molecule()
        mol.add_atom(Atom(symbol=atom1, coords=(0.0, 0.0, 0.0)))
        mol.add_atom(Atom(symbol=atom2, coords=(d, 0.0, 0.0)))
        job = AMSJob(molecule=mol, settings=sett, name=f"dist_{d:.2f}")
        job.run()
        energies.append(job.results.get_energy(unit="kcal/mol"))
=======
    jobs = []
    for d in distances:
        mol = Molecule()
        mol.add_atom(Atom(symbol=atom1, coords=(0.0, 0.0, 0.0)))
        mol.add_atom(Atom(symbol=atom2, coords=(d, 0.0, 0.0)))
        job = AMSJob(molecule=mol, settings=sett, name=f"dist_{d:.2f}")
        jobs.append(job)
        job.run()
>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

::

<<<<<<< HEAD
   [10.02|15:18:00] JOB dist_2.20 STARTED
   [10.02|15:18:00] JOB dist_2.20 RUNNING
   [10.02|15:18:06] JOB dist_2.20 FINISHED
   [10.02|15:18:06] JOB dist_2.20 SUCCESSFUL
   [10.02|15:18:06] JOB dist_2.40 STARTED
   [10.02|15:18:06] JOB dist_2.40 RUNNING
   [10.02|15:18:10] JOB dist_2.40 FINISHED
   [10.02|15:18:10] JOB dist_2.40 SUCCESSFUL
   [10.02|15:18:10] JOB dist_2.60 STARTED
   [10.02|15:18:10] JOB dist_2.60 RUNNING
   ... (PLAMS log lines truncated) ...
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
.. parsed-literal::

    [26.08|15:26:10] JOB dist_2.20 STARTED
    [26.08|15:26:10] JOB dist_2.20 RUNNING
    [26.08|15:26:11] JOB dist_2.20 FINISHED
    [26.08|15:26:12] JOB dist_2.20 SUCCESSFUL
    [26.08|15:26:12] JOB dist_2.40 STARTED
    [26.08|15:26:12] JOB dist_2.40 RUNNING
    [26.08|15:26:13] JOB dist_2.40 FINISHED
    [26.08|15:26:13] JOB dist_2.40 SUCCESSFUL
    [26.08|15:26:13] JOB dist_2.60 STARTED
    [26.08|15:26:13] JOB dist_2.60 RUNNING
    [26.08|15:26:15] JOB dist_2.60 FINISHED
    [26.08|15:26:15] JOB dist_2.60 SUCCESSFUL
    [26.08|15:26:15] JOB dist_2.80 STARTED
    [26.08|15:26:15] JOB dist_2.80 RUNNING
    [26.08|15:26:16] JOB dist_2.80 FINISHED
    [26.08|15:26:16] JOB dist_2.80 SUCCESSFUL
    [26.08|15:26:16] JOB dist_3.00 STARTED
    [26.08|15:26:16] JOB dist_3.00 RUNNING
    [26.08|15:26:18] JOB dist_3.00 FINISHED
    [26.08|15:26:18] JOB dist_3.00 SUCCESSFUL
    [26.08|15:26:18] JOB dist_3.20 STARTED
    [26.08|15:26:18] JOB dist_3.20 RUNNING
    [26.08|15:26:20] JOB dist_3.20 FINISHED
    [26.08|15:26:20] JOB dist_3.20 SUCCESSFUL
    [26.08|15:26:20] JOB dist_3.40 STARTED
    [26.08|15:26:20] JOB dist_3.40 RUNNING
    [26.08|15:26:21] JOB dist_3.40 FINISHED
    [26.08|15:26:21] JOB dist_3.40 SUCCESSFUL
    [26.08|15:26:21] JOB dist_3.60 STARTED
    [26.08|15:26:21] JOB dist_3.60 RUNNING
    [26.08|15:26:23] JOB dist_3.60 FINISHED
    [26.08|15:26:23] JOB dist_3.60 SUCCESSFUL
    [26.08|15:26:23] JOB dist_3.80 STARTED
    [26.08|15:26:23] JOB dist_3.80 RUNNING
    [26.08|15:26:24] JOB dist_3.80 FINISHED
    [26.08|15:26:24] JOB dist_3.80 SUCCESSFUL
    [26.08|15:26:24] JOB dist_4.00 STARTED
    [26.08|15:26:24] JOB dist_4.00 RUNNING
    [26.08|15:26:26] JOB dist_4.00 FINISHED
    [26.08|15:26:26] JOB dist_4.00 SUCCESSFUL

=======
.. parsed-literal::

    [04.02|17:48:15] JOB dist_2.20 STARTED
    [04.02|17:48:15] JOB dist_2.20 RUNNING
    [04.02|17:48:19] JOB dist_2.20 FINISHED
    [04.02|17:48:19] JOB dist_2.20 SUCCESSFUL
    [04.02|17:48:19] JOB dist_2.40 STARTED
    [04.02|17:48:19] JOB dist_2.40 RUNNING
    [04.02|17:48:22] JOB dist_2.40 FINISHED
    [04.02|17:48:22] JOB dist_2.40 SUCCESSFUL
    [04.02|17:48:22] JOB dist_2.60 STARTED
    [04.02|17:48:22] JOB dist_2.60 RUNNING
    [04.02|17:48:25] JOB dist_2.60 FINISHED
    [04.02|17:48:25] JOB dist_2.60 SUCCESSFUL
    [04.02|17:48:25] JOB dist_2.80 STARTED
    [04.02|17:48:25] JOB dist_2.80 RUNNING
    [04.02|17:48:29] JOB dist_2.80 FINISHED
    [04.02|17:48:29] JOB dist_2.80 SUCCESSFUL
    [04.02|17:48:29] JOB dist_3.00 STARTED
    [04.02|17:48:29] JOB dist_3.00 RUNNING
    [04.02|17:48:31] JOB dist_3.00 FINISHED
    [04.02|17:48:32] JOB dist_3.00 SUCCESSFUL
    [04.02|17:48:32] JOB dist_3.20 STARTED
    [04.02|17:48:32] JOB dist_3.20 RUNNING
    [04.02|17:48:34] JOB dist_3.20 FINISHED
    [04.02|17:48:34] JOB dist_3.20 SUCCESSFUL
    [04.02|17:48:34] JOB dist_3.40 STARTED
    [04.02|17:48:34] JOB dist_3.40 RUNNING
    [04.02|17:48:37] JOB dist_3.40 FINISHED
    [04.02|17:48:37] JOB dist_3.40 SUCCESSFUL
    [04.02|17:48:37] JOB dist_3.60 STARTED
    [04.02|17:48:37] JOB dist_3.60 RUNNING
    [04.02|17:48:40] JOB dist_3.60 FINISHED
    [04.02|17:48:40] JOB dist_3.60 SUCCESSFUL
    [04.02|17:48:40] JOB dist_3.80 STARTED
    [04.02|17:48:40] JOB dist_3.80 RUNNING
    [04.02|17:48:43] JOB dist_3.80 FINISHED
    [04.02|17:48:43] JOB dist_3.80 SUCCESSFUL
    [04.02|17:48:43] JOB dist_4.00 STARTED
    [04.02|17:48:43] JOB dist_4.00 RUNNING
    [04.02|17:48:47] JOB dist_4.00 FINISHED
    [04.02|17:48:47] JOB dist_4.00 SUCCESSFUL

>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

Results
~~~~~~~

Print table of results of the distance against the calculated energy.

.. code:: ipython3

<<<<<<< HEAD
   print("== Results ==")
   print("d[A]    E[kcal/mol]")
   for d, e in zip(distances, energies):
       print(f"{d:.2f}    {e:.3f}")
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
    print("== Results ==")
    print("d[A]    E[kcal/mol]")
    for d, e in zip(distances, energies):
        print(f"{d:.2f}    {e:.3f}")
=======
    print("== Results ==")
    try:
        # For AMS2025+ can use JobAnalysis class to perform results analysis
        from scm.plams import JobAnalysis
    
        ja = (
            JobAnalysis(jobs=jobs, std_fields=None)
            .add_field("Dist", lambda j: j.molecule[2].x, display_name="d[A]", fmt=".2f")
            .add_field("Energy", lambda j: j.results.get_energy(unit="kcal/mol"), display_name="E[kcal/mol]", fmt=".3f")
        )
    
        # Pretty-print if running in a notebook
        if "ipykernel" in sys.modules:
            ja.display_table()
        else:
            print(ja.to_table())
    
        energies = ja.Energy
    
    except ImportError:
    
        energies = [j.results.get_energy(unit="kcal/mol") for j in jobs]
    
        print("d[A]    E[kcal/mol]")
        for d, e in zip(distances, energies):
            print(f"{d:.2f}    {e:.3f}")
>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

::

<<<<<<< HEAD
   == Results ==
   d[A]    E[kcal/mol]
   2.20    0.230
   2.40    -0.054
   2.60    -0.127
   2.80    -0.122
   3.00    -0.094
   3.20    -0.066
   3.40    -0.045
   3.60    -0.030
   3.80    -0.020
   4.00    -0.013
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)
.. parsed-literal::

    == Results ==
    d[A]    E[kcal/mol]
    2.20    0.230
    2.40    -0.054
    2.60    -0.127
    2.80    -0.122
    3.00    -0.094
    3.20    -0.066
    3.40    -0.045
    3.60    -0.030
    3.80    -0.020
    4.00    -0.013

=======
.. parsed-literal::

    == Results ==



.. raw:: html

    <div style="max-width: 100%; overflow-x: auto;">
    <table border="1" style="border-collapse: collapse; width: auto; ">
    <thead><tr><th>d[A]<th>E[kcal/mol]</th></tr></thead>
    <tbody>
    <tr><td>2.20</td><td>0.230      </td></tr>
    <tr><td>2.40</td><td>-0.054     </td></tr>
    <tr><td>2.60</td><td>-0.127     </td></tr>
    <tr><td>2.80</td><td>-0.122     </td></tr>
    <tr><td>3.00</td><td>-0.094     </td></tr>
    <tr><td>3.20</td><td>-0.066     </td></tr>
    <tr><td>3.40</td><td>-0.045     </td></tr>
    <tr><td>3.60</td><td>-0.030     </td></tr>
    <tr><td>3.80</td><td>-0.020     </td></tr>
    <tr><td>4.00</td><td>-0.013     </td></tr>
    </tbody>
    </table>
    </div>

>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)

.. code:: ipython3

   import matplotlib.pyplot as plt

   fig, ax = plt.subplots(figsize=(3, 3))
   ax.plot(distances, energies, ".-")
   ax.set_xlabel("He-He distance (Å)")
   ax.set_ylabel("Energy (kcal/mol)");

<<<<<<< HEAD
.. figure:: He2DissociationCurve_files/He2DissociationCurve_12_0.png
||||||| parent of dd6913c (Update some existing examples to use the job analysis tool SO107)

.. image:: He2DissociationCurve_files/He2DissociationCurve_12_0.png

=======

.. image:: He2DissociationCurve_files/He2DissociationCurve_12_0.png


>>>>>>> dd6913c (Update some existing examples to use the job analysis tool SO107)
