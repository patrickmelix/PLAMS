Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

    import numpy as np
    from scm.plams import Settings, Molecule, Atom, AMSJob

Setup Dimer
~~~~~~~~~~~

Create Helium atoms and an array of interatomic distances at which to
run calculation.

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

For each interatomic distance, create a Helium dimer molecule with the
required geometry then the single point energy calculation job. Run the
job and extract the energy.

.. code:: ipython3

    energies = []
    for d in distances:
        mol = Molecule()
        mol.add_atom(Atom(symbol=atom1, coords=(0.0, 0.0, 0.0)))
        mol.add_atom(Atom(symbol=atom2, coords=(d, 0.0, 0.0)))
        job = AMSJob(molecule=mol, settings=sett, name=f"dist_{d:.2f}")
        job.run()
        energies.append(job.results.get_energy(unit="kcal/mol"))


.. parsed-literal::

    [13.08|11:10:14] JOB dist_2.20 STARTED
    [13.08|11:10:14] JOB dist_2.20 RUNNING
    [13.08|11:10:16] JOB dist_2.20 FINISHED
    [13.08|11:10:16] JOB dist_2.20 SUCCESSFUL
    [13.08|11:10:16] JOB dist_2.40 STARTED
    [13.08|11:10:16] JOB dist_2.40 RUNNING
    [13.08|11:10:18] JOB dist_2.40 FINISHED
    [13.08|11:10:18] JOB dist_2.40 SUCCESSFUL
    [13.08|11:10:18] JOB dist_2.60 STARTED
    [13.08|11:10:18] JOB dist_2.60 RUNNING
    [13.08|11:10:20] JOB dist_2.60 FINISHED
    [13.08|11:10:20] JOB dist_2.60 SUCCESSFUL
    [13.08|11:10:20] JOB dist_2.80 STARTED
    [13.08|11:10:20] JOB dist_2.80 RUNNING
    [13.08|11:10:22] JOB dist_2.80 FINISHED
    [13.08|11:10:22] JOB dist_2.80 SUCCESSFUL
    [13.08|11:10:22] JOB dist_3.00 STARTED
    [13.08|11:10:22] JOB dist_3.00 RUNNING
    [13.08|11:10:24] JOB dist_3.00 FINISHED
    [13.08|11:10:24] JOB dist_3.00 SUCCESSFUL
    [13.08|11:10:24] JOB dist_3.20 STARTED
    [13.08|11:10:24] JOB dist_3.20 RUNNING
    [13.08|11:10:25] JOB dist_3.20 FINISHED
    [13.08|11:10:25] JOB dist_3.20 SUCCESSFUL
    [13.08|11:10:25] JOB dist_3.40 STARTED
    [13.08|11:10:25] JOB dist_3.40 RUNNING
    [13.08|11:10:27] JOB dist_3.40 FINISHED
    [13.08|11:10:27] JOB dist_3.40 SUCCESSFUL
    [13.08|11:10:27] JOB dist_3.60 STARTED
    [13.08|11:10:27] JOB dist_3.60 RUNNING
    [13.08|11:10:29] JOB dist_3.60 FINISHED
    [13.08|11:10:29] JOB dist_3.60 SUCCESSFUL
    [13.08|11:10:29] JOB dist_3.80 STARTED
    [13.08|11:10:29] JOB dist_3.80 RUNNING
    [13.08|11:10:31] JOB dist_3.80 FINISHED
    [13.08|11:10:31] JOB dist_3.80 SUCCESSFUL
    [13.08|11:10:31] JOB dist_4.00 STARTED
    [13.08|11:10:31] JOB dist_4.00 RUNNING
    [13.08|11:10:32] JOB dist_4.00 FINISHED
    [13.08|11:10:32] JOB dist_4.00 SUCCESSFUL


Results
~~~~~~~

Print table of results of the distance against the calculated energy.

.. code:: ipython3

    print("== Results ==")
    print("d[A]    E[kcal/mol]")
    for d, e in zip(distances, energies):
        print(f"{d:.2f}    {e:.3f}")


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

