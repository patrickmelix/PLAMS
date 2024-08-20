Worked Example
--------------

Initial Imports
~~~~~~~~~~~~~~~

.. code:: ipython3

    from scm.plams import AMSResults, Units, add_to_class, Settings, read_molecules, AMSJob

Helper Functions
~~~~~~~~~~~~~~~~

Set up a couple of useful functions for extracting results.

.. code:: ipython3

    @add_to_class(AMSResults)
    def get_excitations(results):
        """Returns excitation energies (in eV) and oscillator strenghts (in Debye)."""
        if results.job.ok():
            exci_energies_au = results.readrkf("Excitations SS A", "excenergies", file="engine")
            oscillator_str_au = results.readrkf("Excitations SS A", "oscillator strengths", file="engine")
            # The results are stored in atomic units. Convert them to more convenient units:
            exci_energies = Units.convert(exci_energies_au, "au", "eV")
            oscillator_str = Units.convert(oscillator_str_au, "au", "Debye")
            return exci_energies, oscillator_str
        else:
            return [], []

.. code:: ipython3

    @add_to_class(AMSResults)
    def has_good_excitations(results, min_energy, max_energy, oscillator_str_threshold=1e-4):
        """Returns True if there is at least one excitation with non-vanishing oscillator strenght
        in the energy range [min_energy, max_energy]. Unit for min_energy and max energy: eV."""
        exci_energies, oscillator_str = results.get_excitations()
        for e, o in zip(exci_energies, oscillator_str):
            if min_energy < e < max_energy and o > oscillator_str_threshold:
                return True
        return False

Calculation settings
~~~~~~~~~~~~~~~~~~~~

Configure the settings for the various jobs.

.. code:: ipython3

    # Settings for geometry optimization with the AMS driver:
    go_sett = Settings()
    go_sett.input.ams.Task = "GeometryOptimization"
    go_sett.input.ams.GeometryOptimization.Convergence.Gradients = 1.0e-4

.. code:: ipython3

    # Settings for single point calculation with the AMS driver
    sp_sett = Settings()
    sp_sett.input.ams.Task = "SinglePoint"

.. code:: ipython3

    # Settings for the DFTB engine (including excitations)
    dftb_sett = Settings()
    dftb_sett.input.dftb.Model = "SCC-DFTB"
    dftb_sett.input.dftb.ResourcesDir = "QUASINANO2015"
    dftb_sett.input.dftb.Properties.Excitations.TDDFTB.calc = "singlet"
    dftb_sett.input.dftb.Properties.Excitations.TDDFTB.lowest = 10
    dftb_sett.input.dftb.Occupation.Temperature = 5.0

.. code:: ipython3

    # Settings for the geometry optimization with the ADF engine
    adf_sett = Settings()
    adf_sett.input.adf.Basis.Type = "DZP"
    adf_sett.input.adf.NumericalQuality = "Basic"

.. code:: ipython3

    # Settings for the excitation calculation using the ADF engine
    adf_exci_sett = Settings()
    adf_exci_sett.input.adf.Basis.Type = "TZP"
    adf_exci_sett.input.adf.XC.GGA = "PBE"
    adf_exci_sett.input.adf.NumericalQuality = "Basic"
    adf_exci_sett.input.adf.Symmetry = "NoSym"
    adf_exci_sett.input.adf.Excitations.lowest = 10
    adf_exci_sett.input.adf.Excitations.OnlySing = ""

Load Molecules
~~~~~~~~~~~~~~

Import all xyz files in the folder ‘molecules’.

.. code:: ipython3

    molecules = read_molecules("molecules")

DFTB Prescreen
~~~~~~~~~~~~~~

Perform an initial prescreen of all molecules with DFTB.

.. code:: ipython3

    promising_molecules = {}

.. code:: ipython3

    for name, mol in molecules.items():
        dftb_job = AMSJob(name="DFTB_" + name, molecule=mol, settings=go_sett + dftb_sett)
        dftb_job.run()
    
        if dftb_job.results.has_good_excitations(1, 6):
            promising_molecules[name] = dftb_job.results.get_main_molecule()


.. parsed-literal::

    [13.08|15:32:02] JOB DFTB_H2O STARTED
    [13.08|15:32:02] JOB DFTB_H2O RUNNING
    [13.08|15:32:03] JOB DFTB_H2O FINISHED
    [13.08|15:32:03] JOB DFTB_H2O SUCCESSFUL
    [13.08|15:32:03] JOB DFTB_NH3 STARTED
    [13.08|15:32:03] JOB DFTB_NH3 RUNNING
    [13.08|15:32:11] WARNING: Job DFTB_NH3 finished with nonzero return code
    [13.08|15:32:11] JOB DFTB_NH3 CRASHED
    [13.08|15:32:11] JOB DFTB_S2Cl2 STARTED
    [13.08|15:32:11] JOB DFTB_S2Cl2 RUNNING
    [13.08|15:32:12] JOB DFTB_S2Cl2 FINISHED
    [13.08|15:32:12] JOB DFTB_S2Cl2 SUCCESSFUL


::


    ---------------------------------------------------------------------------

    NameError                                 Traceback (most recent call last)

    Cell In [10], line 6
          3 dftb_job.run()
          5 if dftb_job.results.has_good_excitations(1, 6):
    ----> 6     promising_molecules[name] = dftb_job.results.get_main_molecule()


    NameError: name 'promising_molecules' is not defined


.. code:: ipython3

    print(f"Found {len(promising_molecules)} promising molecules with DFTB")

Optimization and excitations calculation with ADF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each of the molecules identified in the prescreen, run a further
calculation with ADF.

.. code:: ipython3

    for name, mol in promising_molecules.items():
        adf_go_job = AMSJob(name="ADF_GO_" + name, molecule=mol, settings=go_sett + adf_sett)
        adf_go_job.run()
    
        optimized_mol = adf_go_job.results.get_main_molecule()
    
        adf_exci_job = AMSJob(name="ADF_exci_" + name, molecule=optimized_mol, settings=sp_sett + adf_exci_sett)
        adf_exci_job.run()
    
        if adf_exci_job.results.has_good_excitations(2, 4):
            print(f"Molecule {name} has excitation(s) satysfying our criteria!")
            print(optimized_mol)
            exci_energies, oscillator_str = adf_exci_job.results.get_excitations()
            print("Excitation energy [eV], oscillator strength:")
            for e, o in zip(exci_energies, oscillator_str):
                print(f"{e:8.4f}, {o:8.4f}")
