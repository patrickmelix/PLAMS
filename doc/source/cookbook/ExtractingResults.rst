Extracting Results
==================

You can use the following snippets to retrieve results after running the required calculations:

Directly from Functions
-----------------------

Results can be either red from previous calculations (see `How to load old jobs <LoadOldJobs.html>`_) or from an AMSResults instance of a computation just executed within the same workflow.
In either case an AMSResults object should be present at runtime::

   myAMSJob.run()
   myAMSResults = myAMSJob.results if myAMSJob.ok() else None

.. warning::
   Access to any results data should only occur under the condition that `AMSJob.ok()` indicate a successful termination of the computation

Examples: Total Energy and Final Structure
++++++++++++++++++++++++++++++++++++++++++
Multiple functions of the AMSResults API allow for simple access of the most common results

::

   myAMSEnergy = myAMSResults.get_energy(unit='au')

   myAMSStructure = myAMSResults.get_main_molecule()

AMSResults API Functions
++++++++++++++++++++++++
The following members of an AMSResults instance can be used as shown in the above examples to read results

.. list-table::
   :widths: 25 25 50 100
   :header-rows: 1

   * - Property
     - Function
     - Return Type
     - Details
   * - Structure
     - `get_molecule(section)`
     - `Molecule`
     - Structure from `section`
   * - 
     - `get_input_molecule()`
     - `Molecule`
     - Input structure
   * - 
     - `get_main_molecule()`
     - `Molecule`
     - Final structure from any AMS task
   * -
     - `get_history_molecule(step)`
     - `Molecule`
     - Structure from history section at step # `step`
   * - Energy
     - `get_energy()`
     - `Float`
     - Final energy
   * - Gradients
     - `get_gradients()`
     - `Array` (numpy)
     - Gradients from engine calculation
   * - Stress tensor
     - `get_stresstensor()`
     - `Array` (numpy)
     - Stress tensor from periodic engine calculation
   * - Hessian
     - `get_hessian()`
     - `Array` (numpy)
     - Hessian from frequency calculation (AMS/engine)
   * - Elastic tensor
     - `get_elastictensor()`
     - `Array` (numpy)
     - Elastic tensor from periodic calculation
   * - Frequencies
     - `get_frequencies()`
     - `Array` (numpy)
     - Vibrational frequencies
   * - Atomic Charges
     - `get_charges()`
     - `Array` (numpy)
     - Atomic partial charges
   * - Dipole vector
     - `get_dipolemoment()`
     - `Array` (numpy)
     - Electric dipole moment
   * - Nuclear gradients of dipole vector
     - `get_dipolegradients()`
     - `Array` (numpy)
     - Nuclear Gradients of Electric dipole moment
   
From the RKF Interface
----------------------
Other properties not listed in the table above should be retrieved as follows::

   myProperty = myAMSResults.readrkf(section, variable)

It is the responsibility of the user to provide the correct names for `section` and `variable` under which the required result is stored in the rkf file.

Finding Section/Variable Pairs
------------------------------
Looking up the names of the needed sections and variable within rkf files is typically needed for more intricate properties when writing a new PLAMS workflow.
There are two main approaches to search for this information.

From Python Directories
+++++++++++++++++++++++
The AMSResults member function::

   get_rkf_skeleton()

returns a dictionary containing the available sections as keys and the containing varible names as values

KFBrowser
+++++++++
KFBrowser is a GUI module used to inspect rkf files.

.. rst-class:: steps

   \
     | **1.** Open KFBrowser in the GUI via **SCM → KFBrowser**
     | **2.** By default KFBrowser opens the `ams.rkf` file. Where neccessary, switch to **File → open → <engine>.rkf**
     | **3.** Press **ctrl + e** or select **File → Expert Mode** to display the stored file contents
     | **4.** Find the entry of interest. While this is a sometimes not trivial step, most often the required variable is found in either the ``Properties`` or ``AMSresults`` sections.
     | **5.** Once found, the names for `section` and `variable` listed in the rkf file directly corresponds to the `section`/`variable` pair to be used in the `readrkf` function as shown above. 

.. note::
   When reading results from a different rkf file than `ams.rkf` the filename has to be specified as::

     myEngineProperty = myAMSResults.readrkf(section, variable, file=<engine>)

   whereas `<engine>` corresponds to the file `<engine>.rkf` present in the calculation directory.
