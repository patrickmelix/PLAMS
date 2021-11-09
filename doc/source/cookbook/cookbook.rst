Cookbook
========

This is a collection of code snippets showing how to perform recurrent PLAMS tasks.


Settings and input
******************


Create an input block with an header
------------------------------------

This settings

.. code-block:: python

    sett = Settings()
    sett.input.ams.SomeInputBlock['_h'] = 'MyHeader'
    sett.input.ams.SomeInputBlock.SomeOption = 2

will generate the following text input when used as setting for an |AMSJob|:

::

    SomeInputBlock MyHeader
        SomeOption 2
    End


Create an empty input block
---------------------------

This settings

.. code-block:: python

    sett = Settings()
    sett.input.ams.SomeInputBlock

will generate the following text input when used as setting for an |AMSJob|:

::

    SomeInputBlock
    End


Convert an AMS text input into settings object
----------------------------------------------

.. code-block:: python

    from scm.input_parser import InputParser

    text_input="""
    Task SinglePoint
    Engine Band
        Basis
            Type DZ
            Core None
        End
    EndEngine
    """

    s = Settings()
    with InputParser() as parser:
        s.input  = parser.to_settings('ams', text_input)


Molecules
*********

Generate a molecule from a SMILES string
----------------------------------------

.. code-block:: python

    # Compute 10 conformers, optimize with UFF and pick the lowest in energy.
    ethane = from_smiles('C-C', nconfs=10, forcefield='uff')[0]


Extracting Results
******************

You can use the following snippets to retrieve results after running the required calculations:

Directly from Functions
-----------------------

Results can be either red from previous calculations (see :ref:`accessing_old_jobs`) or from an AMSResults instance of a computation just executed within the same workflow.
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

.. _accessing_old_jobs:

Accessing Old Jobs
******************

The following illustrate how to load data from previously executed jobs:

Binding Native PLAMS Jobs
-------------------------

.. warning::
   The jobs should be loaded with a version of PLAMS that is consistent with the version originally used to run the jobs.


From an existing PLAMS working directory with the contents

::

   OLDDIR/
   ├── OLDJOB1/
   |   ├── ams.log
   |   ├── ams.rkf
   |   ├── OLDJOB1.dill
   |   ├── OLDJOB1.err
   |   ├── OLDJOB1.in
   |   ├── OLDJOB1.out
   |   ├── OLDJOB1.run
   |   ├── engine.rkf
   |   ├── output.xyz
   ├── input
   └── logfile

we can bind an instance of the AMSJob class by making use of the `.dill` file.
The AMSJob object in turn contains a results object, which gives access to the data previously calculated.
This can be achieved with the following snippet::

   path       = "OLDDIR/OLDJOB1/OLDJOB1.dill"
   single_JOB = load(path)                                       # AMSJob instance
   if single_JOB.ok():
      energy     = single_JOB.results.get_energy()               # load the desired properties
      structure  = single_JOB.results.get_main_molecule()
      propertyX  = single_JOB.results.readrkf('AMSResults', 'DipoleMoment', file='engine')

More often than not, the working directory will include multiple individual subdirectories, each containing individual PLAMS job.

::

   OLDDIR/
   ├── OLDJOB1/
   |   ├── ams.log
   |   ├── ams.rkf
   |   ├── OLDJOB1.dill
   |   ├── OLDJOB1.err
   |   ├── OLDJOB1.in
   |   ├── OLDJOB1.out
   |   ├── OLDJOB1.run
   |   ├── engine.rkf
   |   ├── output.xyz
   ├── OLDJOB2/
   |   ├── ams.log
   |   ├── ams.rkf
   |   ├── OLDJOB2.dill
   |   ├── OLDJOB2.err
   |   ├── OLDJOB2.in
   |   ├── OLDJOB2.out
   |   ├── OLDJOB2.run
   |   ├── engine.rkf
   |   ├── output.xyz
   ├── OLDJOB3/
   |   ├── ams.log
   |   ├── ams.rkf
   |   ├── OLDJOB3.dill
   |   ├── OLDJOB3.err
   |   ├── OLDJOB3.in
   |   ├── OLDJOB3.out
   |   ├── OLDJOB3.run
   |   ├── engine.rkf
   |   ├── output.xyz
   ├── input
   └── logfile

These can be loaded using the `load_all` function and by providing only the path to the top-level directory::

   path       = "OLDDIR"
   all_JOBS   = load_all(path)

Note that `load_all` wraps the `load` function used above and therefore requires existing `.dill` files in each of the loaded subdirectories.
The `load_all` function yields a dictionary with the paths of the `.dill` files as keys and the corresponding job object as values::

   print(all_JOBS)

::

   {'/home/user/OLDDIR/OLDJOB1/OLDJOB1.dill': <scm.plams.interfaces.adfsuite.ams.AMSJob object at 0x7f0baad340b8>,
    '/home/user/OLDDIR/OLDJOB2/OLDJOB2.dill': <scm.plams.interfaces.adfsuite.ams.AMSJob object at 0x7f0baacf24a8>,
    '/home/user/OLDDIR/OLDJOB3/OLDJOB3.dill': <scm.plams.interfaces.adfsuite.ams.AMSJob object at 0x7f0baad06cf8>}

We can now access these AMSJob instances::

   for this_JOB in all_JOBS.values():
      if this_JOB.ok():
         energy     = this_JOB.results.get_energy()
         structure  = this_JOB.results.get_main_molecule()
         propertyX  = this_JOB.results.readrkf('AMSResults', 'DipoleMoment', file='engine')


Binding old RKF Files
---------------------
In cases where the `.dill` files are not available any more, it is still possible to load the contents of previously generated `.rkf` files into a PLAMS workflow::

   path       = "OLDDIR/OLDJOB1/"
   ext_JOB    = AMSJob.load_external(path)
   if ext_JOB.ok():
      energy     = ext_JOB.results.get_energy()
      structure  = ext_JOB.results.get_main_molecule()

If the `.rkf` file does originate from some other source than any of the direct AMS engines, also an instance of the more generic `SingleJob` class can be used::

   path       = "OLDDIR/OLDJOB1/ams.rkf"
   ext_JOB    = SingleJob.load_external(path)

The downside of this latter approach is that the accessibility to the data is very limited and has to be implemented mostly in terms of pattern-matching searches in the output files.

An alternative way is to make use of the `KFReader` class::

   path       = "OLDDIR/OLDJOB1/ams.rkf"
   rkf_reader = KFReader(path)
   n_steps    = rkf_reader.read("History", "nEntries")
   energy     = rkf_reader.read("History", "Energy({})".format(n_steps))
   structure  = rkf_reader.read("History", "Coords({})".format(n_steps))

Note that also the KFReader class lacks most of the shortcut functions of a proper `AMSResults` object so that the access to the data has to be specified manually.




