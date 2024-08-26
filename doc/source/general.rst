.. _General:

General
============

This section contains general information about installing and running PLAMS.

For a quick-start guide see the :ref:`examples`.

What is PLAMS
-------------

PLAMS (Python Library for Automating Molecular Simulation) is a collection of tools that aims to provide powerful, flexible and easily extendable Python interface to molecular modeling programs.
It takes care of input preparation, job execution, file management and output processing as well as helps with building more advanced data workflows.

Usually the daily work of a computational chemist consists of running a number of calculations.
Those calculations are done using one or more molecular modeling programs like ADF, BAND, Turbomole or Dirac (we will call such programs *external binaries*).
Obtaining results with one of such programs requires a series of steps.
First, the subject of the problem (description of a molecular system, set of desired simulation parameters) has to be presented in the format understandable by molecular modeling program and written to an *input file* which is usually a text file.
Then the program is executed, it runs and produces *output* which is a collection of text or binary files.
That output usually contains more information than is required for a particular problem so data of interest has to be extracted and (possibly) postprocessed.
That different computational tools use different input and output formats and are executed differently.
In most cases many *single calculations* need to be performed to solve the problem of interest.
That requires significant effort to be put into data hygiene to avoid confusing or overwriting input and output files from distinct calculations.

Each of the above steps, apart from actual calculation done by a molecular modeling program, needs to be performed by a human.
Preparing and editing text files, creating folders in the filesystem, copying files between them and reading data from output are tedious, repetitive and highly error-prone work.
Some users deal with it using automation, usually in form of ad hoc shell scripts.
A few programs, like Amsterdam Modeling Suite, offer graphical user interface to help with this kind of work, but again, input preparation and output examination, even though assisted with convenient tools, have to be done by a human.
Quite often it turns out to be a performance bottleneck to create big  automatic computational workflows, where output data of one calculation is used (usually after some processing) as an input to another calculation, sometimes done with different program on a different machine.

PLAMS was created to solve these problems.

* It takes responsibility of tiresome and monotonous technical details allowing you to focus on real science and your problem.
* It lets you do all the things mentioned above (and many more) using simple Python scripts.
* It gives you a helping hand with automating repetitive or complicated tasks while still leaving you with 100% control over what is really happening with your files, disks and CPUs.


What can be done with PLAMS
---------------------------

The key design principle of PLAMS is *flexibility*.
If something (and by something we mean: adjusting an input parameter, executing some program with particular options, extracting a value from output etc.) can be done by hand, it can be done with PLAMS.
The internal structure of the library was designed in highly modular, object-oriented manner.
Thanks to that it takes very little effort to adjust its behavior to one's personal needs or to extend its functionality.


The most important features of PLAMS:

*   preparing, running and examining results of a molecular modeling jobs from within a single Python script
*   convenient automatic file and folder management
*   running jobs in parallel without a need to prepare a special parallel script
*   integration with popular job schedulers (OGE, SLURM, TORQUE)
*   molecular coordinates manipulation using user-friendly |Molecule| class, supporting various formats (``xyz``, ``mol``, ``mol2``, ``pdb``), as well as interfaces to ASE and RDKit.
*   prevention of multiple runs of the same job
*   easy data transfer between separate runs
*   efficient restarting in case of crash
*   full coverage of all input options and output data in :ref:`Amsterdam Modeling Suite <AMSSuite>` programs
*   easy extendable for other programs, job schedulers, file formats etc.
*   interfaces for Dirac, ORCA, CP2K, DFTB+, Crystal, and more.


What PLAMS is *not*
-------------------------

It should be stressed here that PLAMS is not a *program*, it's a *library*.
That means it's not a standalone tool, it doesn't run or do anything by itself.
To work properly, it needs both an external binary on one side and a properly written Python script on the other.
Being a library means that PLAMS is in fact just a collection of commands and objects that can be used from within a regular Python script to perform common molecular modeling tasks.

Because of the above, PLAMS won't take your hand and guide you, it won't detect and warn you if you are about to do something stupid and it won't do anything except the things you explicitly asked for.
You have to understand what you are doing, you have to know how to use the binary you want PLAMS to work with and you have to have at least some basic knowledge of Python programming language.


Contributing
-----------------

Users are warmly welcome to help with enriching PLAMS, as well as to provide
any kind of feedback regarding either PLAMS itself or this documentation to
support@scm.com or directly on PLAMS `GitHub page
<https://github.com/SCM-NV/PLAMS>`_.


Library contents
-------------------------

PLAMS comes in a form of Python 3 package (earlier versions were also compatible with Python 2, but as Python 3 becomes more and more popular, we decided to drop Python 2 compatibility).
The root folder of the package contains the following:

*   ``core``: core subpackage containing all the essential modules defining the skeleton of the library
*   ``interfaces``: subpackage with interfaces to external binaries
*   ``tools``: subpackage with small utilities like unit converter, periodic table, file readers etc.
*   ``recipes``: subpackage with examples of simple job types built on top basic PLAMS elements
*   ``doc``: Sphinx source of this documentation


Installing PLAMS
-------------------------

You can install PLAMS on your computer using one of the following ways:

1.  If you are using Amsterdam Modeling Suite, PLAMS is included as a part of ``scm`` Python package (``$AMSHOME/scripting/scm/plams``) and configured to work with the `Python Stack <../Scripting/Python_Stack/Python_Stack.html>`__ coming with AMSuite (you can access it with ``$AMSBIN/amspython`` command).

2.  The latest PLAMS stable release can be installed directly from PyPi by typing ``pip install plams`` in your command line.

3.  Any current or historic version can be downloaded or cloned from PLAMS `GitHub page <https://github.com/SCM-NV/PLAMS>`_.
    The ``release`` branch points to the latest stable release, while the ``trunk`` branch is the most recent development snapshot.

4.  You can combine methods 2 and 3 and fetch PLAMS from GitHub using ``pip``: ``pip install git+https://github.com/SCM-NV/PLAMS.git@master`` (make sure to have Git installed and to choose the proper branch)

PLAMS requires the following Python packages as dependencies:

*   `numpy <http://www.numpy.org>`_
*   `dill <https://pypi.python.org/pypi/dill>`_ (enhanced pickling)
*   `ase <https://wiki.fysik.dtu.dk/ase>`_ (optional dependency)
*   `rdkit <https://pypi.org/project/rdkit>`_ (optional dependency)

If you are using Amsterdam Modeling Suite, all the above packages are already included in our Python stack.
When you install PLAMS using ``pip``, the required packages (numpy and dill) will be installed automatically.
For optional dependencies, or in any other case you can install them with ``pip install [package name]``.



Running PLAMS
-------------------------

Inside your Python interpreter or in Python scripts, PLAMS is visible as a subpackage of the ``scm`` package, so you can import it with one of the following commands

.. code-block:: python

    # myscript.py
    import scm.plams
    from scm import plams
    from scm.plams import Settings, Molecule, Atom, AMSJob  #  ... other required components


The script can then be run using a Python interpreter e.g. ``$AMSBIN/amspython myscript.py``.


.. _plams-defaults:

PLAMS defaults
-------------------------

PLAMS has a global ``config`` object which contains all the configuration settings for the PLAMS script.
On startup, these are initialized to default values, the details of which are explained in the description of each property on ``config`` (see |global-settings|).
It is recommended to have a look at that these options, to give an overview of what behaviour can be configured.

To change a setting for a script, just set the relevant option on the config to the preferred value, after the import statements.
For example:

.. code-block:: python

    config.log.stdout = 1
    config.job.pickle = False
    config.default_jobrunner = JobRunner(parallel=True, maxjobs=8)


.. _working-folder:

Working folder location
~~~~~~~~~~~~~~~~~~~~~~~~~

All files produced by PLAMS and other programs executed by it are saved in the main working folder (usually in some of its subfolders).
Each separate run of PLAMS has a separate main working folder.
By default the main working folder is located in the directory where your script was executed and is called ``plams_workdir`` (``plams_workdir.002`` if ``plams_workdir`` already existed).
The name and location for the main working folder can be altered by calling the |init| function.
For example:


.. code-block:: python

    init(path="my/path", folder="my_folder")


.. note::

    Each PLAMS run creates a fresh, empty directory for its main working folder.
    If you try to use an existing folder (or don't pick any and ``plams_workdir`` already exists in the current directory), a unique folder is going to be created anyway, by appending ``.002`` (or ``.003``, ``.004`` and so on) to the name of your folder.

What's new in PLAMS for AMS2025?
--------------------------------------

* Call to |init| at the start of a PLAMS script is no longer required. In addition, a call to |finish| is automatically registered at exit. For more information see :ref:`public-functions`.
* Deprecated ``BandJob``, ``DFTBJob``, ``UFFJob``, ``MOPACJob``, ``ReaxFFJob`` and ``ADFJob`` jobs have been removed. These were deprecated since AMS2019, and replaced by |AMSJob|.

What's new in PLAMS for AMS2024?
--------------------------------------

* :ref:`Packmol interface <PackmolInterface>` has been extended to pack in crystal voids and to get the total system charge from the sum of the constituent molecules

* Additions to |AMSResults|: get_normal_modes(), get_polarizability(), get_ir_spectrum(), get_ir_spectrum_md(), get_frequency_spectrum(), get_force_constants()

* Additions to |Molecule|: get_moments_of_inertia(), get_gyration_radius(), align2mol()

What's new in PLAMS for AMS2023?
--------------------------------------

* The :ref:`AMSCalculator` class for running any AMS engine with ASE (see: :ref:`ASECalculatorExample`)

* Classes for calculating :ref:`reduction and oxidation potentials  <RedoxExample>` with ADF and optionally COSMO-RS

* The :ref:`ADFCOSMORSCompoundJob <ADFCOSMORSCompound>` class for running jobs equivalent to "Task COSMO-RS Compound" in the AMS GUI. Such a job generates a .coskf file for use with COSMO-RS.

* The calculation of the :ref:`vibronic density of states<fcf_dos>` has been added to PLAMS.

* Classes for running and restarting :ref:`molecular dynamics (MD) jobs with AMS <AMSMDJob>`

* A class for generating and analyzing :ref:`conformers <conformers_interface>`

* :ref:`Quick jobs <Quickjobs>`, like for example the ``preoptimize()`` function let you quickly optimize a Molecule

* :ref:`Packmol interface <PackmolInterface>` for generating liquid and gas mixtures, solid-liquid interfaces, and microsolvation spheres

* :ref:`FileFormatConversionTools` for converting VASP, Gaussian, or Quantum ESPRESSO output to ams.rkf and engine.rkf files that can be opened with the AMS GUI

* :ref:`PlottingTools` for plotting a molecule or ASE Atoms inside a Jupyter notebook

* :ref:`PlottingTools` for plotting the :ref:`electronic band structure <BandStructureExample>`

* Additions to |AMSResults|: get_homo_energies(), get_lumo_energies, get_smallest_homo_lumo_gap()

* Additions to |Molecule|: guess_atomic_charges(), set_density(), get_unique_bonds(), get_unique_angles()

* Many new :ref:`examples`
