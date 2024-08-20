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
*   ``scripts``: folder with executable scripts

An imporant part of PLAMS worth mentioning here is the executable script used for running or restarting your workflows.
It is called simply ``plams`` and it's located in the ``scripts`` folder.
We will refer to it as the *launch script* or simply *launcher*.
Further in this section you can find the dedicated chapter explaining the usage of the launch script.



Installing PLAMS
-------------------------

You can install PLAMS on your computer using one of the following ways:

1.  If you are using Amsterdam Modeling Suite, PLAMS is included as a part of ``scm`` Python package (``$AMSHOME/scripting/scm/plams``) and configured to work with a built-in Python coming with AMSuite (you can access it with ``amspython`` command).
    The launch script is added to ``$AMSBIN``, so it should be directly visible from your command line (as long as ``$AMSBIN`` is in your ``$PATH``).

2.  The latest PLAMS stable release can be installed directly from PyPi by typing ``pip install plams`` in your command line.
    The launch scipt will be installed along other global system executables (platform dependent) and should be visible from your command line.

3.  Any current or historic version can be downloaded or cloned from PLAMS `GitHub page <https://github.com/SCM-NV/PLAMS>`_.
    The ``release`` branch points to the latests stable release, while the ``master`` branch is the most recent development snapshot.

4.  You can combine methods 2 and 3 and fetch PLAMS from GitHub using ``pip``: ``pip install git+https://github.com/SCM-NV/PLAMS.git@master`` (make sure to have Git installed and to choose the proper branch)

PLAMS requires the following Python packages as dependencies:

*   `numpy <http://www.numpy.org>`_
*   `dill <https://pypi.python.org/pypi/dill>`_ (enhanced pickling)
*   `ase <https://wiki.fysik.dtu.dk/ase>`_ (optional dependency)
*   `rdkit <https://pypi.org/project/rdkit>`_ (optional dependency)

If you are using Amsterdam Modeling Suite, all the above packages are already included in our Python stack.
When you install PLAMS using ``pip``, the required packages (numpy and dill)will be installed automatically.
For optional dependencies, or in any other case you can install them with ``pip install [package name]``.



Running PLAMS
-------------------------

Inside your Python interpreter or in Python scripts PLAMS is visible as a subpackage of the ``scm`` package, so you can import it with one of the following commands

.. code-block:: python

    import scm.plams
    from scm import plams
    from scm.plams import *

.. note::

    Usually in Python ``import *`` is considered a bad practice and discouraged.
    However, PLAMS internally takes care of the namespace tidiness and imports only necessary things with ``import *``.
    Importing with ``import *`` allows you to use identifiers like ``Molecule`` or ``AMSJob`` instead of ``scm.plams.Molecule`` or ``scm.plams.AMSJob`` which makes your scripts shorter and more readable.
    Throughout this documentation it is assumed that ``import *`` is used so identifiers are not prefixed with ``scm.plams.`` in any example.

A PLAMS script is in fact a general Python script that makes use of classes and functions defined in the PLAMS library.
Of course PLAMS can be also run interactively.
After starting your favorite Python interpreter you need to manually import and initialize the environment with ``from scm.plams import *``.
Then you can interactively run any Python command relying on PLAMS.


.. _plams-defaults:

PLAMS Defaults
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



.. _master-script:

The launch script
-------------------------

The launch script is an executable file called simply ``plams`` located in the ``scripts`` folder.
If the ``$PATH`` variable is configured properly, you can type in your command line ``plams -h`` or ``plams --help`` for a short help message.

The launch script provides another way of executing PLAMS scripts.
Note that the launch script handles the PLAMS imports so this does not need to be done in the the script itself.
The launch script offers several command line arguments allowing you to tune the behavior of your script without a need to edit the script itself.


Working folder location
~~~~~~~~~~~~~~~~~~~~~~~~~

The launch script allows you to pick custom name and location for the main working folder.
All files produced by PLAMS and other programs executed by it are saved in the main working folder (usually in some of its subfolders).
Each separate run of PLAMS has a separate main working folder.

By default the main working folder is located in the directory where your script was executed and is called ``plams_workdir`` (``plams_workdir.002`` if ``plams_workdir`` already existed).
The name and location for the main working folder can be altered by calling the |init| function.
For example:


.. code-block:: python

    init(path="my/path", folder="my_folder")


The launch script also allows you to pick custom name and location for the main working folder.
You can change that by supplying ``-p`` and ``-f`` (or ``--path`` and ``--folder``) arguments to the launcher to choose the location and the name of the main working folder.
For example the command::

    plams -p /home/user/science -f polymers myscript.plms


will use ``/home/user/science/polymers`` as the main working folder regardless where this command was executed.

.. note::

    Each PLAMS run creates a fresh, empty directory for its main working folder.
    If you try to use an existing folder (or don't pick any and ``plams_workdir`` already exists in the current directory), a unique folder is going to be created anyway, by appending ``.002`` (or ``.003``, ``.004`` and so on) to the name of your folder.


Passing variables
~~~~~~~~~~~~~~~~~~~~~~~~~

When using the launcher you can pass variables to your script directly from the command line.
This can be done with ``-v`` (or ``--var``) parameter that follows the syntax ``-v variable=value`` (mind the lack of spaces around equal sign, it is a must).
For a script executed that way, there is an additional global string variable with the name ``variable`` and the value ``'value'`` visible from within the script.
For example if the script in file ``script1.plms`` looks like this::

    print('Chosen basis: ' + basis)
    print('Number of points: ' + n)
    print(type(n))
    # do something depending on n and basis

and you execute it with::

    plams -v n=10 -v basis=DZP script1.plms

the standard output will be:

.. code-block:: none

    Chosen basis: DZP
    Number of points: 10
    str
    [rest of the output]

Three important things to keep in mind about ``-v`` parameter:

*   no spaces around equal sign,
*   each variable requires separate ``-v``,
*   the type of the variable is **always** string (like in the example above).
    If you want to pass some numerical values, make sure to convert them from strings to numbers inside your script.


Importing past jobs
~~~~~~~~~~~~~~~~~~~~~~~~~

You can instruct the launcher to load the results of some previously run jobs by supplying the path to the main working folder of a finished PLAMS run with ``-l`` (or ``--load``) parameter.
To find out why this could be useful, please see |pickling| and |RPM|.

This mechanism is equivalent to using |load_all| function at the beginning of your script.
That means executing your script with ``plams -l /some/path myscript.plms`` works just like putting ``load_all('/some/path')`` at the beginning of ``myscript.plms`` and running it with ``plams myscript.plms``.
The only difference is that, when using |load_all| inside the script, you can access each of the loaded jobs separately by using the dictionary returned by |load_all|.
This is not possible with ``-l`` parameter, but all the loaded jobs will be visible to |RPM|.

Multiple different folders can be supplied with ``-l`` parameter, but each of them requires a separate ``-l`` flag::

    plams -l /some/path -l /other/path myscript.plms


Restarting failed script
~~~~~~~~~~~~~~~~~~~~~~~~~

The launch script can be called with an additional argumentless ``-r`` parameter (or ``--restart``).
In such a case the launcher enters "restart mode".
In the restart mode the folder specified by ``-f`` (or the latest ``plams_workdir[.xxx]`` if ``-f`` is not used) is first renamed by appending ``.res`` to folder's original name (let's call it ``foldername``).
Successful jobs from ``foldername.res`` are loaded at the beginning of the current run, which is executed in a new, empty main working folder called ``foldername``.
Whenever the new run encounters a job identical to a successful job present in ``foldername.res``, the new job execution is skipped and the whole job folder is linked (hardlinked) from ``foldername.res`` to ``foldername``.
That way the restart run will not redo any work present in old ``foldername``, but rather back it up to ``foldername.res`` and restart from the point when the old run was terminated.
For example, after::

    $ plams -f stuff myscript.plms
    [17:28:40] PLAMS working folder: /home/user/stuff
    # [some successful work]
    [17:56:22] Execution interrupted by the following exception:
    # [exception details]

you can edit ``myscript.plms``, remove the cause of crash and restart your script with::

    $ plams -f stuff -r myscript.plms
    RESTART: Moving stuff to stuff.res and restarting from it
    [18:03:34] PLAMS working folder: /home/user/stuff

(the above command needs to be executed in ``/home/user``.
Otherwise, you need to add ``-p /home/user`` to tell the master script where to look for ``stuff``).
The same example with the default folder name::

    $ plams myscript.plms
    [17:28:40] PLAMS working folder: /home/user/plams_workdir
    # [some successful work]
    [17:56:22] Execution interrupted by the following exception:
    # [exception details]

    [...debug the script...]

    $ plams -r myscript.plms
    RESTART: Moving plams_workdir to plams_workdir.res and restarting from it
    [18:03:34] PLAMS working folder: /home/user/plams_workdir

For more detailed explanation of the restart mechanism, please see |RPM|, |pickling| and |restarting|.


Multiple input scripts
~~~~~~~~~~~~~~~~~~~~~~~~~

The launch script can be called with more than one positional argument, like for example::

    plams script1.plms script2.plms script3.plms

All files supplied that way are concatenated into one script and then executed (that means things declared in script1 are visible in script2 and script3).
Using this feature for completely unrelated scripts is probably not a good idea, but it can be useful, for example, when first files contain just definitions of your own functions, derived classes, settings tweaks etc. that are then used in the last file::

    plams config/debug_run.plms settings/adf/adf_fde.plms actual_script.plms

That way you can build your own library of reusable code snippets for tasks that are most frequently occurring in your daily work, customize PLAMS according to your personal preferences and make your working environment truly modular.

.. note::

    The ``.plms`` file extension for PLAMS scripts is just a convention.
    Scripts can be any text files.


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
