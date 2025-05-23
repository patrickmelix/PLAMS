AMS driver and engines
----------------------

.. currentmodule:: scm.plams.interfaces.adfsuite.ams

The AMS driver is a new program introduced in the 2018 release that unifies the way in which different computational engines of Amsterdam Modelling Suite are called.
You can find more information about the AMS driver in the `corresponding part of the documentation <../../AMS/General.html>`_.

.. _AMS_preparing_input:

Preparing input
~~~~~~~~~~~~~~~

.. tip::

    Starting with AMS2024, you can also use `PISA <../../pisa/index.html>`_
    (Python Input System for AMS) to specify the input to AMS.

    However, almost all PLAMS examples still use the input description described
    on this page.


.. note::

    Input files handling in the AMS driver is case insensitive.

The input file for the AMS driver consists of keys and values organized in blocks and subblocks:

.. code-block:: none

    Task GeometryOptimization

    GeometryOptimization
        Convergence
            Gradients 1.0e-4
        End
    End

    Properties
       NormalModes true
    End

    System
        Atoms
            C       0.00000000       0.00000000       0.00000000
            H       0.63294000      -0.63294000      -0.63294000
            H      -0.63294000       0.63294000      -0.63294000
            H       0.63294000       0.63294000       0.63294000
            H      -0.63294000      -0.63294000       0.63294000
        End
    End

    Engine DFTB
        Model DFTB3
        ResourcesDir DFTB.org/3ob-3-1
    EndEngine

Such a structure can be reflected in a natural way by a multi-level character of |Settings|.
The example input file presented above can be generated by::

    s = Settings()
    #AMS driver input
    s.input.ams.Task = 'GeometryOptimization'
    s.input.ams.GeometryOptimization.Convergence.Gradients = 1.0e-4
    s.input.ams.Properties.NormalModes = 'true'
    #DFTB engine input
    s.input.DFTB.Model = 'DFTB3'
    s.input.DFTB.ResourcesDir = 'DFTB.org/3ob-3-1'

    m = Molecule('methane.xyz')
    j = AMSJob(molecule=m, settings=s)
    j.run()

If an entry is a regular key-value pair it is printed in one line (like ``Task GeometryOptimization`` above).
If an entry is a nested |Settings| instance it is printed as a block and entries inside this instance correspond to the contents of that block.

One of the blocks is special: the engine block.
It defines the computational engine used to perform the task defined in all other blocks.
The contents of the engine block are not processed by the AMS driver, but rather passed to the corresponding engine instead.
Because of that every AMS input file can be seen as composed of two distinct parts: the engine input (everything **inside** the engine block) and the driver input (everything apart from the engine block).
That distinction is reflected by how a |Settings| instance for |AMSJob| is structured.
As we can see, the ``input`` branch of job settings is divided into two branches: the ``ams`` branch for the driver input and the ``DFTB`` branch for the engine input.

.. note::

    In general, PLAMS will use all the contents of the ``ams`` branch (spelling not case-sensitive) to construct the driver input and the contents of every other branch to construct a separate engine block with the same name as the branch (like ``DFTB`` in the example above).
    In the present moment only applications with a single engine block are implemented in the AMS driver, but that will most likely change in the near future.

The contents of each branch of ``myjob.settings.input`` are translated to a string using the same logic:

*   Entries within each block (including the top level) are listed in the alphabetical order.
*   Both keys and values are kept in their original case.
*   Strings used as values can contain spaces and all kinds of special characters, including new lines.
    They are printed in an unchanged form in the input file.
*   If you need to put a key without any value, you can use ``True`` or an empty string as a value::

        s.input.ams.block.key = True
        s.input.ams.otherkey = ''

        ### translates to:

        block
          key
        end

        otherkey

*   If a value of a key is ``False`` or ``None`` the key is omitted.
*   En empty |Settings| instance produces an empty block::

        s.input.ams.emptyblock = Settings()
        s.input.ams.otherblock #short syntax equivalent to the line above

        ### translates to:

        emptyblock
        end

        otherblock
        end

*   More instances of the same key within one block can be achieved by using a list of values instead of a single value::

        s.input.ams.constraints.atom = [1,5,4]
        s.input.ams.constraints.block = ['ligand', 'residue']

        ### translates to:

        constraints
          atom 1
          atom 5
          atom 4
          block ligand
          block residue
        end

*   Some blocks require (or allow) something to be put in the header line, next to the block name.
    Special key ``_h`` is helpful in these situations::

        s.input.ams.block._h = 'header=very important'
        s.input.ams.block.key1 = 'value1'
        s.input.ams.block.key2 = 'value2'

        ### translates to:

        someblock header=very important
          key1 value1
          key2 value2
        end

*   Another kind of special key can be used to override the default alphabetic ordering of entries within a block, or just to insert arbitrary strings into the block::

        s.input.ams.block._1 = 'entire line that has to be the first line of block'
        s.input.ams.block._2 = 'second line'
        s.input.ams.block._4 = 'I will not be printed'
        s.input.ams.block.key1 = 'value1'
        s.input.ams.block.key2 = 'value2'

        ### translates to:

        block
          entire line that has to be the first line of block
          second line
          key1 value1
          key2 value2
        end

*   If a value of a key needs to be a path to some KF file with results of a previous AMS calculation, an instance of |AMSJob| or |AMSResults| (or directly |KFFile|) can be used (see :meth:`AMSJob.get_input` for details)::

        oldjob = AMSJob(...)
        oldjob.run()
        newjob = AMSJob(...)
        newjob.settings.input.ams.loadsystem.file = oldjob
        newjob.settings.input.ams.loadengine = (oldjob, 'dftb')

        ### translates to:

        loadengine /home/user/plams_workdir/oldjob/dftb.rkf

        loadsystem
          file = /home/user/plams_workdir/oldjob/ams.rkf
        end

* Convert AMS text-style input to a Settings object (this requires that the SCM python package is installed)::

        text = '''
        Task GeometryOptimization
        Engine DFTB
           Model GFN1-xTB
        EndEngine
        '''

        sett = AMSJob.from_input(text).settings
        print(sett)

        # output:

        input:
              dftb:
                   model: 	GFN1-xTB
              ams:
                   task: 	GeometryOptimization



.. note::

    The algorithm translating |Settings| contents into an input file does not check the correctness of the given data - it simply takes keys and values from |Settings| and prints them in the text file.
    Due to that you are not going to be warned if you make a typo, use a wrong keyword or improper syntax.



Preparing runscript
~~~~~~~~~~~~~~~~~~~

Runscripts for the AMS driver are very simple (see :meth:`AMSJob.get_runscript`).
The only adjustable option (apart from usual ``pre``, ``post``, ``shebang`` and ``stdout_redirect`` which are common for all single jobs) is ``myjob.settings.runscript.nproc``, indicating the number of parallel processes to run AMS with (like with ``-n`` flag or ``NSCM`` environmental variable).



.. _AMSMoleculeHandling:

Molecule handling
~~~~~~~~~~~~~~~~~

There are several ways in which the description of the simulated system can be supplied to |AMSJob|.
The most convenient one is simply by passing a |Molecule| instance::

    mol = Molecule('/path/to/some/file.xyz')
    myjob = AMSJob(name='test', molecule=mol, settings=...)

or::

    mol = Molecule('/path/to/some/file.xyz')
    myjob = AMSJob(...)
    myjob.molecule = mol


.. note::
    
    Instead of passing a |Molecule| object to |AMSJob|, you have the option to use a `Chemical System <../../Scripting/LibBase/ChemicalSystem.html>`_ as well.


A |Molecule| instance stored as the ``molecule`` attribute is automatically processed during the input file preparation and printed in the proper format (see `AMS manual <../../AMS/System.html>`_ for details).
Various details of this process can be adjusted based on attributes of the supplied |Molecule|.
If ``mol.lattice`` is nonempty, the information about periodicity vectors is printed to the ``lattice`` subblock of the ``system`` block.
If the supplied lattice consists of 1 or 2 vectors that do not follow the convention required by AMS (1D -- vector aligned with X axis; 2D -- vectors aligned with XY plane) the whole system is rotated to meet these criteria.
If ``mol.properties.charge`` exists, it is used as the ``charge`` key in the ``system`` block.

Moreover, each |Atom| present in the supplied |Molecule| has its own ``properties`` attribute that can be used to adjust the details of the line generated for this atom in the ``atoms`` block:

*   The atomic symbol is generated based on the atomic number stored in the ``atnum`` attribute of the |Atom|.
    The atomic number of 0 corresponds to the "dummy atom" for which the symbol is empty.
*   If ``atom.properties.ghost`` exists and is ``True``, the atomic symbol is prefixed with ``Gh.``.
*   If ``atom.properties.name`` exists, the name is added after the atomic symbol, separated by a single dot.
*   Leading, trailing and double dots are removed from the atomic symbol.
*   If ``atom.properties.suffix`` exists, it is placed at the end of the line, after the numerical coordinates (it should ba a string)

Example::

    mol = Molecule('xyz/Ethanol.xyz')
    mol[1].properties.ghost = True
    mol[2].properties.name = 'D'
    mol[3].properties.ghost = True
    mol[3].properties.name = 'T'
    mol[4].properties.atnum = 0
    mol[4].properties.name = 'J.XYZ'
    mol[5].properties.atnum = 0
    mol[5].properties.name = 'J.ASD'
    mol[5].properties.ghost = True
    mol[6].properties.suffix = 'whatever text'
    myjob = AMSJob(molecule=mol)

The corresponding fragment of the input file produced by the above code:

.. code-block:: none

    system
      atoms
          1      Gh.C       0.01247       0.02254       1.08262
          2       C.D      -0.00894      -0.01624      -0.43421
          3    Gh.H.T      -0.49334       0.93505       1.44716
          4     J.XYZ       1.05522       0.04512       1.44808
          5  Gh.J.ASD      -0.64695      -1.12346       2.54219
          6         H       0.50112      -0.91640      -0.80440   whatever text
          7         H       0.49999       0.86726      -0.84481
          8         H      -1.04310      -0.02739      -0.80544
          9         O      -0.66442      -1.15471       1.56909
      end
    end

Another, more cumbersome way to provide the system information to |AMSJob| is to manually populate the ``system`` block in job settings::

    s = Settings()
    s.input.ams.system.atoms._1 = 'H     0.0    0.0     0.0'
    s.input.ams.system.atoms._2 = 'O     1.0    0.0     0.0'
    s.input.ams.system.charge = 1.0
    #other settings adjustments
    myjob = AMSJob(settings=s)

An alternative way of supplying molecular coordinates is to use the ``GeometryFile`` key in the ``system`` block::

    s = Settings()
    s.input.ams.system.geometryfile = '/path/to/some/file.xyz'
    #other settings adjustments
    myjob = AMSJob(settings=s)

Currently only the `extended XYZ format <../../AMS/Appendices.html#extendedxyz>`_ is supported.

Finally, one could use the ``LoadSystem`` top-level key and point to an existing ``.rkf`` file with results of some previous calculation::

    s = Settings()
    s.input.loadsystem = '/path/to/some/ams.rkf'
    #other settings adjustments
    myjob = AMSJob(settings=s)


Multiple molecules
++++++++++++++++++

The AMS driver allows multiple occurrences of the ``system`` block in the input file.
Different ``system`` blocks are distinguished by their names defined in the header of the block::

    system protein
        atoms
        ...
        end
    end
    system ligand
        atoms
        ...
        end
    end

The system without such a name is considered the main system.

Multiple systems can be used in |AMSJob| by setting the ``molecule`` attribute to a dictionary, instead of a single |Molecule|.
Such a dictionary should have strings as keys and |Molecule| instances as values.
The main system should have ``''`` (an empty string) as a key.

Other methods of providing the contents of the ``system`` block mentioned above can also be used to provide multiple ``system`` blocks.
``myjob.settings.input.ams.system`` can be a list containing multiple |Settings| instances, one for each system.
Each such instance can be have manually filled ``atoms`` block or use the ``geometryfile`` key.
Special header ``_h`` key can be used to set headers and hence names of different ``system`` blocks.
Multiple instances of the ``LoadSystem`` key (also provided as a list, also with ``_h`` headers) can also be used.

All the methods mentioned above (``molecule`` attribute, ``GeometryFile``, ``LoadSystem``, manual ``system`` block preparation) can be combined in any configuration.
In case of a conflict, the data stored in ``settings.input.ams.system`` takes precedence over ``molecule``.
It is, however, the user's responsibility to make sure that among all the systems provided there is exactly one main system (without a name).



AMSJob API
~~~~~~~~~~

.. autoclass:: AMSJob(name='plamsjob', molecule=None, settings=None, depend=None)
    :exclude-members: _result_type
    :no-private-members:



AMSResults API
~~~~~~~~~~~~~~

.. autoclass:: AMSResults
    :exclude-members: __init__
    :no-private-members:



Other functions
~~~~~~~~~~~~~~~~~~

.. autofunction:: hybrid_committee_engine_settings
    

