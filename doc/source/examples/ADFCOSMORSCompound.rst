.. _ADFCOSMORSCompound:

ADF: Task COSMO-RS Compound
===========================

The ``ADFCOSMORSCompound`` class generates results identical to the "Task COSMO-RS Compound" in the AMS ADF graphical user interface.  This python interface allows users to easily generate the .coskf files for one or many structures.  A possible usage is given in :ref:`ams_crs_workflow`.


Example: generating .coskf files for a set of compounds using parallel calculation
-----------------------------------------------------------------------------------

First, we'll import the necessary classes and enable the parallel calculation through ``JobRunner``. Here, we'll assign one core to each job, and we can have up to eight jobs running all at once

.. code-block:: python

    from scm.plams import init, finish, from_smiles, JobRunner, config
    from scm.plams.recipes.adfcosmorscompound import ADFCOSMORSCompoundJob

    init()

    config.default_jobrunner = JobRunner(parallel=True, maxjobs=8)        # Set the default jobrunner to be parallel 
    config.default_jobmanager.settings.hashing = None                     # Disable rerun prevention
    config.job.runscript.nproc = 1                                        # Number of cores for each job
    config.log.stdout = 1    

Now, we will specify the smiles and name of a set of compounds and generate the intial geometry of each compound using ``from_smiles`` function. With the setting, ``nconfs=100`` and ``forcefield='uff'``, we will generate 100 conformers and find the one with the loweset energy using 'uff' forcefield.

.. code-block:: python

    rd_smiles = ['O'  ,'CO']
    rd_names  = ['H2O','CO']
    rd_mol={}
    for name, smiles in zip(rd_names, rd_smiles):
        rd_mol[name] = molecule = from_smiles(smiles, nconfs=100, forcefield='uff')[0] #lowest energy in 100 conformers

Lastly, we give this information to the ``ADFCOSMORSCompound`` class, including the name of the coskf files as well as the directory in which we'll find them after the calculations complete.  Using the setting, ``preoptimization='GFN1-xTB'`` and ``singlepoint=False``, it will utilize the DFTB for a quick pre-optimization. Subsequently, it will execute a gas phase optimization using ADF, followed by the solvation calculation.

.. code-block:: python

    results = []
    for name, molecule in mol.items():
        job = ADFCOSMORSCompoundJob(
            molecule        = molecule,           # The initial structure
            coskf_name      = name,               # a name to be used for coskf file
            coskf_dir       = 'test_coskfs',      # a directory to put the .coskf files generated
            preoptimization = 'GFN1-xTB',         # perform preoptimize or not 
            singlepoint     = False,              # run a singlepoint in gasphase and solvation calculation without geometry optimization. Cannot be combined with ``preoptimization``
            name            = name              ) # an optional name for the calculation directory
        results.append(job.run())

    finish()

In the ``test_coskfs`` directory, we will find the `H2O.coskf` and `CO.coskf` files.


The entire script can be seen/copied when expanded below.

.. raw:: html

   <details>
   <summary style="color:#008000;cursor:pointer">[show/hide code]</summary>

.. literalinclude:: /examples/CRScosmorscompound/parallel_COSMORScompound.py
    :language: python
    
.. raw:: html

     </details>



Source code for ``ADFCOSMORSCompound``
--------------------------------------

.. raw:: html

   <details>
   <summary style="color:#008000;cursor:pointer">[show/hide code]</summary>

.. literalinclude:: ../../../recipes/adfcosmorscompound.py

.. raw:: html

     </details>

Brief API Documentation
-----------------------

.. automodule:: scm.plams.recipes.adfcosmorscompound
    :no-special-members:
    :exclude-members: _result_type, __init__, new_children, postrun, _get_radii, adf_settings