.. _ADFCOSMORSConformers:

.. |ConformersJob| replace:: :class:`ConformersJob<scm.conformers.ConformersJob>`


Generating multiple conformers for use with COSMO-RS 
----------------------------------------------------

The ``ADFCOSMORSConfJob`` is a customizable class that allows the user to design a conformer generation workflow for COSMO-RS.  The default instance of this class generates a set of conformers and then performs ADF and subsequent COSMO calculations (equivalent to the AMS Task:COSMO-RS Compound) to generate ``.coskf`` files for *all* unique conformers.  This class can be customized by adding both filters (to limit the number of conformers) and additional calculation steps (to improve the final geometry given to ADF and/or to increase the accuracy of energy calculations for filters).


Example usage
=============

First, we'll import the necessary classes:

.. code-block:: python

    from scm.plams.recipes.adfcosmorsconformers import ADFCOSMORSConfJob, ADFCOSMORSConfFilter
    from scm.plams import Molecule, from_smiles, init, finish, Settings
    from scm.conformers import ConformersJob

Now, we'll input the acetic acid molecule with the ``from_smiles`` function and specify the initial number of strucutres that generated with the default RDKitGenerator. 

.. code-block:: python

    mol = from_smiles("CC(=O)O")
    InitialConformers = 7

Let's also specify an additional step to add to the default workflow.  Here, we'll add a DFTB geometry optimization.

.. code-block:: python

    dftb_sett = Settings()
    dftb_sett.input.AMS.Task = "Optimize"
    dftb_sett.input.DFTB

The final thing we need to specify are filters.  Let's make two filters, the first to take a maximum of 5 conformers with a maximum energy range of 10 kcal/mol and the second with 2 conformers and 10 kcal/mol.

.. code-block:: python

    fil1 = ADFCOSMORSConfFilter(5,10)
    fil2 = ADFCOSMORSConfFilter(2,10)

Finally, we give this information to the ``ADFCOSMORSConfJob`` class.  We also specify the name of the coskf files as well as the directory in which we'll find them after the calculations complete.

.. code-block:: python

    a = ADFCOSMORSConfJob(
        mol,
        initial_conformers = InitialConformers
        first_filter = fil1,
        additional   = [(dftb_sett,fil2)],
        coskf_name   = "acetic_acid",
        coskf_dir    = "test_coskfs"
        )
    a.run()


In the ``test_coskfs`` directory, we find two ``.coskf`` files with the following geometries:

.. image:: /examples/CRSConformers/acetic_0.png
   :width: 50%

.. image:: /examples/CRSConformers/acetic_1.png
   :width: 45%

The entire script can be seen/copied when expanded below.

.. raw:: html

   <details>
   <summary style="color:#008000;cursor:pointer">[show/hide code]</summary>

.. literalinclude:: /examples/CRSConformers/test_confs.py
    :language: python
    
.. raw:: html

     </details>


Important Information
======================

Specifying Filters
##################

Filters are given as instances of the ``ADFCOSMORSConfFilter`` class.  This class takes two arguments concerning the number of conformers and maximum energy range.  These are documented in more detail below.  There are 3 places to specify filters:

+ in the ``first_filter`` keyword in ``ADFCOSMORSConfJob``.  This filter will be applied to the initial set of sampled conformers.
+ in the ``final_filter`` keyword in ``ADFCOSMORSConfJob``.  This filter will be applied to the results of the ADF gas phase calculation.
+ in the ``additional`` keyword as the second element of a (|Settings|, ``ADFCOSMORSConfFilter`` ) tuple.  The filter will be applied on the results of the calculation done with the |Settings| from the *same* tuple.  To be perfectly clear, this means a (|Settings|, ``ADFCOSMORSConfFilter`` ) tuple adds a step to the workflow where a calculation is done with the |Settings| object and then those results are immediately filtered with the corresponding ``ADFCOSMORSConfFilter``.  A list of these tuples can be given to perform a series of calculations/filtering steps.

Specifying Settings for additional jobs
#######################################

|Settings| objects provided to the ``additional`` keyword should be valid.  ``ADFCOSMORSConfJob`` will likely fail if there are errors in the |Settings| objects.  In the case of unexpected failures, take a look in the calculation subdirectories called ``additional_i`` where ``i`` is an index representing the index of the additional calculation.

.. note::
    The tasks for **Settings objects should correspond to Tasks for AMSConformers**.  The most important two types are **Optimize** to perform a geometry optimization and **Score** to perform a single point calculation.  `The tasks section <../../AMS/Utilities/Conformers.html#tasks>`__ in the Conformers doc contains all additional information.

Specifying conformer sampling strategies
########################################

The default conformer sampling strategy is based on RDKit.  This is likely sufficient for more use cases, but the user can specify another sampling strategy by providing a |ConformersJob| instance to the ``conf_gen`` keyword.  See also `The generator section <../../AMS/Utilities/Conformers.html#generators>`__ of the conformer documentation for more information about specifying alternative generators.  

Brief API Documentation
=======================

.. automodule:: scm.plams.recipes.adfcosmorsconformers
    :no-special-members:
    :exclude-members: _result_type, __init__, new_children, postrun
