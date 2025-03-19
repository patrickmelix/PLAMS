.. _ADFCOSMORSConformers:

.. |ConformersJob| replace:: :class:`ConformersJob<scm.conformers.ConformersJob>`

Generating multiple conformers for use with COSMO-RS
====================================================

The |ADFCOSMORSConfJob| is a customizable class that allows the user to design a conformer generation workflow for COSMO-RS.  The default instance of this class generates a set of conformers and then performs ADF and subsequent COSMO calculations (equivalent to the AMS Task:COSMO-RS Compound) to generate ``.coskf`` files for *all* unique conformers.  This class can be customized by adding both filters (to limit the number of conformers) and additional calculation steps (to improve the final geometry given to ADF and/or to increase the accuracy of energy calculations for filters). The generated ``.coskf`` can subsequently be utilized by the multi-species COSMO-RS, which accounts for compounds existing in multiple conformations. The application of multi-species COSMO-RS can be explored through  tutorials available in both the `GUI <../../Tutorials/COSMO-RS/COSMO-RS_multispecies.html#cosmo-rs-with-multi-species-components>`__ and `python scripting with pyCRS <../../COSMO-RS/Examples/conformer_CRSManager.html#adding-conformers-to-the-database>`__.

This example will produce two ``.coskf`` files with the following geometries:

.. image:: ../CRSConformers/acetic_0.png
   :width: 50%

.. image:: ../CRSConformers/acetic_1.png
   :width: 45%

.. include:: COSMORSConformers.common_header.rst
.. include:: cosmors_conformers.ipynb.rst
.. include:: COSMORSConformers.common_footer.rst

Important Information
---------------------

Specifying Filters
~~~~~~~~~~~~~~~~~~

Filters are given as instances of the ``ADFCOSMORSConfFilter`` class.  This class takes two arguments concerning the number of conformers and maximum energy range.  These are documented in more detail below.  There are 3 places to specify filters:

+ in the ``first_filter`` keyword in |ADFCOSMORSConfJob|.  This filter will be applied to the initial set of sampled conformers.
+ in the ``final_filter`` keyword in |ADFCOSMORSConfJob|.  This filter will be applied to the results of the ADF gas phase calculation.
+ in the ``additional`` keyword as the second element of a (|Settings|, ``ADFCOSMORSConfFilter`` ) tuple.  The filter will be applied on the results of the calculation done with the |Settings| from the *same* tuple.  To be perfectly clear, this means a (|Settings|, ``ADFCOSMORSConfFilter`` ) tuple adds a step to the workflow where a calculation is done with the |Settings| object and then those results are immediately filtered with the corresponding ``ADFCOSMORSConfFilter``.  A list of these tuples can be given to perform a series of calculations/filtering steps.

Specifying Settings for additional jobs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|Settings| objects provided to the ``additional`` keyword should be valid.  |ADFCOSMORSConfJob| will likely fail if there are errors in the |Settings| objects.  In the case of unexpected failures, take a look in the calculation subdirectories called ``additional_i`` where ``i`` is an index representing the index of the additional calculation.

.. note::
    The tasks for **Settings objects should correspond to Tasks for AMSConformers**.  The most important two types are **Optimize** to perform a geometry optimization and **Score** to perform a single point calculation.  `The tasks section <../../AMS/Utilities/Conformers.html#tasks>`__ in the Conformers doc contains all additional information.

Specifying conformer sampling strategies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The default conformer sampling strategy is based on RDKit.  This is likely sufficient for more use cases, but the user can specify another sampling strategy by providing a |ConformersJob| instance to the ``conf_gen`` keyword.  See also `The generator section <../../AMS/Utilities/Conformers.html#generators>`__ of the conformer documentation for more information about specifying alternative generators.

Brief API Documentation
-----------------------

.. automodule:: scm.plams.recipes.adfcosmorsconformers
    :no-special-members:
    :exclude-members: _result_type, __init__, new_children, postrun
