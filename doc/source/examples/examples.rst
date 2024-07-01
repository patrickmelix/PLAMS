.. _examples:

Examples
========

In this chapter we present example PLAMS scripts covering various applications, from very simple tasks (like running the same calculation for multiple molecules) to more advanced dynamic workflows.

The example scripts use computational engines from the Amsterdam Modeling Suite, and you will need a license to run them. Contact license@scm.com for further questions.

In order to run the examples, the ``AMSBIN`` environment variable should be properly set. You can test this by typing ``$AMSBIN/plams -h`` in a terminal: this should print PLAMS' help message. If this is not the case (e.g. you get 'No such file or directory'), you need to set up the environmental variable ``$AMSBIN`` (see the `Linux Quickstart guide <../../Installation/Linux_Quickstart_Guide.html>`__ for details).

Getting Started
----------------

.. toctree::
 
   WaterOptimization
   Settings/AMSSettingsSystem
   He2DissociationCurve
   ManyJobsInParallel

Molecule analysis
---------------------

.. toctree::

   MoleculesFromRKFTrajectory
   MoleculesTable
   MoleculeSubstitution/MoleculeSubstitutionExample
   ConvertToAMSRKFTrajectory
   PlotCorrelation/PlotCorrelation
   MapMoleculesAndConvertToDCD
   HydrogenBondsFromMD

MD trajectory analysis
----------------------

.. toctree::

  IonicConductivityFromMD.rst
  IRSpectrumFromMD/IRSpectrumFromMD.rst
  IRSpectrumFromMDH2ODimer.rst

Benchmarks
-----------------

.. toctree::

   BasisSetBenchmark
   ReactionEnergyBenchmark

Workflows
------------------

.. toctree::
   :maxdepth: 1

   RedoxPotential
   ExcitationsWorkflow
   AMSTSWorkflow/AMSTSWorkflow
   ChargeTransferIntegralsADF
   gammascan
   ConformersGeneration/ConformersGeneration
   ConformersMultipleMolecules/ConformersMultipleMolecules
   ReactionsDiscovery/ReactionsDiscovery

COSMO-RS and property prediction
-----------------------------------

For more examples, see the `COSMO-RS documentation <../../COSMO-RS/Python_Examples.html>`__.

.. toctree::

   PropertyPrediction/PropertyPrediction
   ams_crs

Packmol and AMS-ASE interfaces
-------------------------------

.. toctree::

   PackMolExample/PackMolExample
   CustomASECalculator
   AMSCalculator/ASECalculator
   AMSCalculator/ChargedAMSCalculator
   i-PI-AMS
   SellaTransitionStateSearch

ParAMS, Simple Active Learning, and pyZacros
---------------------------------------------------

See the respective documentation pages:

* `ParAMS <../../params/index.html>`__ 
* `Simple Active Learning <../../Workflows/SimpleActiveLearning/SimpleActiveLearning.html>`__
* `pyZacros <../../pyzacros/index.html>`__ 

Other AMS calculations
------------------------

.. toctree::

   BAND_NiO_HubbardU
   BandStructure/BandStructure
   AMSPlumedMD/AMSPlumedMD
   QE_AMS_AFM_HubbardU
   BasicMDPostanalysis
   UseLowestEnergy
   M3GNet

Pymatgen
----------------

.. toctree::
    
    XRD/XRD

.. _recipes:

Pre-made recipes
----------------

The examples presented in here are simple job types built using basic PLAMS elements.
They are shipped with PLAMS in the ``recipes`` subpackage and can be directly used in your scripts.
In other words, the code presented there is already included in PLAMS and (unlike examples from two other sections) does not need to be copied to your script.
The source code of ``recipes`` modules is presented here to demonstrate how easy it is to build on top of existing PLAMS elements and create your own fully customized job types.

.. toctree::

    ADFCOSMORSCompound
    ADFCOSMORSConformers
    MDJobs
    adffragment
    ReorganizationEnergy
    adfnbo
    numgrad
    numhess
    pyAHFCDOS
    fcf_dos
    ReuseForceFieldParams

