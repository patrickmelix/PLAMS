.. _ASECalculatorExample:
.. _AMSCalculatorExample:

ASE Calculator for AMS
======================

.. note::

    Follow this example only if you need to use ASE. If you just want to run a
    normal AMS single-point or geometry optimization, do not go through ASE but
    instead see the :ref:`GeoOptWaterExample` example.

This example illustrates how set up an AMS engine for use with ASE. PLAMS is used in four ways:

* ``init()`` is called to initialize the PLAMS workdir (``plams_workdir/``). Outputs from individual AMS jobs are stored there.

* The settings for the AMS job (e.g. which engine to use and its settings) are defined with a PLAMS ``Settings`` object.

* The ``from_smiles()`` function is used to generate an initial structure

* The ``toASE()`` function is used to convert from a PLAMS ``Molecule`` to the ASE ``Atoms``. You can also do the reverse: ``fromASE()`` converts from ASE ``Atoms`` to PLAMS ``Molecule``.

The example shows

* **singlepoint**: A single-point calculation run with the ASE calculator, storing output on disk

* **ams_geoopt**: An AMS geometry optimization job (``Task GeometryOptimization``) run with the ASE calculator, storing output on disk.

* **ase_geoopt**: A geometry optimization through ASE (``ase.optimize.BFGS``), storing output for each frame on disk. AMS is launched once for each frame. This is quite slow and discouraged.

* **ase_geoopt_workermode**: A geometry optimization through ASE (``ase.optimize.BFGS``) using the "AMSworker" mode. In AMSworker mode, AMS is only launched at the start. No output is stored on disk. This minimizes the overhead. See also the :ref:`iPIExample` and :ref:`SellaExample` examples.

.. note::

    In this example, the AMS *driver* is replaced by ASE tools, that use the
    AMS *engines* (ADF, BAND, DFTB, ForceField, ...) to calculate energies and
    forces.

    In the :ref:`CustomASECalculatorExample` example, the AMS *engines* are
    instead replaced with external ASE calculators, that can be coupled to AMS
    *driver* tasks (GeometryOptimization, TransitionStateSearch,
    MolecularDynamics, ...).


**Example usage:** (:download:`Download ASECalculator.py <../../../examples/ASECalculator.py>`)

.. literalinclude:: ../../../examples/ASECalculator.py
	:language: python

