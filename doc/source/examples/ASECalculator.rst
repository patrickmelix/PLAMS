.. _ASECalculatorExample:

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

* A single-point calculation run with the ASE calculator, storing output on disk

* An AMS geometry optimization job (``Task GeometryOptimization``) run with the ASE calculator, storing output on disk

* A geometry optimization through ASE (``ase.optimize.BFGS``), storing output for each frame on disk. AMS is launched once for each frame.

* A geometry optimization through ASE (``ase.optimize.BFGS``) using the "AMSworker" mode. In AMSworker mode, AMS is only launched at the start. No output is stored on disk. This minimizes the overhead.



**Example usage:** (:download:`Download ASECalculator.py <../../../examples/ASECalculator.py>`)

.. literalinclude:: ../../../examples/ASECalculator.py
	:language: python

