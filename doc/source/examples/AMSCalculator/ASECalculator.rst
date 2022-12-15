.. _ASECalculatorExample:
.. _AMSCalculatorExample:

ASE geometry optimizer with AMS forces
===============================================================

.. note::

    Follow this example only if you need to use ASE. If you just want to run a
    normal AMS single-point or geometry optimization, do not go through ASE but
    instead see the :ref:`GeoOptWaterExample` example.

Example illustrating how to use the :ref:`AMSCalculator`. The ``BFGS`` geometry optimizer from ASE is used together with AMS-calculated forces (negative gradients).

.. note::

    This example launches AMS in "AMSworker" mode. This means that AMS is only started at the beginning of the calculation.

To follow along, either

* Download :download:`AMSCalculatorWorkerMode.py` (run as ``$AMSBIN/amspython AMSCalculatorWorkerMode.py``).
* Download :download:`AMSCalculatorWorkerMode.ipynb` (see also: how to install `Jupyterlab <../../../Scripting/Python_Stack/Python_Stack.html#install-and-run-jupyter-lab-jupyter-notebooks>`__ in AMS)


.. include:: AMSCalculatorWorkerMode.rst.include


Complete Python code
----------------------------

.. literalinclude:: AMSCalculatorWorkerMode.py
	:language: python





