.. _CustomASECalculatorExample:

Couple any ASE calculator to AMS Driver
==============================================

This example illustrates how set up an AMS calculation using any external ASE calculator as the engine.

In this example, the EMT calculator included with ASE is used.

.. note::

    In this example, the AMS *engine* is replaced by an ASE calculator.

    In the :ref:`AMSCalculatorExample`, :ref:`iPIExample`, and
    :ref:`SellaExample` examples, the AMS *driver* is replaced by external
    codes that use the AMS *engine* energies and forces.

**Example usage:** (:download:`Download CustomASECalculator.py <../../../examples/CustomASECalculator.py>`)

.. literalinclude:: ../../../examples/CustomASECalculator.py
	:language: python

