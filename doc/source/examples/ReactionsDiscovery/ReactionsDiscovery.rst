.. _rd_example:

Reactions Discovery
=======================

.. seealso::
   
   * The `Reactions Discovery documentation <../../../Workflows/ReactionsDiscovery/ReactionsDiscovery.html>`__ in AMS.

Example illustrating how to use Reactions Discovery with AMS.

**Question to be answered**:  NH₂-CH₂-CH₂-OH + CO₂ + H₂O → ???

**Answer from this example**: Side products include NH₃, NH₂-CH₂-CH=O, OH-NH-CH₂-CH₂-OH, ...

To follow along, either

* Download :download:`reactions_discovery.py` (run as ``$AMSBIN/amspython reactions_discovery.py``).
* Download :download:`reactions_discovery.ipynb` (see also: how to install `Jupyterlab <../../../Scripting/Python_Stack/Python_Stack.html#install-and-run-jupyter-lab-jupyter-notebooks>`__ in AMS)

.. note::

    Reactions Discovery depends on randomly filling a space with molecules. You are likely to get somewhat different results if you run this example!

.. include:: reactions_discovery.rst.include

Complete Python code
----------------------------

.. literalinclude:: reactions_discovery.py
    :language: python
