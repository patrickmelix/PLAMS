.. _IRSpectrumFromMD:

IR spectrum from MD
=======================

Example illustrating how to calculate an IR spectrum from the dipole derivative autocorrelation function with AMS and PLAMS.

The dipole moment is stored at every time step (more often than the structure) by setting the BinLog DipoleMoment option.

* Download :download:`ir_spectrum_from_md.py` (run as ``$AMSBIN/amspython ir_spectrum_from_md.py``).
* Download :download:`ir_spectrum_from_md.ipynb` (see also: how to install `Jupyterlab <../../../Scripting/Python_Stack/Python_Stack.html#install-and-run-jupyter-lab-jupyter-notebooks>`__ in AMS)

.. include:: ir_spectrum_from_md.rst.include

Complete Python code
----------------------------

.. literalinclude:: ir_spectrum_from_md.py
    :language: python
