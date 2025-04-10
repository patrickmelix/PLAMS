.. MapMoleculesAndConvertToDCD:

Convert RKF to DCD while mapping all atoms into molecules
=========================================================

A .rkf file can also be converted to a .dcd trajectory file.
DCD is a compact binary file format for trajectories, which can be used by a range of simulation, visualization and analysis programs.

The example below shows how to generate a DCD file from a completed molecular dynamics calculation in AMS.

**Example usage:** Download script :download:`rkf_to_wrapped_dcd.py <../../../examples/rkf_to_wrapped_dcd.py>` and example :download:`ams.rkf <../../../examples/BasicMDAnalysis/ams.rkf>` from a toy water MD simulation (or choose one of your own output files); then run as ``$AMSBIN/amspython rkf_to_wrapped_dcd.py``.

.. literalinclude:: ../../../examples/rkf_to_wrapped_dcd.py
	:language: python

