.. _IonicConductivityFromMD:

Ionic conductivity from ams.rkf trajectory
==========================================

First, an ams.rkf trajectory with ions in the system needs to be calculated, as is done, for example, in the
`Ionic Conductivity Tutorial <../../Tutorials/MolecularDynamicsAndMonteCarlo/IonicConductivity.html>`__.
Next, the ionic conductivity for the ions in the solution can be computed. The below script first identifies the ions in the molecular system, and guesses their charges (in the case of metals it uses a simple charge equilibration scheme). It then runs the AMS trajectory analysis tool to compute the ionic conductivity for these ions.

**Example usage:** (:download:`Download get_ionic_conductivity.py <../../../examples/IonicConductivity/get_ionic_conductivity.py>`)

.. code-block:: none

   $AMSBIN/amspython get_ionic_conductivity.py /path/to/ams.rkf

To create an example ams.rkf file, a PLAMS script can be found here (:download:`Download NaClwater.py <../../../examples/IonicConductivity/NaClwater.py>`). This script is designed for efficiency. For a more reliable setup, follow the `Ionic Conductivity Tutorial <../../Tutorials/MolecularDynamicsAndMonteCarlo/IonicConductivity.html>`__.

.. literalinclude:: ../../../examples/IonicConductivity/get_ionic_conductivity.py
	:language: python

