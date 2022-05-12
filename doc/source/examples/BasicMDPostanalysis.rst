.. _BasicMDPostanalysisExample:

Basic molecular dynamics analysis
=====================================

This example illustrates how to calculate the basic

* velocity autocorrelation function (VACF)
* diffusion coefficient from the integral of the VACF
* power spectrum (Fourier transform of the VACF)
* viscosity from the Green-Kubo relation (integral of the off-diagonal pressure tensor autocorrelation function)

For details about the functions, see the |AMSResults| API.

.. note::

    More advanced analysis possibilities of these quantities is possible by setting up an :class:`~scm.plams.interfaces.adfsuite.amsanalysis.AMSAnalysisJob` job.


**Example usage:** (:download:`Download BasicMDPostanalysis.py <../../../examples/BasicMDPostanalysis.py>`)

.. literalinclude:: ../../../examples/BasicMDPostanalysis.py
	:language: python

