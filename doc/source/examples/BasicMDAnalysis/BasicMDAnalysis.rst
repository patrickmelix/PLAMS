.. _BasicMDPostanalysisExample:

Basic molecular dynamics analysis
=====================================

**Note**: This example requires AMS2023 or later.

This example illustrates how to calculate the basic

* velocity autocorrelation function (VACF)
* diffusion coefficient from the integral of the VACF
* power spectrum (Fourier transform of the VACF)
* viscosity from the Green-Kubo relation (integral of the off-diagonal pressure tensor autocorrelation function)
* density along the *z* axis

For details about the functions, see the |AMSResults| API.

.. important::

    The example only shows how to technically calculate the numbers.
    For real simulations, run longer MD simulations and carefully converge any
    calculated quantities.

    The viscosity requires especially long MD simulations.

.. note::

    More advanced analysis is possible by setting up an :class:`~scm.plams.interfaces.adfsuite.amsanalysis.AMSAnalysisJob` job.

    See also: `Molecular Dynamics with Python <../../../Tutorials/MolecularDynamicsAndMonteCarlo/MDintroPython/intro.html>`__

.. include:: BasicMDAnalysis.common_header.rst
.. include:: BasicMDPostanalysis.ipynb.rst
.. include:: BasicMDAnalysis.common_footer.rst