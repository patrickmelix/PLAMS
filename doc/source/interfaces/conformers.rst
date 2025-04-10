.. _conformers_interface:

Conformers
==========

AMS's `Conformers <../../AMS/Utilities/Conformers.html>`_ is a flexible tool for conformers generation.

This page documents the PLAMS interface to **Conformers**. For a description of the capabilities and options of the Conformers tool, see the `documentation <../../AMS/Utilities/Conformers.html>`_ in the AMS user manual.

.. seealso::

   The :ref:`conformers generation example <ConformersGenerationExample>` in the :ref:`examples` section of the PLAMS manual.

ConformersJob
-------------

.. important::

    Import these classes from ``scm.conformers``, not ``scm.plams`` !

    .. code-block::

        from scm.conformers import ConformersJob, ConformersResults

The ``ConformersJob`` class, which derives from |SingleJob| class, can be used to set up and run a Conformers calculation.

The input options for the Conformers tool (described `here <../../AMS/Utilities/Conformers.html>`_) can be specified in the ``input.ams`` branch of a setting object. See the :ref:`ConformersGenerationExample` example.

.. autoclass :: scm.conformers.ConformersJob

.. autoclass :: scm.conformers.ConformersResults

The ``plot_conformers()`` function (``from scm.conformers.plams.plot import plot_conformers``) lets you plot some example conformers in a Jupyter notebook:

.. autofunction:: scm.conformers.plams.plot.plot_conformers

