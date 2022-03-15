.. _PackmolInterface:

Packmol interface
~~~~~~~~~~~~~~~~~~

Packmol (`Packmol website <http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml>`__) is a program for creating liquid or gas mixtures. The PLAMS interface only supports

* uniform mixtures
* solid/liquid interfaces

There are four main functions:

* ``packmol_liquid`` (for liquids with only 1 component)
* ``packmol_mixture`` (for liquid mixtures)
* ``packmol_solid_liquid`` (for solid/liquid interface with 1 component in the liquid)
* ``packmol_solid_liquid_mixture`` (for solid/liquid interface with a liquid mixture)

See the :ref:`Packmol example <PackMolExample>` for all the ways these functions can be used.

The above four functions accept an ``executable`` argument, which should
contain the path to the packmol program. If it is not specified, the path to
the packmol program included with the Amsterdam Modeling Suite will be used.

.. currentmodule:: scm.plams.interfaces.molecule.packmol

.. autofunction:: packmol_liquid

.. autofunction:: packmol_mixture

.. autofunction:: packmol_solid_liquid

.. autofunction:: packmol_solid_liquid_mixture

.. autoclass:: PackMol

.. autoclass:: PackMolStructure


