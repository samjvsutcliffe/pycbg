Welcome to pycbg's documentation!
=================================

This module helps configuring MPM simulations for `CB-Geo MPM <https://github.com/cb-geo/mpm>`_: a simple python script can generate the input files required by CB-Geo (see :ref:`Preprocessing`). The results of the simulation can also be imported in python using pycbg (see :ref:`Postprocessing`). This documentation should be used alongside `CB-Geo's documentation <https://mpm.cb-geo.com/#/>`_.

.. Preprocessing:

----

Preprocessing
=============

Preprocessing a MPM simulation for CB-Geo consist in creating several input files :

 - a mesh file, where the positions of all nodes and their interconnections are described. Pycbg saves it under the name `mesh.msh`. Can be created using the :py:class:`Mesh<pycbg.preprocessing.Mesh>` object.
 - a particles file, where the initial positions of all material points are specified. Pycbg saves it under the name `particles.txt`. Can be created using the :py:class:`Particles<pycbg.preprocessing.Particles>` object.
 - an entity sets file (if entity sets are defined), where all entity sets are defined using entities' ids. An entity can be a node, a particle or a cell. Pycbg saves it under the name `entity_sets.txt`. Can be created using the :py:class:`EntitySets<pycbg.preprocessing.EntitySets>` object.

The :py:class:`Simulation<pycbg.preprocessing.Simulation>` object creates :py:class:`Mesh<pycbg.preprocessing.Mesh>`, :py:class:`Particles<pycbg.preprocessing.Particles>`, :py:class:`Materials<pycbg.preprocessing.Materials>` and :py:class:`EntitySets<pycbg.preprocessing.EntitySets>` objects and should be enough to prepare a simulation.

Classes overview
----------------

.. currentmodule:: pycbg.preprocessing

.. autosummary::
   :nosignatures:
   :recursive:
   :template: custom-class-template.rst
   :toctree: stubs

   Mesh
   Particles
   EntitySets
   Materials
   Simulation

----

Postprocessing
==============

**TODO**: *explain postprocessing here*

.. currentmodule:: pycbg.postprocessing

.. autosummary::
   :nosignatures:
   :recursive:
   :template: custom-class-template.rst
   :toctree: stubs

   ResultsReader
   
----

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
