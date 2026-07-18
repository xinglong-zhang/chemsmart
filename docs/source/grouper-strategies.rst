#####################
 Grouping Strategies
#####################

This page provides detailed documentation for each molecular structure grouping strategy available in chemsmart.

***********************
 RMSD-based Strategies
***********************

rmsd (Basic RMSD)
=================

Simple Kabsch RMSD alignment without atom reordering.

.. code:: bash

   chemsmart run grouper -f conformers.xyz rmsd

**Default threshold:** 0.5 Å

This is the fastest RMSD method but assumes atoms are already in the same order. Best for structures from the same
source with consistent atom ordering.

hrmsd (Hungarian RMSD)
======================

RMSD with optimal atom mapping using the Hungarian algorithm (*J. Chem. Inf. Model. 2014, 54, 518*).

.. code:: bash

   chemsmart run grouper -f conformers.xyz hrmsd

**Default threshold:** 0.5 Å

Finds the optimal one-to-one mapping between atoms of the same element type that minimizes RMSD. Useful when atom
ordering may differ.

spyrmsd
=======

RMSD calculation using the spyrmsd library with symmetry handling (*J. Cheminformatics 2020, 12, 49*).

.. code:: bash

   chemsmart run grouper -f conformers.xyz spyrmsd

**Default threshold:** 0.5 Å

Handles molecular symmetry and equivalent atom permutations using graph isomorphism.

irmsd (Invariant RMSD)
======================

Invariant RMSD that considers molecular symmetry and equivalent atom permutations (*J. Chem. Inf. Model. 2025, 65,
4501*).

.. code:: bash

   chemsmart run grouper -f conformers.xyz irmsd
   chemsmart run grouper -f conformers.xyz -T 0.3 irmsd --inversion on

**Default threshold:** 0.125 Å

Strategy-specific Options
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--inversion``
      -  choice
      -  Coordinate inversion checking: auto (default), on, off

**Note:** Requires the external ``irmsd`` package installed in a separate conda environment with numpy>=2.0. See
installation instructions in the README.

pymolrmsd
=========

PyMOL-based RMSD alignment using the ``align`` command (*Schrödinger, L. The PyMOL Molecular Graphics Development
Component, Version 1.8; 2015*).

.. code:: bash

   chemsmart run grouper -f conformers.xyz pymolrmsd

**Default threshold:** 0.5 Å

Uses PyMOL's powerful alignment algorithm. Note that PyMOL only supports single-threaded operation (``-np 1``).

******************************
 Fingerprint-based Strategies
******************************

tfd (Torsion Fingerprint Deviation)
===================================

Compare conformers based on their torsion angle fingerprints (*J. Chem. Inf. Model. 2012, 52, 1499*).

.. code:: bash

   chemsmart run grouper -f conformers.xyz tfd
   chemsmart run grouper -f conformers.xyz -T 0.2 tfd --no-use-weights
   chemsmart run grouper -f conformers.xyz tfd --max-dev spec

**Default threshold:** 0.1

Particularly effective for flexible molecules where conformational differences are primarily due to bond rotations.

Strategy-specific Options
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--use-weights/--no-use-weights``
      -  bool
      -  Use torsion weights in TFD calculation (default: True)

   -  -  ``--max-dev``
      -  choice
      -  Normalization method: equal (default) or spec

tanimoto
========

Group structures using Tanimoto similarity with molecular fingerprints.

.. code:: bash

   chemsmart run grouper -f conformers.xyz tanimoto
   chemsmart run grouper -f conformers.xyz -T 0.85 tanimoto --fingerprint-type morgan

**Default threshold:** 0.9

Strategy-specific Options
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-ft, --fingerprint-type``
      -  choice
      -  Fingerprint type (default: rdkit)

Available fingerprint types:

-  ``rdkit`` / ``rdk``: RDKit topological fingerprints
-  ``morgan``: Morgan circular fingerprints (ECFP-like)
-  ``maccs``: MACCS keys
-  ``atompair``: Atom pair fingerprints
-  ``torsion``: Topological torsion fingerprints
-  ``usr``: Ultrafast Shape Recognition (3D shape-based)
-  ``usrcat``: USR with Credo Atom Types (3D shape-based)

*************************
 Energy-based Strategies
*************************

energy
======

Group molecules based on their energy differences.

.. code:: bash

   chemsmart run grouper -f conformers.xyz energy
   chemsmart run grouper -f conformers.xyz -T 0.5 energy

**Default threshold:** 1.0 kcal/mol

**Note:** All molecules must have energy information. For xyz files, energy is read from the comment line. For log
files, Gibbs free energy is calculated as SCF energy + thermal correction.

******************
 Other Strategies
******************

formula
=======

Group molecules by their molecular formula.

.. code:: bash

   chemsmart run grouper -f molecules.xyz formula

This is a quick way to separate different compounds in a mixed set.

connectivity
============

Group molecules by their molecular connectivity/topology.

.. code:: bash

   chemsmart run grouper -f molecules.xyz connectivity

Groups molecules that have the same bond connectivity graph, regardless of 3D geometry.

isomorphism
===========

Group molecules based on graph isomorphism using RDKit.

.. code:: bash

   chemsmart run grouper -f molecules.xyz isomorphism

Similar to connectivity but uses RDKit's isomorphism checking.

******************
 Complete Linkage
******************

All threshold-based grouping strategies use **complete linkage clustering**, which ensures that all members within a
group are within the threshold distance of each other. This prevents the "chaining effect" where dissimilar structures
could end up in the same group through intermediate structures.
