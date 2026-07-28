###############################
 Conformational Search (CREST)
###############################

This page covers conformational sampling workflows using `CREST <https://crest-lab.github.io/crest-docs/>`_
(Conformer-Rotamer Ensemble Sampling Tool).

For a complete list of CLI options, see :doc:`crest-cli-options`.

****************************
 Free Conformational Search
****************************

Sample the conformational landscape of a molecule using CREST.

.. code:: bash

   chemsmart sub [OPTIONS] crest [CREST_OPTIONS] conformers

Basic Usage
===========

Standard conformational search:

.. code:: bash

   chemsmart sub -s server crest -p project -f struct.xyz -c 0 -m 1 conformers

***********************************
 Constrained Conformational Search
***********************************

Sample conformations while maintaining fixed distance constraints. This is useful for transition state conformer
searches or other applications where certain geometric parameters must be held constant.

.. code:: bash

   chemsmart sub [OPTIONS] crest [CREST_OPTIONS] conformers [SUBCMD_OPTIONS]

Constrained Search Options
==========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--constraints``
      -  string
      -  Coordinate indices to constrain (1-based indexing).

   -  -  ``--force-constant``
      -  float
      -  Harmonic restraint force constant in Hartree/Bohr².

Basic Usage
===========

Constrained search with distance constraints:

.. code:: bash

   chemsmart sub -s server crest -p project -f molecule.log conformers \
      --constraints [[1,2],[2,3],[3,5]]

With explicit force constant:

.. code:: bash

   chemsmart sub -s server crest -p project -f molecule.log conformers \
     --constraints [[1,2],[2,3]] \
     --force-constant 0.25

How Constraints Work
====================

When you specify ``--constraints [[1,2],[2,3],[3,5]]``, CHEMSMART:

#. Reads the input geometry and measures the distances for the specified atom pairs
#. Generates a ``constraints.inp`` file with those distances and force constant if available
#. Passes this to CREST via ``--cinp constraints.inp``

Example generated ``constraints.inp``:

.. code:: text

   $constrain
   force constant=0.25
     distance: 1, 2, 1.4923
     distance: 2, 3, 1.9267
     distance: 3, 5, 1.6804
   $end

***************
 Output Layout
***************

After a successful CREST run, the output directory typically includes:

.. list-table::
   :header-rows: 1

   -  -  File
      -  Description
   -  -  ``crest_conformers.xyz``
      -  Conformer ensemble, sorted by energy
   -  -  ``crest_best.xyz``
      -  Lowest-energy conformer, typically used for subsequent calculations.
   -  -  ``crest_rotamers.xyz``
      -  Rotamer ensemble
   -  -  ``crest.energies``
      -  Relative energy values for all conformers (in kcal/mol)
   -  -  ``{label}.out``
      -  Main CREST log output with sampling details
   -  -  ``constraints.inp``
      -  Constraint file (only for constrained searches)

The generated conformers can be visualized or reused as input for subsequent calculations (e.g. ``-f
path/to/crest_conformers.xyz``).

************
 Next Steps
************

-  :doc:`crest-cli-options` — Full CLI option reference
-  :doc:`grouper-crest-or-traj-workflow` — Process conformer ensembles
-  :doc:`gaussian-conformational-sampling` — Run Gaussian calculations on conformers
-  :doc:`pymol-visualization` — visualize CREST conformer or rotamer ensembles
