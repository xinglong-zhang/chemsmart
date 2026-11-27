#####################
 Sterimol Parameters
#####################

Sterimol parameters describe the shape of a substituent with three lengths measured around the bond that attaches it to
the rest of the molecule. CHEMSMART calculates them from a geometry alone, with no quantum chemistry program. Run the
calculation on the current machine with ``chemsmart run``, or submit it to a compute node with ``chemsmart sub``.

*****************************
 What the Parameters Measure
*****************************

Picture the substituent as a set of van der Waals spheres, and the attachment bond as an axis running from the parent
atom into the substituent.

.. list-table::
   :header-rows: 1
   :widths: 10 50 40

   -  -  Parameter
      -  What it measures
      -  What it tells you

   -  -  ``L``
      -  Furthest extent of the substituent along the axis
      -  How far the group reaches away from the parent

   -  -  ``B1``
      -  Smallest width perpendicular to the axis: over all directions around the axis, the shortest distance from the
         axis to the outer surface
      -  How thin the group is on its narrowest side

   -  -  ``B5``
      -  Largest distance from the axis to the outer surface
      -  How far the group sticks out sideways at its widest

All three are in Å.

***************
 Quick Example
***************

.. code:: bash

   chemsmart run sterimol --file molecule.xyz --axis 2 4 --auto-fragment

.. code:: text

   Sterimol parameters: molecule.xyz (C5H10ClF)
   Axis:         C2 -> C4
   Substituent:  4-7,15-17 (auto-fragment)
   L  = 5.758 Å
   B1 = 1.763 Å
   B5 = 3.856 Å
   Radii:        RDKit periodic table (van der Waals radii)
   L convention: Measured along the axis from the first axis atom, without Verloop's +0.40 Angstrom correction

This measures the group on the atom 4 side of the bond between atoms 2 and 4. Always check the ``Substituent`` line: it
lists the atoms that were actually measured.

The same values are saved in ``molecule_2_4_sterimol.json`` (see `JSON Output`_). To keep only the report on the
terminal, turn off log messages with ``chemsmart run --no-stream sterimol ...``.

Submitting to a Compute Node
============================

Use ``chemsmart sub`` with the same options:

.. code:: bash

   chemsmart sub -s SERVER sterimol --file molecule.xyz --axis 2 4 --auto-fragment

The selection is checked before anything is submitted, so a wrong axis or atom list fails straight away. The job then
runs in the submission folder, writes the same JSON report there, and prints the values to the scheduler's output file,
for example ``molecule_2_4.slurmout``.

Your ``server.yaml`` needs no extra section, because no external program is involved. Submit from the Python environment
where CHEMSMART is installed: the job uses that same environment, so it must also be available on the compute node. Add
``--test`` after ``sub`` to write the scripts without submitting them.

*********
 Options
*********

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filename, --file``
      -  string
      -  Geometry or calculation output file (required)

   -  -  ``--axis ATOM1 ATOM2``
      -  two integers
      -  Parent atom and first substituent atom, 1-based (required). The axis points from ``ATOM1`` to ``ATOM2``

   -  -  ``--auto-fragment``
      -  flag
      -  Select the substituent by cutting the ``ATOM1``-``ATOM2`` bond

   -  -  ``--atoms``
      -  string
      -  Explicit 1-based substituent atoms, including ``ATOM2``, for example ``"6,7,12-18"``

   -  -  ``-i, --index``
      -  string
      -  1-based index of the structure to use (default: ``-1``, the last one)

   -  -  ``-l, --label``
      -  string
      -  Job label (default: the file name, the axis, and any ``--index``)

   -  -  ``-o, --output``
      -  string
      -  JSON report to write (default: ``{label}_sterimol.json``)

Use exactly one of ``--auto-fragment`` and ``--atoms``. There is no default, because measuring the wrong atoms gives
believable numbers.

********************
 Choosing the Atoms
********************

Automatic selection
===================

``--auto-fragment`` measures "the group on this side of this bond". It infers the connectivity from the distances, cuts
the ``ATOM1``-``ATOM2`` bond, and takes everything still connected to ``ATOM2``. The command stops with an error if:

-  the two atoms are not bonded in the inferred connectivity, or
-  the bond is in a ring, so cutting it does not separate anything.

Atoms connected to neither side, such as a counterion, are left out with a warning in the log.

Explicit selection
==================

Use ``--atoms`` for ring bonds, and whenever the inferred connectivity may be wrong, as it can be around metals and in
transition states:

.. code:: bash

   chemsmart run sterimol -f complex.log --axis 1 12 --atoms "12-20,34-41"

The list must include ``ATOM2`` and must not include ``ATOM1``. Every other atom is left out of the measurement.

Trajectories
============

``--index`` selects one structure, and the default is the last one:

.. code:: bash

   chemsmart run sterimol -f scan.xyz -i 3 --axis 2 4 --auto-fragment -o frame3.json

*************
 JSON Output
*************

Every run saves the values together with what is needed to reproduce and compare them:

.. code:: json

   {
     "file": "molecule.xyz",
     "index": "-1",
     "formula": "C5H10ClF",
     "axis": {"atoms": [2, 4], "symbols": ["C", "C"]},
     "selection": "auto-fragment",
     "substituent_atoms": [4, 5, 6, 7, 15, 16, 17],
     "L": 5.758,
     "B1": 1.763,
     "B5": 3.856,
     "units": "Angstrom",
     "radius_source": "RDKit periodic table (van der Waals radii)",
     "length_convention": "Measured along the axis from the first axis atom, without Verloop's +0.40 Angstrom correction"
   }

The real file holds the values at full precision. ``index`` is the option as given, and atom numbers are 1-based.

**************************
 Interpreting the Results
**************************

Read the three numbers together
===============================

No single number describes bulk, so CHEMSMART does not combine them into a score. These values were calculated with
CHEMSMART, with a hydrogen atom of R-H as ``ATOM1``, from one force-field (MMFF) geometry of each group:

.. list-table::
   :header-rows: 1
   :widths: 25 15 15 15

   -  -  Substituent
      -  L (Å)
      -  B1 (Å)
      -  B5 (Å)

   -  -  Methyl
      -  2.79
      -  1.72
      -  2.23

   -  -  Ethyl
      -  3.92
      -  1.74
      -  3.35

   -  -  Isopropyl
      -  3.90
      -  2.11
      -  3.36

   -  -  *tert*-Butyl
      -  3.87
      -  2.96
      -  3.37

   -  -  Trifluoromethyl
      -  3.08
      -  2.13
      -  2.76

   -  -  Phenyl
      -  6.16
      -  1.70
      -  3.35

Some patterns to notice:

-  **Branching shows up in B1.** From ethyl to isopropyl to *tert*-butyl, ``B5`` and ``L`` barely change, while ``B1``
   rises from 1.74 to 2.96 Å. ``B5`` alone cannot tell these groups apart, because one atom reaching outward sets its
   value.

-  **Flat groups have a small B1.** Phenyl has about the same ``B1`` as methyl, because it is only one atom thick on its
   narrow face, yet it has by far the largest ``L``. It can be large in one direction and small in another.

-  **Small differences need care.** These are lengths from one geometry. Treat a difference of a few hundredths of an Å
   as no difference, especially for ``B1``, which is sampled in 0.1° steps.

Compare like with like
======================

A value only means something next to another calculated the same way:

-  **L depends on where you measure from.** It starts at ``ATOM1``, so it depends on the kind of parent atom and on the
   bond length to it. ``B1`` and ``B5`` depend only on the direction of the axis. Use the same kind of parent atom
   across a series.

-  **Literature values usually use a different convention.** Tabulated values use hydrogen as the reference and add a
   correction of 0.40 Å to ``L``. CHEMSMART reports ``L`` without it. Adding 0.40 Å only makes sense when ``ATOM1`` is a
   hydrogen in that position, and not when it is an arbitrary parent carbon or metal.

-  **Radii differ between programs.** CHEMSMART uses RDKit's radii (H = 1.20 Å, C = 1.70 Å). morfeus, for example, uses
   H = 1.10 Å by default, so its values differ even for the same geometry. With the same radii and convention, CHEMSMART
   reproduces morfeus (see ``tests/test_sterimol_cli.py``).

Consider the conformation
=========================

The values describe the geometry you supply. A flexible group, such as *n*-butyl, can give very different ``L`` and
``B5`` values in different conformers, so a single structure may not represent it. Calculate each relevant conformer
with ``--index``, and report how the values vary. It is often more useful to treat the conformers explicitly (Brethomé
et al., 2019).

.. warning::

   ``--auto-fragment`` depends on bonds inferred from distances. For metal complexes and transition states, check the
   ``Substituent`` line, or use ``--atoms``.

******************
 Python Interface
******************

The command wraps two :class:`~chemsmart.io.molecules.structure.Molecule` methods, which also use 1-based atom numbers:

.. code:: python

   from chemsmart.io.molecules.structure import Molecule

   mol = Molecule.from_filepath("molecule.xyz")

   # All atoms on the atom 4 side of the 2-4 bond
   substituent = mol.get_substituent_atoms(atom1=2, atom2=4)

   # Everything that is neither the parent atom nor the substituent
   excluded = [i for i in range(1, mol.num_atoms + 1) if i != 2 and i not in substituent]

   result = mol.calculate_sterimol_parameters(atom1=2, atom2=4, excluded_atoms=excluded)
   print(result["L"], result["B1"], result["B5"])

``calculate_sterimol_parameters`` does not find the substituent by itself: every atom except ``atom1`` and
``excluded_atoms`` is measured. See :doc:`api-molecule-properties` for details.

************
 References
************

-  A. Verloop, W. Hoogenstraaten, J. Tipker, "Development and application of new steric substituent parameters in drug
   design", in *Drug Design*, Vol. VII, Academic Press, 1976.
-  K. C. Harper, E. N. Bess, M. S. Sigman, "Multidimensional steric parameters in the analysis of asymmetric catalytic
   reactions", *Nature Chemistry*, 2012.
-  A. V. Brethomé, S. P. Fletcher, R. S. Paton, "Conformational effects on physical-organic descriptors: the case of
   Sterimol steric parameters", *ACS Catalysis*, 2019.
