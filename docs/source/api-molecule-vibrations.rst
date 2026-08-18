Vibrational Data
================

Vibrational properties are typically populated when reading frequency calculation output files.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   -  -  Property
      -  Description

   -  -  :attr:`~Molecule.vibrational_frequencies`
      -  Vibrational frequencies in cm⁻¹

   -  -  :attr:`~Molecule.vibrational_reduced_masses`
      -  Reduced masses per mode in amu

   -  -  :attr:`~Molecule.vibrational_force_constants`
      -  Force constants per mode in mDyne/Å

   -  -  :attr:`~Molecule.vibrational_ir_intensities`
      -  IR intensities in km/mol

   -  -  :attr:`~Molecule.vibrational_mode_symmetries`
      -  Irreducible representation labels (e.g. ``"A1"``)

   -  -  :attr:`~Molecule.vibrational_modes`
      -  Mass-weighted normal-mode displacements, each ``(N, 3)``

.. code:: python

   # Load a molecule from a frequency calculation output (co2.log has vibrational data)
   mol = Molecule.from_filepath("tests/data/GaussianTests/outputs/co2.log")

   if mol.has_vibrations:
       print(mol.num_vib_frequencies)        # 4
       print(mol.vibrational_frequencies)     # [653.74, 653.74, 1389.03, 2472.77]

       # Generate a geometry displaced along mode 1
       displaced = mol.vibrationally_displaced(mode_idx=1, amp=0.1)

       # Generate a 20-frame trajectory for animation
       frames = mol.vibrationally_displaced(mode_idx=1, amp=0.1, nframes=20)

**Attributes:**

-  ``vibrational_frequencies`` (list) — Vibrational frequencies in cm⁻¹.
-  ``vibrational_reduced_masses`` (list) — Reduced masses per mode in amu.
-  ``vibrational_force_constants`` (list) — Force constants per mode in mDyne/Å.
-  ``vibrational_ir_intensities`` (list) — IR intensities in km/mol.
-  ``vibrational_mode_symmetries`` (list) — Irreducible representation labels (e.g. ``"A1"``).
-  ``vibrational_modes`` (list) — Each entry is an ``(N, 3)`` mass-weighted normal-mode displacement.
-  ``has_vibrations`` (bool) — ``True`` if any vibrational frequencies are present.
-  ``num_vib_frequencies`` (int) — Number of vibrational frequencies available.
-  ``num_vib_modes`` (int) — Number of normal modes available.

:meth:`~chemsmart.io.molecules.structure.Molecule.vibrationally_displaced(mode_idx, amp=0.1, nframes=None)`
Generate a geometry displaced along a vibrational normal mode.

- ``mode_idx`` (int) — 1-based index of the vibrational mode to displace along.
- ``amp`` (float) — Displacement amplitude (in Ångströms).
- ``nframes`` (int or None) — If given, generates a trajectory of ``nframes`` frames from 0 to ``amp``.
- Returns: :class:`Molecule` or ``list[Molecule]`` — A single displaced molecule, or a list of trajectory frames.