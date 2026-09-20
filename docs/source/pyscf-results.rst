#######################
 Reading PySCF Results
#######################

This page covers the files a PySCF job produces and how to reuse them.

**************
 Result Files
**************

A finished job leaves four files sharing the job label:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   -  -  File
      -  Contents
   -  -  ``LABEL.h5``
      -  Structured results; this is the file chemsmart reads
   -  -  ``LABEL.out``
      -  PySCF's own log, for reading by eye
   -  -  ``LABEL.py``
      -  The generated driver script that was executed
   -  -  ``LABEL.err``
      -  Anything the calculation wrote to standard error

An optimization with the packaged ``test`` project produces:

.. code:: text

   molecule_opt_gas_phase.h5
   molecule_opt_gas_phase.out
   molecule_opt_gas_phase.py
   molecule_opt_gas_phase.err

The ``.out`` log is not parsed for numbers. chemsmart uses it to recognize a PySCF calculation and reads every quantity
from the ``.h5`` file beside it, so keep the two together.

.. warning::

   ``LABEL.py`` is regenerated on every run. Edit the project YAML or the command-line options instead; an edited script
   no longer matches the settings recorded with the job.

********************************
 What the Results File Contains
********************************

.. list-table::
   :header-rows: 1
   :widths: 20 80

   -  -  Group
      -  Contents
   -  -  ``spec``
      -  Settings as applied: method, basis, charge, spin, grid, solvent, engine, stages
   -  -  ``provenance``
      -  PySCF and GPU4PySCF versions, thread count, timings, settings digests
   -  -  ``status``
      -  Per-stage convergence, normal termination, and any failure with its traceback
   -  -  ``results``
      -  The computed quantities

Every job writes the SCF energies in Hartree, the geometry in Angstrom, the atomic numbers, the orbital energies and
occupations, and the detected point group. The last entry of the energy array is the final energy.

A ``hess`` stage adds the Hessian, harmonic frequencies in cm⁻¹, normal modes, reduced masses, and force constants.
Forces in Hartree/Bohr, Mulliken charges, and the dipole moment in Debye are written where the method and solvent
support them, and are absent otherwise.

``spec`` holds the applied settings, not the requested ones, including values PySCF resolves at run time such as the
number of basis functions and the dielectric constant of the chosen solvent.

******************
 Reusing a Result
******************

Pass the ``.out`` file wherever a geometry is expected:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule_opt_gas_phase.out hess

Solvent-phase single point on an optimized structure:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule_opt_gas_phase.out sp

Give the ``.out`` file, not the ``.h5``. Both are needed, but the ``.out`` is the one chemsmart accepts as a geometry
source.

******************
 Previewing a Job
******************

``--fake`` writes the driver script and placeholder results without running any chemistry:

.. code:: bash

   chemsmart run --fake pyscf -p test -f molecule.xyz opt

Preview a GPU job on a machine with no device:

.. code:: bash

   chemsmart run --fake pyscf -p test -f molecule.xyz --gpu sp

The results file reports a PySCF version of ``0.0.0-fake``, so a preview is never mistaken for a calculation.
