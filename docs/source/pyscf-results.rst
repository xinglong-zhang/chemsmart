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
      -  Structured results; this is the file CHEMSMART reads
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

The ``.out`` log is not parsed for numbers. Every quantity is read from the ``.h5`` file, so keep it even if you discard
the log. If the job stops while starting up, for example on an unusable basis set, the ``.out`` can be empty or missing;
the ``.h5`` and ``.err`` files record why.

``chemsmart run`` exits with a non-zero status and prints the reason when a PySCF job fails, including an SCF or
geometry optimization that did not converge. In a multi-structure run, the remaining structures still run first. A rerun
under the same label deletes the old ``.h5`` first, so a run that dies early cannot be mistaken for the earlier result.

A PySCF job stores only its final structure, so ``-i`` can select only that structure from a PySCF output. An ``irc``
job's path is stored separately, as described below, and is not selectable with ``-i``. Results that could not be
computed are listed with the reason under ``status/skipped_properties``.

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

A ``hess`` stage adds the Hessian, harmonic frequencies in cm⁻¹, normal modes, reduced masses, and force constants. With
``--scf-stability``, ``results/scf_stability`` records each stability question asked about the final reference, whether
it is stable, and its lowest Hessian eigenvalue in Hartree. A run that did not ask has no record, which is not read as
stable.

A ``ts`` job records under ``status/stages/ts`` whether the saddle-point search converged and, when a Hessian followed,
``saddle_order``, the number of imaginary modes, with the imaginary frequencies themselves. A ``saddle_order`` other
than one fails the job. Without a Hessian it is null: the structure is not known to be a transition state.

An ``irc`` job records the same two fields for its starting structure under ``status/stages/irc``, with each branch's
convergence and number of points under ``branches`` and the starting structure's position in the path as ``ts_index``.
``results/irc_positions`` holds the path in Angstrom and ``results/irc_energies`` its energies in Hartree. The starting
structure's frequencies are a check, not a result, and are not stored as ``vibrational_frequencies``.

Forces in Hartree/Bohr, Mulliken charges, the dipole moment in Debye, and ⟨S²⟩ of the final determinant are written
where the method and solvent support them, and are absent otherwise.

``spec`` holds the applied settings, not the requested ones, including values PySCF resolves at run time such as the
number of basis functions and the dielectric constant of the chosen solvent. It also records what PySCF actually
computed:

-  ``reference``: the SCF reference, for example ``RKS`` for a closed shell or ``UHF`` for a radical.

-  ``xc_applied``: the exchange-correlation functional PySCF evaluated, as LibXC component ids and coefficients, with
   the short- and long-range exact-exchange fractions, the range-separation parameter, and whether non-local correlation
   was included. ``b3lyp``, for instance, records LibXC functional 402 with 20% exact exchange.

-  ``ecp_core_electrons``: the core electrons an effective core potential replaced, by element, or null for an
   all-electron calculation.

-  ``numerics_applied``: the SCF convergence threshold and cycle limit, the DFT integration grid, and the
   density-fitting basis that were actually used, including PySCF's choices for settings left unset.

A results file that declares a different PySCF result contract (``spec/result_contract_version``) stores its quantities
under other names. CHEMSMART refuses to read it rather than report, for example, the reference energy of a correlated
calculation as its result.

******************
 Reusing a Result
******************

Pass the ``.out`` file wherever a geometry is expected:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule_opt_gas_phase.out hess

Solvent-phase single point on an optimized structure:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule_opt_gas_phase.out sp

The ``.h5`` file can also be given directly, and is all that is needed. If the calculation failed before producing a
geometry, reading either file reports the recorded failure.

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
