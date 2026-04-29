#####################
 Other Gaussian Jobs
#####################

This page covers additional Gaussian job types including multi-step link jobs, custom user-defined calculations, and
direct input file execution.

****************************************
 Minimum Energy Cross Point (MECP) Jobs
****************************************

A Minimum Energy Crossing Point (MECP) is a geometry where two potential energy surfaces of different spin multiplicity
are degenerate. It is the spin-forbidden analogue of a transition state and is relevant to intersystem crossing,
spin-state reactivity, and organometallic reaction mechanisms.

ChemSmart provides a **self-contained MECP optimizer** that drives the search entirely in Python, calling Gaussian only
for single-point energies and Cartesian forces at each step. No external MECP code (e.g. MECP2, Harvey's code) is
required.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] mecp [SUBCMD_OPTIONS]

Algorithm
=========

At each optimization step the following operations are performed:

#. **State A SP+forces** — Gaussian single-point energy and Cartesian forces for spin state A at the current geometry.

#. **State B SP+forces** — Same for spin state B.

#. **Effective displacement** — The displacement vector is computed using the penalty-function / projected-gradient
   approach:

   -  Let :math:`\Delta E = E_A - E_B` and :math:`\mathbf{g}_\Delta = \nabla E_A - \nabla E_B`.

   -  **Seam correction** (drives geometry toward the crossing seam): :math:`\mathbf{d}_\text{seam} = -\frac{\Delta
      E}{\|\mathbf{g}_\Delta\|^2}\,\mathbf{g}_\Delta`

   -  **Seam-tangent gradient** (component of :math:`\nabla E_A` perpendicular to :math:`\mathbf{g}_\Delta`):
      :math:`\mathbf{g}_\perp = \nabla E_A - \frac{\nabla E_A \cdot
      \mathbf{g}_\Delta}{\|\mathbf{g}_\Delta\|^2}\,\mathbf{g}_\Delta`

   -  **Downhill step** (minimizes energy along the seam): :math:`\mathbf{d}_\text{seam\text{-}min} =
      -\alpha\,\mathbf{g}_\perp` where :math:`\alpha` is ``step_size`` (Bohr²/Hartree).

   -  **Total displacement**: :math:`\mathbf{d} = \mathbf{d}_\text{seam} + \mathbf{d}_\text{seam\text{-}min}`

#. **Trust radius** — Each atom's displacement vector is independently scaled so that its Cartesian norm does not exceed
   ``trust_radius`` (Bohr).

#. **Convergence check** — See the convergence criteria table below.

#. **Geometry update** — :math:`\mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{d}`.

MECP Options
============

.. list-table::
   :header-rows: 1
   :widths: 30 10 20 40

   -  -  Option
      -  Type
      -  Default
      -  Description

   -  -  ``--multiplicity-a``
      -  int
      -  1 (singlet)
      -  Spin multiplicity of state A. Falls back to ``-m`` if not set.

   -  -  ``--multiplicity-b``
      -  int
      -  multiplicity-a + 2
      -  Spin multiplicity of state B.

   -  -  ``--charge-a``
      -  int
      -  0
      -  Charge of state A. Falls back to ``-c`` if not set.

   -  -  ``--charge-b``
      -  int
      -  charge-a
      -  Charge of state B.

   -  -  ``--title-a``
      -  string
      -  ``"First"``
      -  Title prefix for state A Gaussian jobs.

   -  -  ``--title-b``
      -  string
      -  ``"Second"``
      -  Title prefix for state B Gaussian jobs.

   -  -  ``--max-steps``
      -  int
      -  500
      -  Maximum number of MECP optimization steps.

   -  -  ``--step-size``
      -  float
      -  0.1 Bohr²/Hartree
      -  Scaling factor :math:`\alpha` for the seam-tangent gradient step.

   -  -  ``--trust-radius``
      -  float
      -  0.3 Bohr
      -  Maximum per-atom Cartesian displacement magnitude (Bohr) per step.

   -  -  ``--energy-diff-tol``
      -  float
      -  5.0×10⁻⁵ Hartree
      -  Convergence threshold for :math:`|E_A - E_B|`.

   -  -  ``--force-max-tol``
      -  float
      -  1.5×10⁻⁵ Hartree/Bohr
      -  Convergence threshold for the maximum element of :math:`|\mathbf{g}_\perp|`.

   -  -  ``--force-rms-tol``
      -  float
      -  5.0×10⁻⁴ Hartree/Bohr
      -  Convergence threshold for the RMS of :math:`\mathbf{g}_\perp`.

   -  -  ``--disp-max-tol``
      -  float
      -  1.2×10⁻⁴ Bohr
      -  Convergence threshold for the maximum element of :math:`|\mathbf{d}|`.

   -  -  ``--disp-rms-tol``
      -  float
      -  3.8×10⁻³ Bohr
      -  Convergence threshold for the RMS of :math:`\mathbf{d}`.

Convergence Criteria
====================

Convergence is declared when **all** five criteria are simultaneously satisfied:

.. list-table::
   :header-rows: 1
   :widths: 35 20 15 30

   -  -  Quantity
      -  Symbol
      -  Default
      -  Unit

   -  -  Energy difference
      -  :math:`|E_A - E_B|`
      -  5.0×10⁻⁵
      -  Hartree

   -  -  Maximum seam-tangent gradient
      -  :math:`\max|\mathbf{g}_\perp|`
      -  1.5×10⁻⁵
      -  Hartree/Bohr

   -  -  RMS seam-tangent gradient
      -  :math:`\text{RMS}(\mathbf{g}_\perp)`
      -  5.0×10⁻⁴
      -  Hartree/Bohr

   -  -  Maximum displacement
      -  :math:`\max|\mathbf{d}|`
      -  1.2×10⁻⁴
      -  Bohr

   -  -  RMS displacement
      -  :math:`\text{RMS}(\mathbf{d})`
      -  3.8×10⁻³
      -  Bohr

Output Files
============

Two output files are produced alongside the Gaussian sub-job input/output files:

``<label>_report.log``
   Step-by-step optimization log. The file header records the run settings; each subsequent line reports:

   .. code::

      step=NNN E_A=<Hartree> E_B=<Hartree> dE=<Hartree>
      grad_max=<H/Bohr> grad_rms=<H/Bohr> disp_max=<Bohr> disp_rms=<Bohr>

   The final line reads ``Converged at step NNN.`` on successful convergence. The presence of this ``Converged`` marker
   is used by ``skip_completed`` to avoid re-running a finished job.

``<label>_traj.xyz``
   Multi-frame XYZ trajectory of the MECP geometry at every optimization step (coordinates in Ångström).

Basic Usage
===========

Singlet/triplet MECP from a Gaussian output structure (charge and multiplicity inferred from file):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log mecp

Set spin states explicitly (singlet ↔ triplet):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp --multiplicity-a 1 --multiplicity-b 3

Doublet/quartet MECP for an open-shell cation:

.. code:: bash

   chemsmart sub gaussian -p project -f radical.log -c 1 -m 2 mecp --multiplicity-a 2 --multiplicity-b 4

Tighten convergence and limit the number of steps:

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp \
       --max-steps 200 --step-size 0.05 --trust-radius 0.1 \
       --energy-diff-tol 1.0e-5

Adjust optimization step parameters only:

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp --max-steps 100 --step-size 0.05 --trust-radius 0.1

.. note::

   Each MECP step generates two Gaussian sub-jobs named ``<label>_step<NNN>_A`` and ``<label>_step<NNN>_B``
   (single-point energy + forces). These sub-jobs are always re-run (``skip_completed=False``), while the outer MECP job
   itself honours ``skip_completed`` via the ``Converged`` marker in the report file.

***********
 Link Jobs
***********

Run multi-step Gaussian calculations with linked job steps.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] link [SUBCMD_OPTIONS]

Link Options
============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  string
      -  Job type: opt, ts, modred, scan, sp, irc

   -  -  ``-st, --stable``
      -  string
      -  Stability test options (default: opt)

   -  -  ``-g, --guess``
      -  string
      -  Guess options (default: mix)

   -  -  ``--route``
      -  string
      -  Route for link section

Basic Usage
===========

Link job with optimization:

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.xyz link -j opt

Link job with single point:

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.xyz -c 0 -m 1 -r scf=qc link -j sp -so iterative

Examples
========

Optimization of singlet open-shell structure:

.. code:: bash

   chemsmart sub -s SLURM gaussian -p project -f dimer.gjf -c 0 -m 1 link -j opt

This creates a multi-step workflow:

.. code:: text

   # um062x def2svp stable=opt guess=mix
   ...
   # opt freq um062x def2svp geom=check guess=read
   ...
   #N Geom=AllCheck Guess=TCheck SCRF=Check GenChk UM062X/def2SVP Freq

******************
 Custom User Jobs
******************

Run custom calculations not built into Chemsmart.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] userjob [SUBCMD_OPTIONS]

Custom Job Options
==================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-r, --route``
      -  string
      -  User-defined route (required)

   -  -  ``-a, --append-info``
      -  string
      -  Information to append after coordinates

Basic Usage
===========

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.com -l custom_job userjob -r 'opt freq b3lyp/6-31g*' -a 'B 1 2 F'

*****************************
 Direct Input File Execution
*****************************

Run a pre-prepared Gaussian input file without modifications.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] com

Basic Usage
===========

Run a ``.com`` file:

.. code:: bash

   chemsmart sub gaussian -p project -f input_file.com com

Run a ``.gjf`` file:

.. code:: bash

   chemsmart sub gaussian -p project -f input_file.gjf com

Modify charge and multiplicity:

.. code:: bash

   chemsmart sub gaussian -p project -f input_file.com -c 1 -m 2 com
