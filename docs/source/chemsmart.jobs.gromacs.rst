############
GROMACS Jobs
############

ChemSmart provides workflow generation and execution support for GROMACS.
The stable interface currently covers prepared energy minimization (EM),
constant-volume equilibration (NVT), and constant-pressure equilibration (NPT)
jobs.

This page distinguishes between:

* the behavior currently implemented by ChemSmart;
* the standard GROMACS workflow described in the official documentation; and
* features that are not automated yet.

The documentation below is aligned with the GROMACS 2026.3 user guide,
command reference, and official introductory molecular-dynamics tutorial.

Workflow overview
=================

The stable prepared workflow is:

.. code-block:: text

   structure + topology + MDP
              |
              v
          gmx grompp
              |
              v
             TPR
              |
              v
          gmx mdrun

``gmx grompp`` reads the molecular topology, coordinates, simulation
parameters, and an optional index file, then produces a binary ``.tpr`` run
input file. ``gmx mdrun`` uses the ``.tpr`` file to run energy minimization or
molecular dynamics.

ChemSmart follows this pattern by:

#. accepting a prepared structure and matching topology;
#. using a user-provided MDP file or generating one when it is absent;
#. assembling the TPR with ``grompp``; and
#. executing the calculation with ``mdrun``.

Supported job types
===================

.. list-table::
   :header-rows: 1
   :widths: 15 27 28 30

   * - Job
     - CLI subcommand
     - Generated MDP
     - Intended stage
   * - EM
     - ``chemsmart run gromacs em``
     - ``em.mdp``
     - Energy minimization
   * - NVT
     - ``chemsmart run gromacs nvt``
     - ``nvt.mdp``
     - Temperature equilibration at fixed volume
   * - NPT
     - ``chemsmart run gromacs npt``
     - ``npt.mdp``
     - Pressure and density equilibration at fixed temperature

The CLI subcommand selects the concrete job class. When a YAML file is used,
keep ``project.job_type`` consistent with the selected subcommand.

Prepared workflow
=================

Required inputs
---------------

A prepared workflow requires:

``structure_file``
   A coordinate file accepted by ``gmx grompp``, normally a ``.gro`` file.

``top_file`` or ``topology_file``
   A matching GROMACS topology.

The MDP file is optional in ChemSmart. When ``mdp_file`` is absent,
``GromacsInputWriter`` writes a stage-specific default file before validation.

Optional inputs
---------------

``mdp_file``
   A custom MDP file. When supplied, ChemSmart uses it without overwriting it.

``tpr_file``
   The desired TPR path. When absent, the job supplies its default output name.

``index_file``
   An optional ``.ndx`` file passed to ``gmx grompp`` with ``-n``.

``itp_files``
   Include-topology files associated with the project.

.. important::

   The current runner passes the main topology to ``gmx grompp`` with ``-p``.
   It does not rewrite topology ``#include`` directives or automatically create
   include search paths. Every referenced ``.itp`` or position-restraint file
   must already be resolvable from the working directory, the topology include
   path, or the GROMACS library path.

Prepared command pattern
------------------------

For a prepared job, ChemSmart builds commands equivalent to:

.. code-block:: bash

   gmx grompp \
     -f stage.mdp \
     -c input.gro \
     -p topol.top \
     -o stage.tpr

   gmx mdrun -deffnm stage

When configured, ChemSmart may additionally pass:

* ``-n index.ndx`` to ``grompp``;
* ``-maxwarn N`` to ``grompp``;
* ``-nt N`` to ``mdrun``;
* ``-ntmpi N`` to ``mdrun``;
* ``-ntomp N`` to ``mdrun``; and
* values from ``mdrun_extra_args``.

``-deffnm`` sets the default filename stem for the GROMACS output files.

.. warning::

   GROMACS documents ``grompp -maxwarn`` as an option that is not intended for
   normal use and may produce unstable systems. Correct and understand warnings
   before setting ``grompp_maxwarn`` above zero.

MDP generation
==============

Behavior
--------

When ``mdp_file`` is supplied:

* the file is used directly;
* ChemSmart does not overwrite it; and
* the user is responsible for protocol validity.

When ``mdp_file`` is absent:

* EM writes ``em.mdp``;
* NVT writes ``nvt.mdp``; and
* NPT writes ``npt.mdp``.

The generated MDP files are minimal implementation defaults. They are useful
for workflow automation and testing, but they are not universal scientific
protocols. Force-field requirements, system composition, restraints, timestep,
coupling groups, equilibration duration, and production settings must be
reviewed for each system.

Generated EM defaults
---------------------

.. list-table::
   :header-rows: 1
   :widths: 34 22 44

   * - Option
     - Value
     - Meaning
   * - ``integrator``
     - ``steep``
     - Steepest-descent energy minimization
   * - ``emtol``
     - ``1000.0``
     - Maximum-force convergence threshold
   * - ``emstep``
     - ``0.01``
     - Initial minimization step in nm
   * - ``nsteps``
     - ``50000``
     - Maximum minimization steps
   * - ``cutoff-scheme``
     - ``Verlet``
     - Verlet cutoff scheme
   * - ``coulombtype``
     - ``PME``
     - Particle-mesh Ewald electrostatics
   * - ``rcoulomb`` / ``rvdw``
     - ``1.0`` / ``1.0``
     - Cutoff distances in nm
   * - ``pbc``
     - ``xyz``
     - Periodic boundaries in all dimensions

Generated NVT defaults
----------------------

.. list-table::
   :header-rows: 1
   :widths: 34 22 44

   * - Option
     - Value
     - Meaning
   * - ``integrator``
     - ``md``
     - Leap-frog molecular dynamics
   * - ``dt``
     - ``0.002``
     - Timestep in ps
   * - ``nsteps``
     - ``50000``
     - 100 ps at the default timestep
   * - ``continuation``
     - ``no``
     - Start a new constrained run
   * - ``constraints``
     - ``h-bonds``
     - Constrain bonds involving hydrogen
   * - ``constraint_algorithm``
     - ``lincs``
     - LINCS constraint solver
   * - ``tcoupl``
     - ``V-rescale``
     - Stochastic velocity-rescaling thermostat
   * - ``tc-grps``
     - ``System``
     - One temperature-coupling group
   * - ``tau_t``
     - ``0.1``
     - Temperature-coupling time constant in ps
   * - ``ref_t``
     - ``300``
     - Reference temperature in K
   * - ``pcoupl``
     - ``no``
     - No pressure coupling
   * - ``gen_vel``
     - ``yes``
     - Generate Maxwell-distributed velocities
   * - ``gen_temp``
     - ``300``
     - Velocity-generation temperature in K
   * - ``gen_seed``
     - ``-1``
     - Use a pseudo-random seed

``temperature``, ``timestep``, ``constraints``, ``constraint_algorithm``, and
``thermostat`` can be supplied through project settings. Other stage-specific
defaults remain internal to the job and writer.

Generated NPT defaults
----------------------

.. list-table::
   :header-rows: 1
   :widths: 34 22 44

   * - Option
     - Value
     - Meaning
   * - ``integrator``
     - ``md``
     - Leap-frog molecular dynamics
   * - ``dt``
     - ``0.002``
     - Timestep in ps
   * - ``nsteps``
     - ``50000``
     - 100 ps at the default timestep
   * - ``continuation``
     - ``yes``
     - Continue from an equilibrated structure
   * - ``constraints``
     - ``h-bonds``
     - Constrain bonds involving hydrogen
   * - ``constraint_algorithm``
     - ``lincs``
     - LINCS constraint solver
   * - ``tcoupl``
     - ``V-rescale``
     - Stochastic velocity-rescaling thermostat
   * - ``tc-grps``
     - ``System``
     - One temperature-coupling group
   * - ``tau_t``
     - ``0.1``
     - Temperature-coupling time constant in ps
   * - ``ref_t``
     - ``300``
     - Reference temperature in K
   * - ``pcoupl``
     - ``Parrinello-Rahman``
     - Current ChemSmart implementation default
   * - ``pcoupltype``
     - ``isotropic``
     - Isotropic pressure coupling
   * - ``tau_p``
     - ``2.0``
     - Pressure-coupling time constant in ps
   * - ``ref_p``
     - ``1.0``
     - Reference pressure in bar
   * - ``compressibility``
     - ``4.5e-5``
     - Water compressibility near 300 K in bar :sup:`-1`
   * - ``gen_vel``
     - ``no``
     - Reuse velocities from the input state

.. important::

   The current official GROMACS introductory tutorial uses ``C-rescale`` for
   NPT equilibration. The GROMACS manual states that C-rescale can be used for
   both equilibration and production, while Parrinello-Rahman can show large
   oscillations when the starting pressure is far from the target.

   ChemSmart currently generates ``Parrinello-Rahman`` unless ``barostat`` is
   overridden. This is an implementation default, not a recommendation for
   every system. Use a reviewed custom MDP file when protocol fidelity matters.

Relationship to the official tutorial
=====================================

The official introductory tutorial follows this general sequence:

.. code-block:: text

   system preparation
        |
        v
       EM
        |
        v
       NVT
        |
        v
       NPT
        |
        v
   production MD

Energy minimization
-------------------

The official example assembles an EM TPR from an MDP file, an ionized
coordinate file, and the topology, then runs ``mdrun``:

.. code-block:: bash

   gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
   gmx mdrun -deffnm em

ChemSmart's prepared EM workflow represents this pattern directly.

NVT equilibration
-----------------

The official introductory example uses the energy-minimized structure and also
passes restraint reference coordinates:

.. code-block:: bash

   gmx grompp \
     -f nvt.mdp \
     -c em.gro \
     -r em.gro \
     -p topol.top \
     -o nvt.tpr

   gmx mdrun -deffnm nvt

The official NVT example generates velocities, uses V-rescale temperature
coupling, and disables pressure coupling.

.. warning::

   When position restraints are active, GROMACS requires restraint coordinates
   through ``grompp -r``. The current ChemSmart prepared runner does not yet
   expose a dedicated restraint-reference input and does not add ``-r``.

   Therefore, a topology that activates position restraints may not work with
   the current generated command. This is a known workflow limitation, not a
   reason to bypass the resulting GROMACS error.

NPT equilibration
-----------------

The official introductory example uses the final NVT coordinates, the
restraint-reference coordinates, and the NVT checkpoint:

.. code-block:: bash

   gmx grompp \
     -f npt.mdp \
     -c nvt.gro \
     -r nvt.gro \
     -t nvt.cpt \
     -p topol.top \
     -o npt.tpr

   gmx mdrun -deffnm npt

The checkpoint preserves velocities and state variables needed for a continuous
simulation. The official tutorial sets ``continuation = yes`` and
``gen_vel = no`` for this stage.

.. warning::

   The current ChemSmart prepared runner does not yet add ``grompp -t`` or
   ``mdrun -cpi``. A prepared NPT job can read velocities stored in a ``.gro``
   file, but it does not provide exact checkpoint continuation of thermostat,
   barostat, and other state variables.

   Do not describe the current NPT command as an exact continuation of the NVT
   run. Checkpoint-aware continuation requires a future runner extension.

CLI usage
=========

Direct prepared jobs
--------------------

Generate an EM MDP automatically:

.. code-block:: bash

   chemsmart run gromacs em \
     --structure input.gro \
     --top topol.top \
     --workflow prepared

Use a custom EM MDP:

.. code-block:: bash

   chemsmart run gromacs em \
     --mdp em.mdp \
     --structure input.gro \
     --top topol.top \
     --workflow prepared

Generate an NVT MDP automatically:

.. code-block:: bash

   chemsmart run gromacs nvt \
     --structure em.gro \
     --top topol.top \
     --workflow prepared

Generate an NPT MDP automatically:

.. code-block:: bash

   chemsmart run gromacs npt \
     --structure nvt.gro \
     --top topol.top \
     --workflow prepared

Use an index file and associated ITP files:

.. code-block:: bash

   chemsmart run gromacs nvt \
     --structure em.gro \
     --top topol.top \
     --index index.ndx \
     --itp protein.itp \
     --itp posre.itp \
     --workflow prepared

The ``--itp`` option can be supplied multiple times. The topology must still
contain correct and resolvable ``#include`` statements.

YAML-driven jobs
----------------

Run a stage from a project YAML file:

.. code-block:: bash

   chemsmart run gromacs -p project.yaml em
   chemsmart run gromacs -p project.yaml nvt
   chemsmart run gromacs -p project.yaml npt

The explicit equivalent is:

.. code-block:: bash

   chemsmart run gromacs --project-yaml project.yaml npt

Project YAML
============

Structure
---------

A GROMACS project YAML can contain:

``project``
   Project metadata, workflow mode, and expected job type.

``inputs``
   Concrete input and output file paths.

``gromacs_settings``
   Reusable GROMACS settings shared by job stages.

``runtime``
   ``grompp`` warning handling and ``mdrun`` resource arguments.

Relative paths are resolved from the directory containing the YAML file.

Minimal prepared EM project
---------------------------

.. code-block:: yaml

   project:
     name: prepared_em
     type: gromacs
     job_type: em
     mode: prepared

   inputs:
     structure_file: input.gro
     topology_file: topol.top

   runtime:
     grompp_maxwarn: 0

Run it with:

.. code-block:: bash

   chemsmart run gromacs -p project.yaml em

Minimal prepared NVT project
----------------------------

.. code-block:: yaml

   project:
     name: prepared_nvt
     type: gromacs
     job_type: nvt
     mode: prepared

   inputs:
     structure_file: em.gro
     topology_file: topol.top

   gromacs_settings:
     timestep: 0.002
     temperature: 300.0
     thermostat: V-rescale
     constraints: h-bonds
     constraint_algorithm: lincs

   runtime:
     grompp_maxwarn: 0
     mdrun_ntmpi: 1
     mdrun_ntomp: 8

Run it with:

.. code-block:: bash

   chemsmart run gromacs -p project.yaml nvt

Minimal prepared NPT project
----------------------------

.. code-block:: yaml

   project:
     name: prepared_npt
     type: gromacs
     job_type: npt
     mode: prepared

   inputs:
     structure_file: nvt.gro
     topology_file: topol.top

   gromacs_settings:
     timestep: 0.002
     temperature: 300.0
     pressure: 1.0
     thermostat: V-rescale
     barostat: C-rescale
     constraints: h-bonds
     constraint_algorithm: lincs

   runtime:
     grompp_maxwarn: 0
     mdrun_ntmpi: 1
     mdrun_ntomp: 8

Run it with:

.. code-block:: bash

   chemsmart run gromacs -p project.yaml npt

The example explicitly selects ``C-rescale`` to follow the current official
introductory tutorial more closely. Omitting ``barostat`` uses the current
ChemSmart writer default.

Custom MDP project
------------------

A custom MDP file takes priority over generated defaults:

.. code-block:: yaml

   project:
     name: custom_npt
     type: gromacs
     job_type: npt
     mode: prepared

   inputs:
     mdp_file: npt-reviewed.mdp
     structure_file: nvt.gro
     topology_file: topol.top
     index_file: index.ndx
     itp_files:
       - protein.itp
       - posre.itp

   runtime:
     grompp_maxwarn: 0
     mdrun_extra_args:
       - -pin
       - "on"

Project settings reference
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 25 47

   * - Key
     - Typical section
     - Purpose
   * - ``workflow``
     - ``project.mode``
     - ``prepared`` or ``full_setup``
   * - ``job_type``
     - ``project.job_type``
     - Expected stage: ``em``, ``nvt``, or ``npt``
   * - ``mdp_file``
     - ``inputs``
     - Optional custom MDP
   * - ``ions_mdp_file``
     - ``inputs``
     - MDP used to assemble the ion-placement TPR
   * - ``structure_file``
     - ``inputs``
     - Prepared coordinate file
   * - ``input_pdb``
     - ``inputs``
     - Raw PDB or structure for ``full_setup``
   * - ``top_file`` / ``topology_file``
     - ``inputs``
     - Main topology
   * - ``tpr_file``
     - ``inputs``
     - Desired run-input output path
   * - ``index_file``
     - ``inputs``
     - Optional index file
   * - ``itp_files``
     - ``inputs``
     - Associated include-topology files
   * - ``force_field``
     - ``gromacs_settings``
     - ``pdb2gmx -ff`` value for ``full_setup``
   * - ``water_model``
     - ``gromacs_settings``
     - ``pdb2gmx -water`` value for ``full_setup``
   * - ``timestep``
     - ``gromacs_settings``
     - Generated NVT/NPT ``dt`` in ps
   * - ``temperature``
     - ``gromacs_settings``
     - Generated NVT/NPT reference temperature in K
   * - ``pressure``
     - ``gromacs_settings``
     - Generated NPT reference pressure in bar
   * - ``thermostat``
     - ``gromacs_settings``
     - Generated NVT/NPT ``tcoupl``
   * - ``barostat``
     - ``gromacs_settings``
     - Generated NPT ``pcoupl``
   * - ``constraints``
     - ``gromacs_settings``
     - Generated NVT/NPT bond constraints
   * - ``constraint_algorithm``
     - ``gromacs_settings``
     - Generated NVT/NPT constraint solver
   * - ``box_type``
     - ``gromacs_settings``
     - ``editconf -bt`` value
   * - ``box_distance``
     - ``gromacs_settings``
     - ``editconf -d`` value in nm
   * - ``solvent_file``
     - ``inputs`` or ``gromacs_settings``
     - Optional ``solvate -cs`` file
   * - ``positive_ion``
     - ``gromacs_settings``
     - ``genion -pname`` value
   * - ``negative_ion``
     - ``gromacs_settings``
     - ``genion -nname`` value
   * - ``neutral``
     - ``gromacs_settings``
     - Add ``genion -neutral`` when true
   * - ``genion_group``
     - ``gromacs_settings``
     - Group name supplied to interactive ``genion``
   * - ``grompp_maxwarn``
     - ``runtime``
     - ``grompp -maxwarn`` value
   * - ``mdrun_threads``
     - ``runtime``
     - ``mdrun -nt`` value
   * - ``mdrun_ntmpi``
     - ``runtime``
     - ``mdrun -ntmpi`` value
   * - ``mdrun_ntomp``
     - ``runtime``
     - ``mdrun -ntomp`` value
   * - ``mdrun_extra_args``
     - ``runtime``
     - Extra arguments appended to ``mdrun``

Experimental full setup
=======================

ChemSmart also contains an initial non-interactive setup path:

.. code-block:: text

   pdb2gmx
      |
   editconf
      |
   solvate
      |
   grompp for ion placement
      |
   genion
      |
   final grompp
      |
   mdrun

This maps to the standard GROMACS preparation tools:

``pdb2gmx``
   Adds hydrogens and creates GROMACS coordinates and topology information from
   a supported structure, using the selected force field and water model.

``editconf``
   Defines or changes the simulation box.

``solvate``
   Adds solvent and updates the topology when ``-p`` is used.

``genion``
   Replaces solvent molecules with monoatomic ions and can update the topology.

Current requirements
--------------------

The ``full_setup`` settings validator requires:

* ``input_pdb`` or ``structure_file``;
* ``force_field``; and
* ``water_model``.

Current defaults include:

* cubic box;
* 1.0 nm solute-to-box distance;
* positive ion ``NA``;
* negative ion ``CL``;
* neutralization enabled; and
* solvent group ``SOL``.

.. warning::

   The current full-setup runner is an initial happy-path implementation. Its
   final state and command labels are still organized around EM setup. Use
   ``full_setup`` with the EM job while this workflow remains experimental.

   Prepared EM, NVT, and NPT jobs are the stable documented interface.

Output files
============

With ``mdrun -deffnm stage``, GROMACS normally writes files using the same
``stage`` stem. Typical outputs include:

``stage.log``
   Text log.

``stage.gro``
   Final coordinates and velocities.

``stage.edr``
   Energy, temperature, pressure, and related quantities.

``stage.trr`` and/or ``stage.xtc``
   Trajectory data, depending on the MDP output settings.

``stage.cpt``
   Checkpoint state when checkpoint output is enabled.

ChemSmart considers a GROMACS job normally terminated when the corresponding
log file contains ``Finished mdrun``.

Validation and troubleshooting
==============================

Structure/topology consistency
------------------------------

The coordinates and topology must describe the same system. A common failure is
a mismatch between the number of coordinates and the molecule counts in the
topology.

Topology includes
-----------------

Use valid GROMACS ``#include`` directives and ensure the referenced files are
available. ``gmx grompp -pp processed.top`` can be useful when debugging
preprocessor and include behavior.

GROMPP warnings
---------------

Read and correct warnings rather than automatically increasing
``grompp_maxwarn``. Inspect ``mdout.mdp`` to see the parameters actually read by
``grompp``. ``gmx dump -s stage.tpr`` can be used to inspect the generated run
input.

Position restraints
-------------------

If the MDP activates position restraints, the current prepared runner lacks the
required ``-r`` reference-coordinate option. This needs workflow support rather
than a warning override.

NPT continuation
----------------

For continuous NVT-to-NPT state transfer, the official tutorial supplies the
NVT checkpoint with ``grompp -t``. Current prepared jobs do not yet model this
checkpoint input.

Scientific validation
---------------------

A successful command is not sufficient evidence that a system is equilibrated.
Evaluate at least:

* maximum force after EM;
* temperature during NVT;
* pressure and density during NPT; and
* warnings, LINCS messages, box behavior, and physical stability.

Pressure fluctuates strongly in finite molecular systems. Assess averages,
density behavior, equilibration length, and system-specific convergence rather
than expecting an instantaneously constant pressure.

Current limitations
===================

The current implementation does not yet:

* generate a complete topology for arbitrary systems in the stable prepared
  workflow;
* rewrite or stage topology include files automatically;
* pass restraint reference coordinates with ``grompp -r``;
* pass checkpoint continuation input with ``grompp -t``;
* restart ``mdrun`` with ``-cpi``;
* automatically chain EM output into NVT and NVT output into NPT;
* validate the scientific suitability of generated MDP defaults;
* select coupling groups based on system composition;
* distinguish protein, solvent, membrane, ligand, or multi-component protocols;
  or
* replace system-specific equilibration and production validation.

Official GROMACS references
===========================

The following official resources were used to align this page:

* `GROMACS introductory molecular-dynamics tutorial
  <https://tutorials.gromacs.org/docs/md-intro-tutorial.html>`_
* `GROMACS getting-started guide
  <https://manual.gromacs.org/current/user-guide/getting-started.html>`_
* `gmx grompp
  <https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html>`_
* `gmx mdrun
  <https://manual.gromacs.org/current/onlinehelp/gmx-mdrun.html>`_
* `gmx pdb2gmx
  <https://manual.gromacs.org/current/onlinehelp/gmx-pdb2gmx.html>`_
* `gmx editconf
  <https://manual.gromacs.org/current/onlinehelp/gmx-editconf.html>`_
* `gmx solvate
  <https://manual.gromacs.org/current/onlinehelp/gmx-solvate.html>`_
* `gmx genion
  <https://manual.gromacs.org/current/onlinehelp/gmx-genion.html>`_
* `GROMACS MDP options
  <https://manual.gromacs.org/current/user-guide/mdp-options.html>`_
* `GROMACS file formats
  <https://manual.gromacs.org/current/reference-manual/file-formats.html>`_
