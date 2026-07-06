##############
 GROMACS Jobs
##############

This module provides the initial GROMACS job and runner structure for executing prepared molecular dynamics workflows in
ChemSmart.

***************
 Current Scope
***************

The current implementation focuses on prepared GROMACS workflows, where the user provides the required GROMACS input
files:

-  MDP file
-  structure file, such as GRO or PDB
-  topology file
-  optional ITP files
-  optional index file

For the prepared workflow, the runner assembles the TPR file using ``grompp`` and then executes the simulation using
``mdrun``:

.. code:: text

   grompp -> mdrun

***************
 Runner Design
***************

The GROMACS runner separates different GROMACS subcommands into individual command builders. This allows different
workflows to reuse and combine commands as needed.

Current command builders include:

-  ``_get_grompp_command()``
-  ``_get_mdrun_command()``
-  ``_get_pdb2gmx_command()``
-  ``_get_editconf_command()``
-  ``_get_solvate_command()``
-  ``_get_genion_command()``

The current stable workflow is the prepared workflow:

.. code:: text

   grompp -> mdrun

Other workflows, such as full system setup, are planned but not fully implemented yet.

**************************
 Executable Configuration
**************************

The GROMACS executable is handled through ``GromacsExecutable`` and used by the GROMACS runner. The default executable
is:

.. code:: text

   gmx

A custom executable path can also be provided to the runner. In future server-specific workflows, this can be further
connected to ChemSmart server YAML settings so that the same workflow can run on local machines, clusters, or servers.

******************
 Project Settings
******************

GROMACS project-level settings are represented by ``GromacsProjectSettings``. These settings are designed to store
reusable project information, including input files, workflow type, and simulation parameters.

The settings can be created from a dictionary:

.. code:: python

   from chemsmart.settings.gromacs import GromacsProjectSettings

   settings = GromacsProjectSettings.from_dict(
       {
           "project_name": "prepared_em",
           "workflow": "prepared",
           "mdp_file": "em.mdp",
           "structure_file": "input.gro",
           "top_file": "topol.top",
           "tpr_file": "em.tpr",
           "itp_files": ["forcefield.itp"],
       }
   )

They can also be created from a YAML file:

.. code:: python

   settings = GromacsProjectSettings.from_yaml("project.yaml")

When settings are loaded from YAML, relative input paths are resolved against the directory containing the YAML file.
This allows a GROMACS project folder to be moved or executed from a different working directory more safely.

Example YAML structure:

.. code:: yaml

   project:
     name: prepared_em
     type: gromacs
     mode: prepared

   gromacs_settings:
     force_field: user_provided
     water_model: user_provided
     timestep: 0.002
     temperature: 300.0
     pressure: 1.0
     thermostat: V-rescale
     barostat: Parrinello-Rahman
     constraints: h-bonds
     constraint_algorithm: LINCS

   inputs:
     mdp_file: em.mdp
     structure_file: input.gro
     topology_file: topol.top
     tpr_file: em.tpr
     index_file: index.ndx
     itp_files:
       - forcefield.itp
       - ligand.itp

****************************
 Job Creation from Settings
****************************

Project settings can be converted into GROMACS jobs:

.. code:: python

   from chemsmart.jobs.gromacs.job import GromacsEMJob

   job = GromacsEMJob.from_project_settings(
       settings=settings,
       molecule=None,
       jobrunner=None,
   )

This provides a direct path from project configuration to executable job objects:

.. code:: text

   project.yaml
       -> GromacsProjectSettings
       -> GromacsEMJob
       -> GromacsJobRunner

*******
 Notes
*******

The project settings do not automatically decide simulation parameters. Instead, they record user-defined settings
explicitly and make the workflow easier to reproduce and automate.

The current stable workflow is the prepared workflow, where users provide existing MDP, structure, and topology files.
Full system setup is planned but not fully implemented yet.
