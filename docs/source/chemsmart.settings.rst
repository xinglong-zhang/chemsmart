chemsmart.settings package
##########################

The settings package resolves user, server and project YAML into validated
objects shared by CLI users and agents.

Configuration modules
=====================

``chemsmart.settings.user``
   User-level locations and defaults under ``~/.chemsmart``.

``chemsmart.settings.server``
   Local and scheduler resource descriptions.

``chemsmart.settings.submitters``
   SLURM, PBS/Torque, LSF and related submission adapters.

``chemsmart.settings.executable``
   Program executable discovery and selection.

``chemsmart.settings.capabilities``
   Maintained program/job capability metadata.

Program project modules
=======================

``gaussian``, ``orca``, ``pyscf`` and ``xtb`` load program-specific project
sections while preserving a common ChemSmart configuration model. Unknown or
incompatible fields fail before native input generation.

See :doc:`configuration-overview`, :doc:`configuration-server-settings` and
:doc:`configuration-project-settings` for the supported user-facing formats.
