chemsmart.jobs.orca package
###########################

The ORCA backend translates project settings into ORCA inputs, runs or safely
previews jobs, validates job-specific completion, and exposes parsed results to
the shared analysis layer.

Core modules
============

``chemsmart.jobs.orca.job``
   ORCA job objects and result-file ownership.

``chemsmart.jobs.orca.settings``
   Validated settings, canonical job options and route construction.

``chemsmart.jobs.orca.writer``
   Native input materialisation from ChemSmart settings.

``chemsmart.jobs.orca.runner``
   Local execution, scratch handling, normal-termination checks and previews.

Scientific job modules
======================

Single-point, optimisation, transition-state, IRC, NEB and excited-state
operations are implemented in ``singlepoint``, ``opt``, ``ts``, ``irc``,
``neb`` and ``td`` respectively. Additional modules cover scans, QMMM, pKa and
QRC workflows.

Use :doc:`orca-cli-options` for the user-facing command surface and
:doc:`configuration-project-settings` for project configuration.
