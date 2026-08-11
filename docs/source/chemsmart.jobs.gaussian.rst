chemsmart.jobs.gaussian package
###############################

The Gaussian backend translates project settings into Gaussian inputs, runs or
previews the selected job, and exposes the resulting files to the shared
analysis layer.

Core modules
============

``chemsmart.jobs.gaussian.job``
   Gaussian job objects and result-file ownership.

``chemsmart.jobs.gaussian.settings``
   Validated settings and Gaussian route construction.

``chemsmart.jobs.gaussian.writer``
   Native input materialisation from ChemSmart settings.

``chemsmart.jobs.gaussian.runner``
   Local execution, scratch handling and fake-preview behavior.

Scientific job modules
======================

Single-point, optimisation, transition-state, IRC, scan and excited-state
operations are implemented in ``singlepoint``, ``opt``, ``ts``, ``irc``,
``scan`` and ``tddft`` respectively. Additional modules cover QMMM, pKa, QRC,
RESP, NCI and other specialised Gaussian workflows.

Use :doc:`gaussian-cli-options` for the user-facing command surface and
:doc:`configuration-project-settings` for project configuration.
