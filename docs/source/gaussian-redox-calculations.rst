.. _gaussian-redox-calculations:

#############################
 Gaussian Redox Calculations
#############################

This page covers **Gaussian redox job submission** — gas-phase opt+freq and solvent single-points for Ox, Red, Ref_ox,
and Ref_red.

.. note::

   Output-file analysis is backend-independent. After calculations finish, use ``chemsmart run redox analyze``. See
   :ref:`redox-calculations` for the exchange scheme, dual-level free energies, and the reference registry. The same
   analysis is also registered under ``chain``; see :doc:`chain-jobs`.

.. contents:: Table of Contents
   :local:
   :depth: 2

*************
 Quick Start
*************

``-f`` is the oxidized target. The reduced target uses the same geometry unless ``--red`` is given. Built-in ``fc_fc+``
does not bundle geometries, so ``--ref-ox`` is required. ``--ref-red`` defaults to the same geometry with charge ``ox −
n``.

.. code:: bash

   chemsmart run gaussian -p my_project -f ox.xyz -c 1 -m 2 redox \
       --ref-ox ref_ox.xyz

   chemsmart sub gaussian -p my_project -f ox.xyz -c 1 -m 2 redox \
       --red red.xyz --ref-ox ref_ox.xyz

Chain submit uses the Gaussian alias from the chain YAML:

.. code:: bash

   chemsmart sub chain -p combined -f ox.xyz -c 1 -m 2 \
       redox --program gaussian --ref-ox ref_ox.xyz

Where:

-  ``-p my_project``: Gaussian project settings (gas opt/freq theory and solv SP theory)
-  ``-f ox.xyz``: Oxidized target geometry
-  ``-c 1 -m 2``: Charge and multiplicity of Ox
-  ``--ref-ox``: Oxidized reference geometry (``--ref-red`` optional)

This runs Opt (Ox, Red) → Ref Opt → SP → Ref SP.

************************
 Job Output File Naming
************************

Sub-job labels determine output filenames. For a job with label ``mol_redox``:

.. code:: text

   mol_redox_ox_opt.log      # oxidized target gas-phase opt+freq
   mol_redox_red_opt.log     # reduced target gas-phase opt+freq
   mol_redox_RefOx_opt.log   # oxidized reference opt+freq
   mol_redox_RefRed_opt.log  # reduced reference opt+freq
   mol_redox_ox_sp.log       # oxidized target solvent SP
   mol_redox_red_sp.log      # reduced target solvent SP
   mol_redox_RefOx_sp.log    # oxidized reference solvent SP
   mol_redox_RefRed_sp.log   # reduced reference solvent SP

********************
 After Calculations
********************

Use the backend-independent analysis command (not ``gaussian redox``):

.. code:: bash

   chemsmart run redox analyze --e-ref 0.0 \
       --ox-gas mol_redox_ox_opt.log \
       --red-gas mol_redox_red_opt.log \
       --ref-ox-gas mol_redox_RefOx_opt.log \
       --ref-red-gas mol_redox_RefRed_opt.log \
       --ox-solv mol_redox_ox_sp.log \
       --red-solv mol_redox_red_sp.log \
       --ref-ox-solv mol_redox_RefOx_sp.log \
       --ref-red-solv mol_redox_RefRed_sp.log

See :ref:`redox-calculations` for option tables and the reference registry.

**********
 See Also
**********

-  :ref:`redox-calculations`
-  :ref:`orca-redox-calculations`
-  :ref:`chain-workflow-subcommands`
-  :doc:`thermochemistry-analysis`
