.. _orca-redox-calculations:

#########################
 ORCA Redox Calculations
#########################

This page covers **ORCA redox job submission**. The command structure mirrors Gaussian redox for a consistent
experience.

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

``-f`` is the oxidized target. Built-in ``fc_fc+`` does not bundle geometries, so ``--ref-ox`` is required.
``--ref-red`` defaults to the same geometry with charge ``ox − n``.

.. code:: bash

   chemsmart run orca -p my_project -f ox.xyz -c 1 -m 2 redox \
       --ref-ox ref_ox.xyz

Chain submit uses the ORCA alias from the chain YAML:

.. code:: bash

   chemsmart sub chain -p combined -f ox.xyz -c 1 -m 2 \
       redox --program orca --ref-ox ref_ox.xyz

Where:

-  ``-p my_project``: ORCA project settings
-  ``-f ox.xyz``: Oxidized target geometry
-  ``-c 1 -m 2``: Charge and multiplicity of Ox
-  ``--ref-ox``: Oxidized reference geometry (``--ref-red`` optional)

This runs Opt (Ox, Red) → Ref Opt → SP → Ref SP.

************************
 Job Output File Naming
************************

For a job with label ``mol_redox``:

.. code:: text

   mol_redox_ox_opt.out
   mol_redox_red_opt.out
   mol_redox_RefOx_opt.out
   mol_redox_RefRed_opt.out
   mol_redox_ox_sp.out
   mol_redox_red_sp.out
   mol_redox_RefOx_sp.out
   mol_redox_RefRed_sp.out

********************
 After Calculations
********************

.. code:: bash

   chemsmart run redox analyze --e-ref 0.0 \
       --ox-gas mol_redox_ox_opt.out \
       --red-gas mol_redox_red_opt.out \
       --ref-ox-gas mol_redox_RefOx_opt.out \
       --ref-red-gas mol_redox_RefRed_opt.out \
       --ox-solv mol_redox_ox_sp.out \
       --red-solv mol_redox_red_sp.out \
       --ref-ox-solv mol_redox_RefOx_sp.out \
       --ref-red-solv mol_redox_RefRed_sp.out

ORCA ``.out`` and Gaussian ``.log`` files can be mixed. See :ref:`redox-calculations`.

**********
 See Also
**********

-  :ref:`redox-calculations`
-  :ref:`gaussian-redox-calculations`
-  :ref:`chain-workflow-subcommands`
-  :doc:`thermochemistry-analysis`
