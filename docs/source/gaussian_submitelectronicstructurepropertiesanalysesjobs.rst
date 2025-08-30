Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

########################################################
 Submit Electronic Structure Properties & Analyses Jobs
########################################################

ChemSmart provides comprehensive electronic structure analysis capabilities using Gaussian. This section covers single
point calculations, excited state properties, bond analysis, and molecular interaction studies.

*******************
 Single Point Jobs
*******************

Run single point energy calculations on optimized geometries.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] sp [SUBCMD_OPTIONS]

SP-Specific OPTIONS
===================

Please set solvent parameters at the end when needed.

.. list-table:: Solvent Model Options
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description
<<<<<<< HEAD

   -  -  ``--remove-solvent/--no-remove-solvent``
      -  bool
      -  Whether to use solvent model in the job. Defaults to project settings

   -  -  ``-sm, --solvent-model``
      -  string
      -  Solvent model to be used for single point (default = none, use project settings)

   -  -  ``-si, --solvent-id``
      -  string
      -  Solvent ID to be used for single point (default = none, use project settings)

   -  -  ``-so, --solvent-options``
      -  string
      -  Additional solvent options in scrf=() route. default=None. E.g., iterative in scrf=(smd,water,iterative)

Basic Usage
===========
=======

   -  -  ``--remove-solvent/--no-remove-solvent``
      -  bool
      -  Whether to use solvent model in the job. Defaults to project settings

   -  -  ``-sm, --solvent-model``
      -  string
      -  Solvent model to be used for single point (default = none, use project settings)

   -  -  ``-si, --solvent-id``
      -  string
      -  Solvent ID to be used for single point (default = none, use project settings)

   -  -  ``-so, --solvent-options``
      -  string
      -  Additional solvent options in scrf=() route. default=None. E.g., iterative in scrf=(smd,water,iterative)

SP Basic Usage
==============
>>>>>>> 8a3e0a1 (Docs setup (#263))

**Basic single point calculation**:

   .. code:: console

      chemsmart sub -s shared gaussian -p sp_calc -f optimized.log sp

**Single point with solvent settings**:

   .. code:: console

      chemsmart sub -s shared gaussian -p sp_solvent -f molecule.log sp -so iterative

-  For single-point job that user wants to test which uses different solvent model and id from that specified in
   ``<project>``, one can do:

   .. code:: console

      chemsmart sub -s shared gaussian -p sp_solvent -f molecule.log sp -sm cpcm -si toluene

   to specify a different solvent model ``<user_solvent_model>`` and solvent ``<user_solvent_id>``.

<<<<<<< HEAD
Examples
========
=======
SP Examples
===========
>>>>>>> 8a3e0a1 (Docs setup (#263))

**Run sp job in gas**

-  To run sp job in gas, one can do:

   .. code:: console

      chemsmart sub -s cu gaussian -p project1 -f ethanol_opt.log -c 0 -m 1 sp --remove-solvent

   output file will be named as ``ethanol_opt_gas_phase.log`` and use no solvent model. Except from this command, the
   output file will be named as ``ethanol_opt_<solv model>_<solv id>.log``

****************
 DI-AS Analysis
****************

Run Distortion-Interaction/Activation-Strain analysis for reaction mechanisms.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] dias [SUBCMD_OPTIONS]

DI-AS Specific OPTIONS
======================

.. list-table:: DI-AS Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description
<<<<<<< HEAD

   -  -  ``-i, --fragment-indices``
      -  string
      -  Indices of one fragment for DI-AS analysis (required)

   -  -  ``-n, --every-n-points``
      -  int
      -  Every nth points along the IRC file to prepare for DI-AS analysis (default=3)

   -  -  ``-s, --solv/--no-solv``
      -  bool
      -  Turn on/off solvent for DI-AS job calculations (default=False)

   -  -  ``-m, --mode``
      -  string
      -  Mode of DI-AS analysis. Options: irc, ts (default=irc)

   -  -  ``-c1, --charge-of-fragment1``
      -  int
      -  Charge of fragment 1 (default=None)

   -  -  ``-m1, --multiplicity-of-fragment1``
      -  int
      -  Multiplicity of fragment 1 (default=None)

   -  -  ``-c2, --charge-of-fragment2``
      -  int
      -  Charge of fragment 2 (default=None)

   -  -  ``-m2, --multiplicity-of-fragment2``
      -  int
      -  Multiplicity of fragment 2 (default=None)

Basic Usage
===========
=======

   -  -  ``-i, --fragment-indices``
      -  string
      -  Indices of one fragment for DI-AS analysis (required)

   -  -  ``-n, --every-n-points``
      -  int
      -  Every nth points along the IRC file to prepare for DI-AS analysis (default=3)

   -  -  ``-s, --solv/--no-solv``
      -  bool
      -  Turn on/off solvent for DI-AS job calculations (default=False)

   -  -  ``-m, --mode``
      -  string
      -  Mode of DI-AS analysis. Options: irc, ts (default=irc)

   -  -  ``-c1, --charge-of-fragment1``
      -  int
      -  Charge of fragment 1 (default=None)

   -  -  ``-m1, --multiplicity-of-fragment1``
      -  int
      -  Multiplicity of fragment 1 (default=None)

   -  -  ``-c2, --charge-of-fragment2``
      -  int
      -  Charge of fragment 2 (default=None)

   -  -  ``-m2, --multiplicity-of-fragment2``
      -  int
      -  Multiplicity of fragment 2 (default=None)

DI-AS Basic Usage
=================
>>>>>>> 8a3e0a1 (Docs setup (#263))

**Basic DI-AS analysis on IRC**

-  For example to run DI-AS job for fragment 1 with atoms numbered from 5-17 at every 10 steps along the irc.log file:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f irc.log dias -i 5-17 -n 10

***********
 RESP Jobs
***********

Run RESP (Restrained Electrostatic Potential) charge fitting calculations.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] resp

<<<<<<< HEAD
Basic Usage
===========
=======
RESP Basic Usage
================
>>>>>>> 8a3e0a1 (Docs setup (#263))

-  **Basic RESP calculation**

      .. code:: console

         chemsmart sub -s shared gaussian -p resp -f molecule.xyz resp

.. note::

   This command will create an input file with fixed route for RESP job: ``HF/6-31+G(d) SCF=Tight Pop=MK
   IOp(6/33=2,6/41=10,6/42=17,6/50=1)``

   .. code:: console

      ------------------------------------------------------------------
      # HF/6-31+G(d) SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=17,6/50=1)
      ------------------------------------------------------------------

**********
 NCI Jobs
**********

Run Non-Covalent Interaction analysis for intermolecular interactions.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] nci

<<<<<<< HEAD
Basic Usage
===========
=======
NCI Basic Usage
===============
>>>>>>> 8a3e0a1 (Docs setup (#263))

-  **Basic NCI calculation**:

      .. code:: console

         chemsmart sub gaussian -p nci_analysis -f complex.xyz nci

*************
 TD-DFT Jobs
*************

Run time-dependent DFT calculations for excited state properties.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] td [SUBCMD_OPTIONS]

TD-DFT Specific Options
=======================

.. list-table:: TD-DFT Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description
<<<<<<< HEAD

   -  -  ``-s, --states``
      -  string
      -  States for closed-shell singlet systems. Options: 'singlets', 'triplets', '50-50' (default=singlets)

   -  -  ``-r, --root``
      -  int
      -  Specifies the "state of interest". The default is the first excited state (N=1) (default=1)

   -  -  ``-n, --nstates``
      -  int
      -  Solve for M states. If 50-50, this gives the number of each type of state to solve (default=3)

   -  -  ``-e, --eqsolv``
      -  string
      -  Whether to perform equilibrium or non-equilibrium PCM solvation (default=None)

Basic Usage
===========
=======

   -  -  ``-s, --states``
      -  string
      -  States for closed-shell singlet systems. Options: 'singlets', 'triplets', '50-50' (default=singlets)

   -  -  ``-r, --root``
      -  int
      -  Specifies the "state of interest". The default is the first excited state (N=1) (default=1)

   -  -  ``-n, --nstates``
      -  int
      -  Solve for M states. If 50-50, this gives the number of each type of state to solve (default=3)

   -  -  ``-e, --eqsolv``
      -  string
      -  Whether to perform equilibrium or non-equilibrium PCM solvation (default=None)

TD-DFT Basic Usage
==================
>>>>>>> 8a3e0a1 (Docs setup (#263))

**Basic TD-DFT calculation**

   .. code:: console

      chemsmart sub gaussian -p td_calc -f molecule.xyz td

**TD-DFT with specific states**

   .. code:: console

      chemsmart sub gaussian -p td_triplets -f molecule.xyz td -s triplets -n 5

**TD-DFT with 50-50 singlet-triplet mix**

   .. code:: console

      chemsmart sub gaussian -p td_mixed -f molecule.xyz td -s 50-50 -n 4

**********
 WBI Jobs
**********

Run Wiberg Bond Index calculations for bond analysis.

<<<<<<< HEAD
Basic Usage
===========
=======
WBI Basic Usage
===============
>>>>>>> 8a3e0a1 (Docs setup (#263))

**Basic WBI calculation**

   .. code:: console

      chemsmart sub gaussian -p wbi_analysis -f opt.log wbi

<<<<<<< HEAD
Examples
========
=======
WBI Examples
============
>>>>>>> 8a3e0a1 (Docs setup (#263))

**Using wbi command for NBO analysis**

   .. code:: console

      chemsmart sub -s cu gaussian -p project -f complex1_opt.log wbi

   Add keyword ``pop=nboread`` to gaussian for NBO3.1 software analysis.
