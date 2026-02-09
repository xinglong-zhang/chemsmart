#######################################################
 Electronic Structure Properties & Analyses (Gaussian)
#######################################################

This page covers electronic structure analysis capabilities using Gaussian, including single point calculations, excited
state properties, bond analysis, and molecular interaction studies.

*******************
 Single Point Jobs
*******************

Run single point energy calculations on optimized geometries.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] sp [SUBCMD_OPTIONS]

SP Options
==========

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description

   -  -  ``--remove-solvent/--no-remove-solvent``
      -  bool
      -  Disable solvent model (uses project settings by default)

   -  -  ``-sm, --solvent-model``
      -  string
      -  Solvent model (default: project settings)

   -  -  ``-si, --solvent-id``
      -  string
      -  Solvent ID (default: project settings)

   -  -  ``-so, --solvent-options``
      -  string
      -  Additional SCRF options (e.g., ``iterative``)

Basic Usage
===========

Standard single point calculation:

.. code:: bash

   chemsmart sub gaussian -p project -f optimized.log sp

With custom solvent settings:

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.log sp -sm cpcm -si toluene

Gas phase calculation (no solvent):

.. code:: bash

   chemsmart sub gaussian -p project -f ethanol_opt.log -c 0 -m 1 sp --remove-solvent

The output file will be named ``ethanol_opt_gas_phase.log``.

****************
 DI-AS Analysis
****************

Run Distortion-Interaction/Activation-Strain (DI-AS) analysis for reaction mechanisms.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] dias [SUBCMD_OPTIONS]

DI-AS Options
=============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-i, --fragment-indices``
      -  string
      -  Atom indices of one fragment (required)

   -  -  ``-n, --every-n-points``
      -  int
      -  Sample every nth point along IRC (default: 3)

   -  -  ``-s, --solv/--no-solv``
      -  bool
      -  Enable/disable solvent (default: disabled)

   -  -  ``-m, --mode``
      -  string
      -  Analysis mode: ``irc`` or ``ts`` (default: irc)

   -  -  ``-c1, --charge-of-fragment1``
      -  int
      -  Charge of fragment 1

   -  -  ``-m1, --multiplicity-of-fragment1``
      -  int
      -  Multiplicity of fragment 1

   -  -  ``-c2, --charge-of-fragment2``
      -  int
      -  Charge of fragment 2

   -  -  ``-m2, --multiplicity-of-fragment2``
      -  int
      -  Multiplicity of fragment 2

Basic Usage
===========

Run DI-AS analysis for atoms 5-17 at every 10th point along an IRC:

.. code:: bash

   chemsmart sub gaussian -p project -f irc.log dias -i 5-17 -n 10

.. note::

   The ``irc.log`` file should be the IRC output from the transition state to the **reactant side**.

***********
 RESP Jobs
***********

Run RESP (Restrained Electrostatic Potential) charge fitting calculations.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] resp

Basic Usage
===========

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.xyz resp

.. note::

   This creates an input file with a fixed route for RESP: ``HF/6-31+G(d) SCF=Tight Pop=MK
   IOp(6/33=2,6/41=10,6/42=17,6/50=1)``

**********
 NCI Jobs
**********

Run Non-Covalent Interaction (NCI) analysis for intermolecular interactions.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] nci

Basic Usage
===========

.. code:: bash

   chemsmart sub gaussian -p project -f complex.xyz nci

*************
 TD-DFT Jobs
*************

Run time-dependent DFT calculations for excited state properties.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] td [SUBCMD_OPTIONS]

TD-DFT Options
==============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-s, --states``
      -  string
      -  State type: ``singlets``, ``triplets``, or ``50-50`` (default: singlets)

   -  -  ``-r, --root``
      -  int
      -  State of interest (default: 1)

   -  -  ``-n, --nstates``
      -  int
      -  Number of states to solve (default: 3)

   -  -  ``-e, --eqsolv``
      -  string
      -  Equilibrium or non-equilibrium PCM solvation

Basic Usage
===========

Standard TD-DFT calculation:

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.xyz td

Calculate triplet states:

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.xyz td -s triplets -n 5

Calculate 50-50 singlet-triplet mix:

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.xyz td -s 50-50 -n 4

.. _wbi-jobs:

**********
 WBI Jobs
**********

Run Wiberg Bond Index (WBI) calculations for bond analysis.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] wbi

Basic Usage
===========

.. code:: bash

   chemsmart sub gaussian -p project -f opt.log wbi

This adds the keyword ``pop=nboread`` for NBO3.1 software analysis.
