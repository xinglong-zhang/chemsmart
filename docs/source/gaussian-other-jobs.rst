#####################
 Other Gaussian Jobs
#####################

This page covers broken-symmetry link jobs, fragment guesses, custom user-defined calculations, and direct input file
execution.

***********
 Link Jobs
***********

Run multi-step Gaussian calculations with linked job steps. Current Link job only supports openshell singlet broken
symmetry guess as first job, followed by link job to run different job types (opt, modred, ts, scan, sp, etc).

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
      -  Guess options (default: mix). Separate multiple options with a comma, e.g. ``mix,always``.

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

To use multiple guess options, separate them with a comma:

.. code:: bash

   chemsmart sub -s SLURM gaussian -p project -f dimer.gjf -c 0 -m 1 link -j opt -g mix,always

This sets ``guess=(mix,always)`` in the route string:

.. code:: text

   # um062x def2svp stable=opt guess=(mix,always)
   ...
   # opt freq um062x def2svp geom=check guess=read
   ...

.. _gaussian-fragment-jobs:

*********************
 Fragment Guess Jobs
*********************

``frag`` generates an unrestricted DFT guess from fragments, checks wavefunction stability, and runs the selected task.
Fragment charges and signed multiplicities apply to the initial guess.

Fragment Options
================

.. list-table::
   :header-rows: 1

   -  -  Option
      -  Description

   -  -  ``-j, --jobtype``
      -  Required continuation. Choose ``opt``, ``ts``, ``sp``, ``modred``, ``scan``, ``irc``, ``td``, ``wbi``, or
         ``nci``.

   -  -  ``-st, --stable``
      -  Defaults to ``opt``.

   -  -  ``-fi, --fragment-ids``
      -  Groups of 1-based atom indices and inclusive ranges, such as ``[[1,4-10],[2,3]]``. Assign each atom once.

   -  -  ``-fc, --fragment-charges``
      -  One integer charge per group. The charges must sum to the total charge.

   -  -  ``-fm, --fragment-multiplicities``
      -  One signed multiplicity per group. A negative sign assigns the fragment's unpaired electrons to beta orbitals.
         The signed spins must agree with the total multiplicity.

Put Gaussian options such as total charge (``-c``), total multiplicity (``-m``), and solvent before ``frag``. Put the
options above after ``frag``.

Input Examples
==============

For two separated hydrogen atoms in an open-shell singlet, use total multiplicity 1 and fragment multiplicities
``[2,-2]``.

.. code:: bash

   chemsmart sub gaussian -p project -f h2.xyz -c 0 -m 1 frag -j sp \
       -fi '[[1],[2]]' -fc '[0,0]' -fm '[2,-2]'

For the triplet, use ``-m 3`` and ``-fm '[2,2]'``.

For ``.com``, ``.gjf``, and ``.inp`` files, ``frag`` reads the ``Element(Fragment=N)`` labels and charge/multiplicity
pairs from the first geometry block. Each explicit fragment option overrides its corresponding field.

.. code:: bash

   chemsmart sub gaussian -p project -f fragments.com frag -j opt

Generated Workflow
==================

#. ``Guess=(Fragment=N,Only)`` builds the fragment guess without solvent.
#. ``Stable=Opt`` checks the wavefunction by default and reoptimizes it if an instability is found.
#. The requested task runs with the selected project and CLI settings.

Both continuation steps use ``Geom=Check Guess=Read`` and the same checkpoint. They use the configured solvent settings.

******************
 Custom User Jobs
******************

Run custom calculations not built into CHEMSMART.

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
