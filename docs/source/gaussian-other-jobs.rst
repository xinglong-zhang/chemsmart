#####################
 Other Gaussian Jobs
#####################

This page covers additional Gaussian job types including multi-step link
jobs, custom user-defined calculations, and direct input file execution.

***********
 Link Jobs
***********

Run multi-step Gaussian calculations with linked job steps.

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
      -  Guess options (default: mix)

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
   #N Geom=AllCheck Guess=TCheck SCRF=Check GenChk UM062X/def2SVP Freq

******************
 Custom User Jobs
******************

Run custom calculations not built into Chemsmart.

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
