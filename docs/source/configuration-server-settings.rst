#################
 Server Settings
#################

Configure server-specific settings for your HPC cluster or local machine. Server configuration files are YAML files
stored in the ``~/.chemsmart/server/`` directory that define how Chemsmart submits and executes computational chemistry
jobs.

**********************
 Server Configuration
**********************

Overview
========

The ``~/.chemsmart/server/`` directory contains YAML files for different computing environments. Each file defines the
server configuration needed to generate submission scripts. Chemsmart provides several template configurations in
``chemsmart/settings/templates/.chemsmart/server/``:

-  ``SLURM.yaml`` - For clusters using the SLURM scheduler
-  ``PBS.yaml`` - For clusters using the PBS/Torque scheduler
-  ``local.yaml`` - For local workstations without a job scheduler
-  ``small.yaml`` - Example for SLURM with specific resource limits

To use a template, copy it to your ``~/.chemsmart/server/`` directory and customize it:

.. code:: bash

   cp ~/.chemsmart/server/SLURM.yaml ~/.chemsmart/server/myserver.yaml

Then use it with: ``chemsmart sub -s myserver <other commands>``

Configuration Structure
=======================

Each server configuration file contains:

#. **SERVER** section - Defines scheduler and resource allocation settings
#. **Program-specific sections** - Configure individual programs (GAUSSIAN, ORCA, NCIPLOT, etc.)

SERVER Section
==============

The SERVER section defines the job scheduler and compute resource settings.

SCHEDULER
---------

**Type:** String

**Options:** ``SLURM``, ``PBS``, ``Null``

**Description:** The job scheduler system used by your cluster. Set to ``Null`` for local execution without a scheduler.

**Examples:**

.. code:: yaml

   SERVER:
       SCHEDULER: SLURM    # For SLURM-based clusters

   SERVER:
       SCHEDULER: PBS      # For PBS/Torque-based clusters

   SERVER:
       SCHEDULER: Null     # For local workstation

QUEUE_NAME
----------

**Type:** String or Null

**Description:** The name of the queue/partition to submit jobs to. This is scheduler-specific and depends on your
cluster configuration. Set to ``Null`` for local execution.

**Examples:**

.. code:: yaml

   QUEUE_NAME: normal        # Generic queue name
   QUEUE_NAME: RM-shared     # SLURM partition name
   QUEUE_NAME: RM-small      # SLURM partition with limited resources
   QUEUE_NAME: Null          # No queue for local execution

NUM_HOURS
---------

**Type:** Integer or Null

**Description:** Maximum wall-clock time for the job in hours. Set to ``Null`` for local execution or when not required
by the scheduler.

**Examples:**

.. code:: yaml

   NUM_HOURS: 24    # 24-hour time limit
   NUM_HOURS: 48    # 48-hour time limit for long jobs
   NUM_HOURS: 8     # 8-hour time limit for small jobs
   NUM_HOURS: Null  # No time limit for local execution

MEM_GB
------

**Type:** Integer

**Description:** Amount of memory to request in gigabytes (GB).

**Examples:**

.. code:: yaml

   MEM_GB: 400   # Request 400 GB memory
   MEM_GB: 100   # Request 100 GB memory
   MEM_GB: 48    # Request 48 GB memory for smaller jobs
   MEM_GB: 40    # Request 40 GB for local workstation

NUM_CORES
---------

**Type:** Integer

**Description:** Number of CPU cores to request for the job.

**Examples:**

.. code:: yaml

   NUM_CORES: 64   # Request 64 cores
   NUM_CORES: 12   # Request 12 cores for smaller jobs

NUM_GPUS
--------

**Type:** Integer or Null

**Description:** Number of GPUs to request. Set to ``0`` or ``Null`` if GPUs are not needed.

**Examples:**

.. code:: yaml

   NUM_GPUS: 0     # No GPUs requested
   NUM_GPUS: 1     # Request 1 GPU
   NUM_GPUS: Null  # No GPU specification

NUM_THREADS
-----------

**Type:** Integer

**Description:** Number of threads to use for parallel execution. This typically matches NUM_CORES but can be set
differently depending on the application's threading model.

**Examples:**

.. code:: yaml

   NUM_THREADS: 64   # Use 64 threads
   NUM_THREADS: 12   # Use 12 threads

SUBMIT_COMMAND
--------------

**Type:** String or Null

**Description:** The command used to submit jobs to the scheduler. Set to ``Null`` for local execution.

**Examples:**

.. code:: yaml

   SUBMIT_COMMAND: sbatch   # For SLURM
   SUBMIT_COMMAND: qsub     # For PBS/Torque
   SUBMIT_COMMAND: Null     # For local execution

PROJECT
-------

**Type:** String or Null (optional)

**Description:** Project or account number for billing/accounting on HPC systems. Comment out or set to ``Null`` if not
required.

**Examples:**

.. code:: yaml

   PROJECT: 13003611
   ##PROJECT: 13002374   # Commented out alternative project
   PROJECT: Null         # No project specification

SCRATCH_DIR
-----------

**Type:** String or Null

**Description:** Path to the scratch directory for temporary files. Set to ``null`` if not using a specific scratch
location.

**Examples:**

.. code:: yaml

   SCRATCH_DIR: /scratch/user
   SCRATCH_DIR: null

USE_HOSTS
---------

**Type:** Boolean

**Description:** Whether to use host-specific configurations. Set to ``true`` to enable host-based settings, ``false``
to disable.

**Examples:**

.. code:: yaml

   USE_HOSTS: true
   USE_HOSTS: false

EXTRA_COMMANDS
--------------

**Type:** Multiline string

**Description:** Additional shell commands to include in the job submission script. This can be used to load modules,
set environment variables, or activate conda environments. Use the pipe (``|``) character to define multiline content.

**Examples:**

.. code:: yaml

   EXTRA_COMMANDS: |
       export PATH=$HOME/bin/chemsmart:$PATH
       export PYTHONPATH=$HOME/bin/chemsmart:$PYTHONPATH
       source ~/miniconda3/etc/profile.d/conda.sh
       conda activate chemsmart

Program-Specific Sections
=========================

Each computational chemistry program (GAUSSIAN, ORCA, NCIPLOT) has its own configuration section. These sections define
program-specific paths, execution settings, and environment variables.

GAUSSIAN Section
----------------

Configuration for Gaussian quantum chemistry software.

EXEFOLDER
^^^^^^^^^

**Type:** String

**Description:** Path to the Gaussian installation directory.

**Example:**

.. code:: yaml

   GAUSSIAN:
       EXEFOLDER: ~/bin/g16

LOCAL_RUN
^^^^^^^^^

**Type:** Boolean

**Description:** Whether to run Gaussian in local/serial mode (``True``) or parallel mode (``False``). When ``True``,
uses serial execution commands; when ``False``, uses parallel execution commands.

**Example:**

.. code:: yaml

   LOCAL_RUN: True   # Use serial execution
   LOCAL_RUN: False  # Use parallel execution

SCRATCH
^^^^^^^

**Type:** Boolean

**Description:** Whether to use a scratch directory for temporary Gaussian files. When ``True``, Gaussian will write
temporary files to the scratch directory specified in ENVARS.

**Example:**

.. code:: yaml

   SCRATCH: True   # Use scratch directory
   SCRATCH: False  # Run in job directory

CONDA_ENV
^^^^^^^^^

**Type:** Multiline string

**Description:** Commands to activate the conda environment for Gaussian. Use the pipe (``|``) character for multiline
content.

**Example:**

.. code:: yaml

   CONDA_ENV: |
       source ~/miniconda3/etc/profile.d/conda.sh
       conda activate chemsmart

MODULES
^^^^^^^

**Type:** Multiline string

**Description:** Module loading commands for HPC systems. Use this to load required modules before running Gaussian.

**Example:**

.. code:: yaml

   MODULES: |
       module purge
       module load craype-x86-rome
       module load libfabric/1.11.0.4.125

SCRIPTS
^^^^^^^

**Type:** Multiline string

**Description:** Additional scripts to run before executing Gaussian, such as initialization scripts.

**Example:**

.. code:: yaml

   SCRIPTS: |
       tcsh -c "source ~/bin/g16/bsd/g16.login"

ENVARS
^^^^^^

**Type:** Multiline string

**Description:** Environment variables required by Gaussian. Essential variables include SCRATCH, GAUSS_EXEDIR, and
g16root.

**Example:**

.. code:: yaml

   ENVARS: |
       export SCRATCH=~/scratch
       export GAUSS_EXEDIR=~/bin/g16
       export g16root=~/bin/g16

ORCA Section
------------

Configuration for ORCA quantum chemistry software.

EXEFOLDER
^^^^^^^^^

**Type:** String

**Description:** Path to the ORCA installation directory.

**Example:**

.. code:: yaml

   ORCA:
       EXEFOLDER: ~/bin/orca_6_0_0

LOCAL_RUN
^^^^^^^^^

**Type:** Boolean

**Description:** Whether to run ORCA in local/serial mode (``True``) or parallel mode (``False``). ORCA typically runs
in parallel mode (``False``).

**Example:**

.. code:: yaml

   LOCAL_RUN: False  # Use parallel execution

SCRATCH
^^^^^^^

**Type:** Boolean

**Description:** Whether to use a scratch directory for temporary ORCA files.

**Example:**

.. code:: yaml

   SCRATCH: True   # Use scratch directory
   SCRATCH: False  # Run in job directory

CONDA_ENV
^^^^^^^^^

**Type:** Multiline string

**Description:** Commands to activate the conda environment for ORCA.

**Example:**

.. code:: yaml

   CONDA_ENV: |
       source ~/miniconda3/etc/profile.d/conda.sh
       conda activate ~/miniconda3/envs/chemsmart

MODULES
^^^^^^^

**Type:** Multiline string

**Description:** Module loading commands required for ORCA, typically including MPI libraries.

**Example:**

.. code:: yaml

   MODULES: |
       module purge
       module load libfabric
       module load openmpi

ENVARS
^^^^^^

**Type:** Multiline string

**Description:** Environment variables required by ORCA, including scratch directory and MPI paths.

**Example:**

.. code:: yaml

   ENVARS: |
       export SCRATCH=~/scratch
       export PATH=$HOME/bin/openmpi-4.1.6/build/bin:$PATH
       export LD_LIBRARY_PATH=$HOME/bin/openmpi-4.1.6/build/lib:$LD_LIBRARY_PATH

NCIPLOT Section
---------------

Configuration for NCIPLOT (Non-Covalent Interactions) software.

EXEFOLDER
^^^^^^^^^

**Type:** String

**Description:** Path to the NCIPLOT installation directory.

**Example:**

.. code:: yaml

   NCIPLOT:
       EXEFOLDER: ~/bin/nciplot

LOCAL_RUN
^^^^^^^^^

**Type:** Boolean

**Description:** Whether to run NCIPLOT in local/serial mode (``True``) or parallel mode (``False``).

**Example:**

.. code:: yaml

   LOCAL_RUN: False

SCRATCH
^^^^^^^

**Type:** Boolean

**Description:** Whether to use a scratch directory for temporary NCIPLOT files.

**Example:**

.. code:: yaml

   SCRATCH: True

CONDA_ENV
^^^^^^^^^

**Type:** Multiline string

**Description:** Commands to activate the conda environment for NCIPLOT.

**Example:**

.. code:: yaml

   CONDA_ENV: |
       source ~/miniconda3/etc/profile.d/conda.sh
       conda activate ~/miniconda3/envs/chemsmart

MODULES
^^^^^^^

**Type:** Multiline string

**Description:** Module loading commands for NCIPLOT.

**Example:**

.. code:: yaml

   MODULES: |
       module purge

ENVARS
^^^^^^

**Type:** Multiline string

**Description:** Environment variables required by NCIPLOT, including NCIPLOT_HOME.

**Example:**

.. code:: yaml

   ENVARS: |
       export SCRATCH=~/scratch
       export NCIPLOT_HOME=~/bin/nciplot

Complete Configuration Examples
===============================

SLURM Cluster Configuration
---------------------------

Complete example for a SLURM-based HPC cluster:

.. code:: yaml

   SERVER:
       SCHEDULER: SLURM
       QUEUE_NAME: normal
       NUM_HOURS: 24
       MEM_GB: 400
       NUM_CORES: 64
       NUM_GPUS: 0
       NUM_THREADS: 64
       SUBMIT_COMMAND: sbatch
       SCRATCH_DIR: null
       USE_HOSTS: true
       EXTRA_COMMANDS: |
           #extra commands to activate chemsmart environment in submission script
   GAUSSIAN:
       EXEFOLDER: ~/bin/g16
       LOCAL_RUN: True
       SCRATCH: True
       CONDA_ENV: |
           source ~/miniconda3/etc/profile.d/conda.sh
           conda activate chemsmart
       MODULES: |
           module purge
           module load craype-x86-rome
           module load libfabric/1.11.0.4.125
       SCRIPTS: |
           tcsh -c "source ~/bin/g16/bsd/g16.login"
       ENVARS: |
           export SCRATCH=~/scratch
           export GAUSS_EXEDIR=~/bin/g16
           export g16root=~/bin/g16
   ORCA:
       EXEFOLDER: ~/bin/orca_6_0_0
       LOCAL_RUN: False
       SCRATCH: True
       CONDA_ENV: |
           source ~/miniconda3/etc/profile.d/conda.sh
           conda activate ~/miniconda3/envs/chemsmart
       MODULES: |
           module purge
           module load libfabric
           module load openmpi
       ENVARS: |
           export SCRATCH=~/scratch
   NCIPLOT:
       EXEFOLDER: ~/bin/nciplot
       LOCAL_RUN: False
       SCRATCH: True
       CONDA_ENV: |
           source ~/miniconda3/etc/profile.d/conda.sh
           conda activate ~/miniconda3/envs/chemsmart
       MODULES: |
           module purge
       ENVARS: |
           export SCRATCH=~/scratch
           export NCIPLOT_HOME=~/bin/nciplot

PBS/Torque Cluster Configuration
--------------------------------

Complete example for a PBS/Torque-based HPC cluster:

.. code:: yaml

   SERVER:
       SCHEDULER: PBS
       QUEUE_NAME: normal
       NUM_HOURS: 24
       MEM_GB: 400
       NUM_CORES: 64
       NUM_GPUS: 0
       NUM_THREADS: 64
       SUBMIT_COMMAND: qsub
       PROJECT: 13003611
       SCRATCH_DIR: null
       USE_HOSTS: true
       EXTRA_COMMANDS: |
           #extra commands to activate chemsmart environment in submission script
   GAUSSIAN:
       EXEFOLDER: ~/bin/g16
       LOCAL_RUN: True
       SCRATCH: True
       CONDA_ENV: |
           source ~/miniconda3/etc/profile.d/conda.sh
           conda activate chemsmart
       MODULES: |
           module purge
           module load craype-x86-rome
           module load libfabric/1.11.0.4.125
       SCRIPTS: |
           tcsh -c "source ~/bin/g16/bsd/g16.login"
       ENVARS: |
           export SCRATCH=~/scratch
           export GAUSS_EXEDIR=~/bin/g16
           export g16root=~/bin/g16
   ORCA:
       EXEFOLDER: ~/bin/orca_6_0_0
       LOCAL_RUN: False
       SCRATCH: True
       CONDA_ENV: |
           source ~/miniconda3/etc/profile.d/conda.sh
           conda activate ~/miniconda3/envs/chemsmart
       MODULES: |
           module purge
           module load libfabric
           module load openmpi
       ENVARS: |
           export SCRATCH=~/scratch
   NCIPLOT:
       EXEFOLDER: ~/bin/nciplot
       LOCAL_RUN: False
       SCRATCH: True
       CONDA_ENV: |
           source ~/miniconda3/etc/profile.d/conda.sh
           conda activate ~/miniconda3/envs/chemsmart
       MODULES: |
           module purge
       ENVARS: |
           export SCRATCH=~/scratch
           export NCIPLOT_HOME=~/bin/nciplot

Local Workstation Configuration
-------------------------------

Complete example for a local workstation without a job scheduler:

.. code:: yaml

   SERVER:
       SCHEDULER: Null
       QUEUE_NAME: Null
       NUM_HOURS: Null
       MEM_GB: 40
       NUM_CORES: 12
       NUM_GPUS: 0
       NUM_THREADS: 12
       SUBMIT_COMMAND: Null
       PROJECT: Null
       SCRATCH_DIR: Null
       USE_HOSTS: True
       EXTRA_COMMANDS: |
           #extra commands to activate chemsmart environment in submission script
   GAUSSIAN:
       EXEFOLDER: ~/bin/g16
       LOCAL_RUN: True
       SCRATCH: True
       CONDA_ENV: |
           source ~/miniconda3/etc/profile.d/conda.sh
           conda activate chemsmart
       MODULES: |
           module purge
       SCRIPTS: |
           tcsh -c "source ~/bin/g16/bsd/g16.login"
       ENVARS: |
           export SCRATCH=~/scratch
           export GAUSS_EXEDIR=~/bin/g16
           export g16root=~/bin/g16
   ORCA:
       EXEFOLDER: ~/bin/orca_6_0_0
       LOCAL_RUN: False
       SCRATCH: False
       CONDA_ENV: |
           source ~/miniconda3/etc/profile.d/conda.sh
           conda activate ~/miniconda3/envs/chemsmart
       MODULES: |
           module purge
       ENVARS: |
           export SCRATCH=~/scratch
   NCIPLOT:
       EXEFOLDER: ~/bin/nciplot
       LOCAL_RUN: False
       SCRATCH: True
       CONDA_ENV: |
           source ~/miniconda3/etc/profile.d/conda.sh
           conda activate ~/miniconda3/envs/chemsmart
       MODULES: |
           module purge
       ENVARS: |
           export SCRATCH=~/scratch
           export NCIPLOT_HOME=~/bin/nciplot

Customization Tips
==================

When customizing server configuration files:

#. **Scheduler-specific settings**: Adjust SCHEDULER, QUEUE_NAME, and SUBMIT_COMMAND based on your cluster's job
   scheduler.
#. **Resource limits**: Set NUM_HOURS, MEM_GB, NUM_CORES to match your cluster's queue limits and job requirements.
#. **Module system**: Update MODULES sections to load the correct versions of libraries and tools available on your
   system.
#. **Software paths**: Update EXEFOLDER paths to point to your actual installations of Gaussian, ORCA, and NCIPLOT.
#. **Scratch directories**: Set SCRATCH environment variables to valid paths on your system. Some HPC systems provide
   node-local scratch (e.g., ``/tmp``) while others use network-attached scratch directories.
#. **Conda environments**: Adjust conda activation commands to match your conda installation path and environment names.
#. **Project accounting**: Add or remove PROJECT field based on whether your cluster requires project/account numbers
   for job submission.
#. **MPI configuration**: For ORCA, ensure the MPI library paths are correctly set in ENVARS to match your system's MPI
   installation.

Using Custom Server Configurations
==================================

After creating a custom server configuration file:

#. Save it in ``~/.chemsmart/server/`` with a descriptive name (e.g., ``myserver.yaml``)

#. Use it with Chemsmart commands via the ``-s`` flag:

   .. code:: bash

      chemsmart sub -s myserver -g opt -i input.xyz

#. Verify the generated submission script before submitting to ensure all paths and settings are correct

#. Test with a small job first to validate the configuration works correctly on your system
