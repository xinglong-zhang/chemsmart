Installation on an HPC Cluster
##############################

CHEMSMART uses the same project YAML and program commands for local execution
and scheduler submission. Install the controller on a login or workflow node,
then describe each compute environment in a server YAML file.

For a standalone x86 Ubuntu CPU machine, start with
:doc:`installation-ubuntu-cpu-server`. This page covers scheduler-specific
configuration after the controller and chemistry programs are available.

Prerequisites
=============

Confirm with the cluster administrator:

- scheduler type and submission command;
- account or project name, partition or queue, and wall-time limits;
- available CPU, GPU, memory, scratch and filesystem policy;
- environment modules for Gaussian, ORCA, MPI, xTB, NCIPLOT or Python;
- whether outbound network access is allowed from login and compute nodes; and
- where persistent program results may be stored.

Do not copy local workstation paths or module names into cluster configuration.

Controller installation
=======================

Follow :doc:`installation-linux-macos` using a user-owned Conda environment.
On a fresh cluster clone:

.. code-block:: bash

   git clone https://github.com/Hongjiseung-ROK/chemsmart.git
   cd chemsmart
   make env
   conda activate chemsmart
   make install
   make configure

If compute nodes cannot access the internet, build the environment and package
cache on an allowed node or use the cluster's supported offline installation
workflow.

Scheduler configuration
=======================

Create a named file under ``~/.chemsmart/server``. Example SLURM shape:

.. code-block:: yaml

   SERVER:
     SCHEDULER: SLURM
     QUEUE_NAME: compute
     NUM_HOURS: 2
     MEM_GB: 32
     NUM_CORES: 8
     NUM_GPUS: 0
     NUM_THREADS: 8
     SUBMIT_COMMAND: sbatch
     PROJECT: my_account
     SCRATCH_DIR: /scratch/USER/chemsmart
     USE_HOSTS: false
     EXTRA_COMMANDS: |
       ulimit -s unlimited

   ORCA:
     EXEFOLDER: /opt/orca
     LOCAL_RUN: false
     SCRATCH: true
     MODULES: |
       module purge
       module load openmpi
     ENVARS: |
       export PATH=/opt/openmpi/bin:$PATH
       export LD_LIBRARY_PATH=/opt/openmpi/lib:$LD_LIBRARY_PATH

Replace every value with cluster facts. For PBS/Torque or LSF, use the matching
``SCHEDULER`` and submission command. See
:doc:`configuration-server-settings` for all fields.

Submission qualification
========================

Generate the script without submitting:

.. code-block:: bash

   chemsmart sub --test --fake -s my_cluster orca \
     -p test -f examples/xtb/water.xyz sp

Inspect modules, environment activation, scratch path, core count, memory,
wall time, account and generated CHEMSMART command. Then submit one small
single-point calculation and verify scheduler exit state, program normal
termination, parsed quantity and result location before qualifying advanced
workflows.

Never use an agent plan or a successful script render as evidence that a
cluster job ran. Scheduler submission, engine completion and scientific
validation are separate observations.
