#######################################
 Server Settings
#######################################

Configure server-specific settings for your HPC cluster.

**********************
 Server Configuration
**********************

The ``~/.chemsmart/server/`` directory contains YAML files for different HPC
clusters. Each file defines the server configuration needed to generate
submission scripts.

Example SLURM server configuration (``~/.chemsmart/server/shared.yaml``):

.. code-block:: yaml

   SERVER:
       SCHEDULER: SLURM
       QUEUE_NAME: RM-shared
       NUM_HOURS: 48
       MEM_GB: 100
       NUM_CORES: 64
       NUM_GPUS: Null
       NUM_THREADS: 64
       SUBMIT_COMMAND: sbatch
       SCRATCH_DIR: null
       USE_HOSTS: true
       EXTRA_COMMANDS: |
           export PATH=$HOME/bin/chemsmart:$PATH
           export PYTHONPATH=$HOME/bin/chemsmart:$PYTHONPATH
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
           tcsh -c "source ~/programs/g16/bsd/g16.login"
       ENVARS: |
           export SCRATCH=/tmp
           export GAUSS_EXEDIR=~/bin/g16
           export g16root=~/bin/g16
   ORCA:
       EXEFOLDER: ~/bin/orca_6_0_1
       LOCAL_RUN: False
       ENVARS: |
           export PATH=$HOME/bin/openmpi-4.1.6/build/bin:$PATH
           export LD_LIBRARY_PATH=$HOME/bin/openmpi-4.1.6/build/lib:$LD_LIBRARY_PATH

Customize this file for your specific HPC system, including:

- Queue names and resource limits
- Module loading commands
- Environment variables
- Scratch directory paths
