Thank you for downloading and using chemsmart. Please confirm the operating system you are using and follow the
corresponding configuration steps.

If you encounter any problems, please feel free to contact us.

#######################################
 Set up the Server Settings on Cluster
#######################################

Once the user has completed the configuration of chemsmart, the server configuration on the cluster can be customized by
modifying corresponding files.

-  The ``~/.chemsmart/server/`` directory contains files related to server setup for a particular HPC cluster that the
   user is using. For example, we can specify a SLURM based server setting as ``~/.chemsmart/server/shared.yaml`` with
   the following information:

      .. code:: console

         SERVER:
             SCHEDULER: SLURM
             QUEUE_NAME: RM-shared
             NUM_HOURS: 48
             MEM_GB: 100
             NUM_CORES: 64
             NUM_GPUS: Null
             NUM_THREADS: 64
             SUBMIT_COMMAND: sbatch
             ##PROJECT: 13003611
             ##PROJECT: 13002374
             SCRATCH_DIR: null
             USE_HOSTS: true
             EXTRA_COMMANDS: |
                 export PATH=$HOME/bin/chemsmart:$PATH
                 export PATH=$HOME/bin/chemsmart/chemsmart/cli:$PATH
                 export PATH=$HOME/bin/chemsmart/chemsmart/scripts:$PATH
                 export PYTHONPATH=$HOME/bin/chemsmart:$PYTHONPATH
         GAUSSIAN:
             EXEFOLDER: ~/bin/g16
             LOCAL_RUN: True
             SCRATCH: True  # set scratch to True to run in scratch folder
             CONDA_ENV: |   # program-specific conda env
                 source ~/miniconda3/etc/profile.d/conda.sh
                 conda activate chemsmart
             MODULES: |
                 module purge
                 # module load craype-x86-rome
                 # module load libfabric/1.11.0.4.125
             SCRIPTS: |
                 tcsh -c "source ~/programs/g16/bsd/g16.login"
             ENVARS: |
                 export SCRATCH=/tmp # required if scratch is true
                 export GAUSS_EXEDIR=~/bin/g16
                 export g16root=~/bin/g16
         ORCA:
             EXEFOLDER: ~/bin/orca_6_0_1
             LOCAL_RUN: False
             ENVARS: |
                 export PATH=$HOME/bin/openmpi-4.1.6/build/bin:$PATH
                 export LD_LIBRARY_PATH=$HOME/bin/openmpi-4.1.6/build/lib:$LD_LIBRARY_PATH

   This file can be customized by user for different submission systems. This file contains the server configuration
   information that is needed for chemsmart to automatically write the submission script for each job.
