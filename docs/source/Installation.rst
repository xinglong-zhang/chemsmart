Installation
====================

Thank you for downloading and using chemsmart. Please follow the steps for installation and testing. If you encounter any problems, please feel free to contact us.

Installation and configuration
-------------------------------

Users are recommended to install conda environments to mange the packages required by this software toolkit. Either Anaconda3 or Miniconda3 may be installed. See more information on conda installation process at https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html.

*   Once conda has been installed successfully, one can clone this package to a local directory via

    .. code-block:: console

        git clone https://github.com/xinglong-zhang/chemsmart.git

*   Cd into chemsmart folder. For linux and MacOS systems which support ``make``, users can run

    .. code-block:: console

        make env

    to create a running environment.

    By default, this will create a conda environment named ``chemsmart``, which installs all the required python packages for this toolkit.

    If conda is not installed, one can run

    .. code-block:: console

        make env USE_CONDA=false

    or

    .. code-block:: console

        make virtualenv

    to install using virtualenv. It is however recommanded that ``conda`` be used.

.. tip::

    Help options are available by typing ``make help``.

*   After the virtual conda environment is created and activated via：

    .. code-block:: console

        conda activate chemsmart

    one can run

    .. code-block:: console

        make install

    which installs the packages and dependencies required for ``chemsmart`` package.

*   Next, one can run

    .. code-block:: console

        make configure

    to sets up the user-specific directory ``~/.chemsmart`` automatically. You will be prompt to enter the paths to g16 and ORCA software, which will then be added automatically. The correct ``conda`` path for the user will also be updated.

    **The configuration also adds the environment variables for chemsmart to the user ``~/.bashrc`` file.**

    The ``~/.chemsmart/usersettings.yaml`` file contains informations such as project number or account number that are required in a typical submission script that specifies the account for use at some HPC servers. It can also contain options specifying user's email to inform user of the job start and job end once a job is submitted. If more features are needed, please submit a request via `Issues`. A typical `~/.chemsmart/usersettings.yaml` file looks like this:

    .. code-block:: console

        PROJECT: 1234567  # alias ACCOUNT FOR SLURM
        EMAIL: abc@gmail.com

    The ``~/.chemsmart/server/`` directory contains files related to server setup for a particular HPC cluster that the user is using. For example, we can specify a SLURM based server setting as ``~/.chemsmart/server/shared.yaml`` with the following information:

    .. code-block:: console

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

    This file can be customized by user for different submission systems. This file contains the server configuration information that is needed for chemsmart to automatically write the submission script for each job.

    The ``~/.chemsmart/gaussian/`` directory contains files related to gaussian project settings, which contain DFT functional and basis set etc, that is required to write the input file for running a gaussian job. For example, we can specify a test project settings in ``~/.chemsmart/gaussian/test.yaml`` with the following information:

    .. code-block:: console

        gas:
          functional: m062x  # quotes required for string with spaces
          basis: def2svp
          solvent_model: smd
          solvent_id: dichloroethane
        solv:
          functional: m062x
          basis: def2tzvp
          freq: False
          solvent_model: smd
          solvent_id: dichloroethane
        td:
          functional: cam-b3lyp
          basis: genecp
          heavy_elements: ['I']
          heavy_elements_basis: def2-SVPD
          light_elements_basis: def2SVP
          freq: False
          ##solvent_model: smd
          ##solvent_id: DiethylEther

    By default, the ``gas`` phase settings are used for all jobs such as geometry optimization, transition state search etc, and the ``solv`` settings are used for single point calculations; the ``td`` settings are used to run TD-DFT calculations. One can specify additional project settings in ``~/.chemsmart/gaussian/`` in a similar way to adapt to each project that one wishes to run. If setting

    .. code-block:: console

        gas: Null

    Then all jobs will use settings specified in ``solv``, i.e., all calculations will be run in implicit solvation model.


    The ``~/.chemsmart/orca/`` directory contains files related to ORCA project settings, which contain DFT functional and basis set etc, that is required to write the input file for running an ORCA job. For example, we can specify a test project settings in ``~/.chemsmart/orca/test.yaml`` with the following information:

    .. code-block:: console

        gas:
          functional: M062X
          basis: def2-SVP
        solv:
          ab_initio: DLPNO-CCSD(T)
          functional: Null
          basis: Extrapolate(2/3,cc)
          aux_basis: AutoAux
          defgrid: DEFGRID3
          freq: False
          scf_tol: TightSCF
          scf_algorithm: KDIIS
          scf_maxiter: 500
          mdci_cutoff: Normal
          mdci_density: None
          dipole: False
          solvent_model: SMD
          solvent_id: "toluene"

    This will run jobs in the gas phase (geometry and TS opt etc) using M062X/def2-SVP method and run single point with solvent correction using DLPNO-CCSD(T)/CBS with cc-pVDZ/cc-pVTZ extrapolation in SMD(toluene), for example. Again, users can customize different settings in different ``~/.chemsmart/orca/*project_settings*.yaml`` files to adapt to different project requirements.


.. warning::

    ``make configure`` would set up ``~/.chemsmart`` mostly correctly, a user should check the contents in ``~/.chemsmart`` to make sure that these match the **server configurations** on which chemsmart is to be used (e.g., modules, scratch directories etc). Depending on the server queue system you are using (e.g., SLURM or TORQUE), one may copy e.g., ``~/.chemsmart/server/SLURM.yaml`` to your own customised server ``~/.chemsmart/server/custom.yaml`` and modify it accordingly, such that the submission becomes ``chemsmart sub -s custom <other commands>``.

.. tip::

    One also need to set up scratch directories where scratch jobs may be run (for Gaussian and ORCA jobs, by default, these are run in scratch folder), one may do ``ls -s /path/to/scratch/ ~/scratch``.

.. note::
    A user can modify the contents in ``~/.chemsmart`` files freely without affecting or needing to know the ``chemsmart`` source code.

*   The ``make configure`` will also add the required paths to the user ``~/.bashrc`` file. User may need to do

    .. code-block:: console

        source ~/.bashrc

    to effect the changes.


*   Once ``make configure`` is done, one can optionally run

    .. code-block:: console

        make fmt

    and

    .. code-block:: console

        make lint

    to format and lint the codes (this should have been handled by the developers). Also optionally, one can run

    .. code-block:: console

        make test

    to make sure that all tests in chemsmart pass.


*   Finally one can clean up by running

    .. code-block:: console

        make clean

Testing Installations
-------------------------------

Installations is deemed successfully if the commands ``make install`` and ``make configure`` do not return any errors. Installation will also create a ``~/.chemsmart`` containing the required files. In addition, the paths for chemsmart packages should be correctly added to the user ``~/.bashrc`` file. Finally, one should be able to run

.. code-block:: console

    chemsmart --help

to get the options for running chemsmart package.