Tutorials
=========

Welcome to the tutorials! We’re thrilled to have you here. Please go through the code examples, and don’t hesitate to contact our team if you have questions or feedback.

How to Set up the Server Configuration on the Cluster
-------------------------------

Once the user has completed the installation, the server configuration on the cluster can be customized by modifying corresponding files.

.. note::
    A user can modify the contents in ``~/.chemsmart`` files freely without affecting or needing to know the ``chemsmart`` source code.

*   The ``~/.chemsmart/usersettings.yaml`` file contains informations such as project number or account number that are required in a typical submission script that specifies the account for use at some HPC servers. It can also contain options specifying user's email to inform user of the job start and job end once a job is submitted. If more features are needed, please submit a request via `Issues`. A typical ``~/.chemsmart/usersettings.yaml`` file looks like this:

        .. code-block:: console

            PROJECT: 1234567  # alias ACCOUNT FOR SLURM
            EMAIL: abc@gmail.com

*   The ``~/.chemsmart/server/`` directory contains files related to server setup for a particular HPC cluster that the user is using. For example, we can specify a SLURM based server setting as ``~/.chemsmart/server/shared.yaml`` with the following information:

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

How to Set up the Project Settings
-------------------------------

Chemsmart allows user to execute projects either locally or on a cluster. Users can refer to the following examples to configure projects as required.

*    The ``~/.chemsmart/gaussian/`` directory contains files related to gaussian project settings, which contain DFT functional and basis set etc, that is required to write the input file for running a gaussian job. For example, we can specify a test project settings in ``~/.chemsmart/gaussian/test.yaml`` with the following information:

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

*   The ``~/.chemsmart/orca/`` directory contains files related to ORCA project settings, which contain DFT functional and basis set etc, that is required to write the input file for running an ORCA job. For example, we can specify a test project settings in ``~/.chemsmart/orca/test.yaml`` with the following information:

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

*   One also need to set up scratch directories where scratch jobs may be run (for Gaussian and ORCA jobs, by default, these are run in scratch folder), one may do ``ls -s /path/to/scratch/ ~/scratch``.



How to Submit a Geometry Optimization Job Using Gaussian
-------------------------------
With setup completed, one is able to run different Gaussian jobs via command-line interface (CLI).

Basic Geometry Optimization
^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit (and run) a geometry optimization job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> opt

    where ``<server_name>`` is the one of the servers specified in ``~/.chemsmart/server/*.yaml`` files, without ``.yaml`` extension; ``<project>`` is one of the project settings specified in ``~/.chemsmart/gaussian/*.yaml`` files, without ``.yaml`` extension; and ``<input_file>`` is an input file the user wishes to run job on. Note that this input file can be any format, such as ``.xyz``, Gaussian ``.com``, ``.gjf`` or ``.log`` file or ORCA ``.inp`` or ``.out`` file.

Geometry Optimization Job with Frozen Atoms
^^^^^^^^^^^^^^^^^^^^^^^^

*   If one wants to submit geometry optimization with frozen atoms in the molecule (such as https://www.researchgate.net/post/Freezing-atoms-in-gaussian-how-to-do-it), one can do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> opt -f <indices_of_atoms_to_freeze>

    or example, to submit the geometry optimization job with atoms numbered 1 to 10 frozen, one can do

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f input.com opt -f 1-10

.. Note::

    1-indexed numbers are used, instead of 0-indexed numbers in Python language, since most visualization softwares for moleculare are 1-indexed.

Optimize a Fixed Number of Lowest Energy Conformers
^^^^^^^^^^^^^^^^^^^^^^^^

*   To optimize a fixed number of lowest energy conformers, ``n_conformers_to_opt``, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j opt -n <n_conformers_to_opt>

.. note::

    If the job terminates before ``<n_conformers_to_opt>`` are all optimized, perhaps due to walltime limit, resubmitting the job will continue crest opt job until all ``<n_conformers_to_opt>`` are optimized. Charge and multiplicity need to be specified, as these cannot be obtained from the supplied .xyz file. In fact, whenever .xyz file is used as input, the charge and multiplicity should be specified via ``-c <charge> -m <multiplicity`` via CLI.

Optimize Unique Structure from an MD Trajectory File
^^^^^^^^^^^^^^^^^^^^^^^^

*   To optimize unique structure from an MD trajectory file, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> saopt

.. note::

    This optimizes all the unique structures available in the MD trajectory ``<input_file>``. Typically, the ``<input_file>`` is a list of all structures on an MD trajectory obtained by ASE MD run and named ``md.traj``. **(TODO: this method is not properly implemented in chemsmart yet.)**

Optimize a Fixed Number of Lowest Energy Structures
^^^^^^^^^^^^^^^^^^^^^^^^

*   To optimize a fixed number of lowest energy structures, ``<num_structures_to_opt>``, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> saopt -n <n_conformers_to_opt>


    If the job terminates before ``<n_conformers_to_opt>`` are all optimized, perhaps due to walltime limit, resubmitting the job will continue crest opt job until all ``<n_conformers_to_opt>`` are optimized. Charge and multiplicity need to be specified, as these cannot be obtained from the supplied .traj file.

    Two grouper types for determining/clustering unique structures are available from CLI option ``-g``:

    1. Sequential grouper (default), selected by option value of seq, which sequentially checks for unique structures in a given list of md structures, and

    2. Self-consistent grouper, selected by option value of sc, which self-consistently checks for unique structures in a given list of md structures using the reverse Cuthill–McKee algorithm for structure clustering. By default, only the last 0.1 proportion of the structures of the md.traj file is considered. This can be changed via cli option ``-x <proportion_structures_to_use>``.

    For example, to consider the last 20% of the structures in md.traj trajectory file, then uses Sequential grouper to group those structures into unique structures and run the 10 lowest energy structures from the list of unique structures found by the grouper:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f imd.traj saopt -x 0.2 -n 10 -g seq

run opt or modred or ts conformers from CREST run output
^^^^^^^^^^^^^^^^^^^^^^^^

*   To run opt or modred or ts conformers from CREST run output, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j opt

    or

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j modred -c [1,4]

    or

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> -c <system_charge> -m <system_multiplicity> crest -j ts

    respectively

.. note::

    This optimizes all the conformers available in the ``<input_file>``. Typically, the ``<input_file>`` is a list of all conformers obtained by CREST program and named ``crest_conformers.xyz``.



How to Submit a Transition State and Reaction Pathway Analysis Job Using Gaussian
-------------------------------

Transition State Modredundant
^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit transition state modredundant job (frozen coordinates optimization), do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> modred -c <list_of_coords_to_constraint>

    For example, to submit a modredundant job with constraints on bond between atom 4 and atom 17 and on bond between atom 9 and atom 10, do:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f input.com modred -c [[4,17],[9,10]]

Transition State Search
^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit transition state search job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> ts

Intrinsic Reaction Coordinate (IRC)
^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit intrinsic reaction coordinate (IRC) job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> irc

Relaxed Potential Energy Surface (PES) Scan
^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit relaxed potential energy surface (PES) scan, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> scan -c <list_of_coords_to_constraint> -s <scan_step_size> -n <num_scan_steps>

*   For example, to submit the PES scan job with along bond between atom 4 and atom 17 for 10 steps with 0.1Å increment per step:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f input.com scan -c [[4,17]] -s 0.1 -n 10



How to Submit a Energy and Thermodynamic Analysis Job Using Gaussian
-------------------------------

Single Point Job
^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit single point job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> sp

    For single-point job that user wants to test which uses different solvent model and id from that specified in ``<project>``, one can do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> sp -sm <user_solvent_model> -si <user_solvent_id>

    to specify a different solvent model ``<user_solvent_model>`` and solvent ``<user_solvent_id>``.

Distortion-Interaction/Activation-Strain (DI-AS)
^^^^^^^^^^^^^^^^^^^^^^^^

*   To run distortion-interaction/activation-strain (DI-AS) job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <irc_output_file_for_dias> dias -i <indices_of_any_one_fragment> -n <number_of_every_n_step_along_irc_to_run>

    For example to run DI-AS job for fragment 1 with atoms numbered from 5-17 at every 10 steps along the irc.log file:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f irc.log dias -i 5-17 -n 10



How to Submit a Electronic Structure Analysis Job Using Gaussian
-------------------------------

RESP Charges Fitting
^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit RESP charges fitting job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> resp

.. note::

    This creates an input file with fix route for RESP job:
    ``HF/6-31+G(d) SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=17,6/50=1)``

Non-covalent Interaction Job
^^^^^^^^^^^^^^^^^^^^^^^^

*   To submit non-covalent interaction job, do:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <project> -f <input_file> nci



How to Submit Custom Job Submission Using Gaussian
-------------------------------

Run a Job with Pre-prepared Gaussian Input File Directly
^^^^^^^^^^^^^^^^^^^^^^^^

*   If a user wants to run a job with pre-prepared Gaussian input file directly, one can run the job directly using:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <any_project_settings> -f <input_file> com

Run a Custom Job
^^^^^^^^^^^^^^^^^^^^^^^^

*   Generally, if a user wants to run job that is currently not present in our package, one can run custom job using:

    .. code-block:: console

        chemsmart sub -s <server_name> gaussian -p <any_project_settings> -f <input_file> userjob -r <user_defined_gaussian_route> -a <appending_information_as_string_at_the_end_of_input_file_after_coordinates_specification>

    For example, to create an input file named ``user_defined_job.com`` with user-specified route ``mnr functional/basis solvent`` etc and ``B 1 2 F\nA 1 2 3 F`` at the end of the input file after the specification of coordinates, run

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -l user_defined_job userjob -r 'mnr functional/basis solvent etc' -a 'B 1 2 F\nA 1 2 3 F'



How to Use General CLI Options for All Jobs
-------------------------------

Specify Name of the File
^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can specify the name of the file to be created for the job, without file extension, they want to run by using the option ``-l``, e.g.:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -l custom_job_name opt

    will create input file named ``custom_job_name.com`` instead of the default ``test_opt.com``.

Append String to Input File Name
^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also simply append a string to the base name of the filename supplied, e.g.:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -a append_string ts

    will create input file named ``test_append_string.com`` instead of the default ``test_ts.com``.

Modify the Charge and Multiplicity
^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also modify the charge and multiplicity from the CLI, e.g.:

    Modify the charge in ``test.com`` to charge of +1 in the newly created input file ``test_charge.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -c 1 -a charge opt

    Modify the multiplicity in ``test.com`` to multiplicity of 3 in the newly created input file ``test_multiplicity.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -m 3 -a multiplicity opt

    Modify the charge to +1 and multiplicity to 2 in the newly created input ``file test_charge_multiplicity.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -c 1 -m 2 -l test_charge_multiplicity opt

.. tip::

    This can be useful when, e.g., using optimized structure of a neutral closed-shell (charge 0, multiplicity 1) system to run a charged radical ion (e.g., charge +1 and multiplicity 2 in radical cation).

Modify the Functional and Basis
^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also modify the functional and basis from the CLI to differ from those in project settings, e.g.:

    Modify the functional to ``b3lyp`` in the newly created input file ``test_functional.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -x b3lyp -a functional opt

    Modify the basis to ``6-31G*`` in the newly created input file ``test_basis.com`` via:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -b "6-31G*" -a basis opt

Specify Additional Optimization Options
^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also specify additional optimization options for opt=() in the route, for example,

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -o maxstep=8,maxsize=12 -a opt_options opt

    will create ``opt=(maxstep=8,maxsize=12)`` as part of the route in the newly created input file ``test_opt_options.com``.

Add in Additional Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

*   Users can also add in additional parameters used in the route, e.g.:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com --r nosymm -a route_params opt

    will add in ``nosymm`` as part of the route in the newly created input file ``test_route_params.com``.

Select the Particular Structure in file
^^^^^^^^^^^^^^^^^^^^^^^^
*   If one has more than one structure in the supplied file for input preparation, one can select the particular structure to perform job on by using the ``-i/--index`` option, e.g.:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f small.db -i 5 -c 0 -m 1 opt

    will take the 5th structure (1-indexed, as in chemsmart) from ase database file, ``small.db``, to create the input file for geometry optimization.



How to Submit Jobs Using ORCA
-------------------------------
Similar commands exists for ORCA job submssions.

*   One can run

        .. code-block:: console

            chemsmart sub orca --help

    to find out more.
