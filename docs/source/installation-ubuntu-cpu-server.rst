################################
 Ubuntu CPU Target-Host Example
################################

This page is one target-host example for running CHEMSMART on an Ubuntu CPU server. It separates the CHEMSMART
controller from external chemistry programs and introduces real execution in stages. It is not a universal reference:
verify architecture, scheduler, resources and program compatibility on the actual machine.

CHEMSMART does not distribute Gaussian, ORCA, xTB, NCIPLOT, PySCF, MPI, or provider credentials. Install and license
each required backend separately.

****************************
 Example host prerequisites
****************************

-  an actively maintained Ubuntu release on a supported architecture;
-  resources sized for the selected molecule and method;
-  local SSD scratch space sized for the selected methods;
-  outbound HTTPS for package installation and authorised model providers; and
-  a non-root account for all calculations.

Large post-HF, numerical Hessian, IRC, NEB, or conformer workflows may require substantially more memory, time, and
scratch space. Resource requirements are part of the scientific plan, not universal defaults.

*****************
 System packages
*****************

.. code:: bash

   sudo apt-get update
   sudo apt-get install -y build-essential ca-certificates curl git \
     libegl1 libgl1 libsm6 libxext6 libxrender1

Install the Conda distribution appropriate for the host architecture if the server does not already provide one.

**********************************
 Install the CHEMSMART controller
**********************************

.. code:: bash

   git clone https://github.com/Hongjiseung-ROK/chemsmart.git
   cd chemsmart
   conda env create -f environment.yml
   conda activate chemsmart
   python -m pip install .
   chemsmart config

The package requires Python 3.10. ``chemsmart config`` creates local configuration under ``~/.chemsmart``. Never commit
that directory because it contains host paths and scheduler settings.

Verify the controller before configuring engines:

.. code:: bash

   chemsmart --version
   chemsmart --help
   chemsmart run --help
   chemsmart agent --help

***********************************
 Configure a local CPU server file
***********************************

Start from the generated ``~/.chemsmart/server/local.yaml`` and replace every example path with a path observed on the
server. A minimal shape is:

.. code:: yaml

   SERVER:
     SCHEDULER: null
     QUEUE_NAME: null
     NUM_HOURS: 1
     MEM_GB: 32
     NUM_CORES: 4
     NUM_GPUS: 0
     NUM_THREADS: 4
     SUBMIT_COMMAND: null
     SCRATCH_DIR: /absolute/user-owned/scratch/chemsmart
     USE_HOSTS: false

   ORCA:
     EXEFOLDER: /opt/orca
     LOCAL_RUN: true
     SCRATCH: true
     ENVARS: |
       export PATH=/opt/openmpi/bin:$PATH
       export LD_LIBRARY_PATH=/opt/openmpi/lib:$LD_LIBRARY_PATH

   XTB:
     EXEFOLDER: /opt/xtb/bin
     LOCAL_RUN: true
     SCRATCH: true

   PYSCF:
     EXEFOLDER: /opt/conda/envs/chemsmart-pyscf/bin
     LOCAL_RUN: true
     SCRATCH: false

Use an absolute user-owned path in ``SCRATCH_DIR``. Create the directory and verify write permission before execution.

****************************
 ORCA and MPI compatibility
****************************

Install ORCA from its authorised distribution. Parallel ORCA builds require a compatible MPI implementation; a binary
linked to Open MPI must not be launched with MPICH tools. Check the binary and launcher on the target server:

.. code:: bash

   /opt/orca/orca --version
   /opt/openmpi/bin/mpirun --version
   ldd /opt/orca/orca | grep -i mpi

First qualify ORCA with one core. Increase the core count only after serial normal termination, parsing, and a small
parallel calculation are all green. A zero process exit status is not enough; CHEMSMART also checks ORCA's normal
termination and job-specific convergence.

*****************
 xTB CPU backend
*****************

Install the maintained xTB release in a dedicated program environment or a server-managed location. Confirm the exact
executable that server YAML selects:

.. code:: bash

   /opt/xtb/bin/xtb --version

The maintained CHEMSMART xTB surface is CPU ``sp``, ``opt``, and ``hess``. Unknown native features are not silently
added to generated commands.

*******************************
 PySCF CPU compute environment
*******************************

PySCF is a Python library backend. Keep its compute interpreter separate from the controller when binary or scientific
dependencies differ:

.. code:: bash

   conda create -n chemsmart-pyscf python=3.10 pip -y
   conda activate chemsmart-pyscf
   python -m pip install \
     "pyscf==2.14.0" "h5py>=3.10,<4" "numpy<2" \
     "geometric==1.1.1" "pyscf-dispersion==1.5.0"
   python -c "import h5py, pyscf; print(pyscf.__version__, h5py.__version__)"

Point ``PYSCF.EXEFOLDER`` at this environment's ``bin`` directory. CHEMSMART generates a standalone Python calculation
and reads its structured HDF5 result; the compute environment does not need to import the CHEMSMART source tree.

GPU4PySCF requires a separately qualified NVIDIA/CUDA stack. Do not select it as an automatic fallback on a CPU-only
host.

*****************************
 Staged server qualification
*****************************

Run stages in order and stop on the first scientific or environment failure.

#. **Controller and schema**

   .. code:: bash

      chemsmart --version
      chemsmart run orca --help
      chemsmart run pyscf --help
      chemsmart run xtb --help

#. **Program-free previews**

   .. code:: bash

      chemsmart run --fake --no-scratch xtb \
        -p test -f examples/xtb/water.xyz sp
      chemsmart run --fake --no-scratch orca \
        -p test -f examples/xtb/water.xyz sp
      chemsmart run --fake --no-scratch pyscf \
        -p test -f examples/xtb/water.xyz sp

#. **Small real single points**

   Execute one backend at a time with one core and a small closed-shell molecule. Inspect generated input, selected
   program version, normal termination, parsed energy, charge, multiplicity, atom order and units.

#. **Geometry and frequencies**

   Qualify OPT and HESS/FREQ separately, then verify the downstream frequency calculation consumes the validated
   optimized geometry.

#. **Advanced workflows**

   Introduce TS, IRC, NEB, numerical derivatives, post-HF composites and thermochemistry only after the underlying
   methods and resources are known to work on this host.

#. **Agent workflow**

   Run ``chemsmart agent plan`` first. Permit real engine execution only after reviewing and approving the exact
   generated YAML, DAG, compiled command, molecular state, environment, and resources.

************************
 Provider configuration
************************

Set up provider profiles and credentials with the guided command; it writes ``~/.chemsmart/agent/agent.yaml`` and stores
the key in the managed store at ``~/.chemsmart/agent/keys.env`` (0600), which is parsed as data and never sourced or
committed.

.. code:: bash

   chemsmart config agent
   chemsmart agent plan \
     --provider PROFILE \
     --task-file task.md \
     --workspace ./agent-workspace

See :doc:`agent-workflows` for the scientific and execution boundary.

*****************
 Scheduler hosts
*****************

For SLURM, PBS/Torque, or LSF, create a separate server YAML and use ``chemsmart sub``. Generate a script with
``--test`` before any submission. Scheduler modules, accounts, queues, wall time and scratch policy are host facts and
must be confirmed with the administrator.

See :doc:`installation-hpc-cluster` and :doc:`configuration-server-settings` for the submission model.
