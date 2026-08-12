Ubuntu x86 CPU Server and Private Container
############################################

This page defines the first portable ChemSmart execution environment for an
x86-64 Ubuntu CPU server. The preferred deployment is the private, pinned CPU
research image. A source installation remains available for development, but
it is not the reference benchmark environment.

Image scope
===========

The ``linux/amd64`` image contains ChemSmart 3.1.4, Python 3.10,
PySCF 2.14.0, xTB 6.7.1, NumPy 1.26.4, h5py 3.16.0,
geomeTRIC 1.1.1, pyscf-dispersion 1.5.0, PyMOL open source 3.1.0,
Open Babel, FFmpeg, ASE, and the remaining locked ChemSmart dependencies.

GPU4PySCF, CUDA, Gaussian, and ORCA are not installed. A PySCF GPU request is
rejected rather than silently executed on the CPU. Gaussian and ORCA require a
future separately qualified deployment.

The candidate workflow builds the exact selected ``main`` commit and records
that revision in the immutable image metadata. The private package is:

.. code-block:: text

   ghcr.io/hongjiseung-rok/chemsmart

``cpu-<12-character-main-commit>`` is an immutable candidate. After human
scientific review, the same digest may be promoted without rebuilding to
``3.1.4-cpu`` and ``cpu-main``.

Host requirements
=================

- Ubuntu 22.04 x86-64 with a current Docker Engine;
- a non-root server account permitted to run Docker;
- local writable workspace and scratch storage;
- sufficient CPU and RAM for the selected method; and
- authenticated read access to the private GHCR package.

The image favors portability over AVX2, AVX-512, or ``-march=native`` tuning.
VM provisioning and scheduler integration are later deployment stages.

Authenticate and pull
=====================

Use a GitHub token with ``read:packages``. Pass it through stdin and clear the
shell variable after login:

.. code-block:: bash

   export GHCR_READ_TOKEN='REDACTED'
   printf '%s' "$GHCR_READ_TOKEN" | \
     docker login ghcr.io --username Hongjiseung-ROK --password-stdin
   unset GHCR_READ_TOKEN

   docker pull ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu

The package owner must verify that the GHCR package is private, remove
inherited public-repository access if GitHub exposes it, and grant only the
intended package access.

Runtime contract
================

The process runs as UID/GID 1000, uses ``/workspace`` for calculations and
``/scratch`` for temporary engine data, and writes caches under ``/tmp``.
Make those locations writable and keep the remaining filesystem read-only:

.. code-block:: bash

   mkdir -p "$PWD/workspace" "$PWD/scratch"
   chmod u+rwx "$PWD/workspace" "$PWD/scratch"

   docker run --rm -it \
     --read-only \
     --ulimit stack=-1:-1 \
     --cpus 8 \
     --memory 32g \
     --tmpfs /tmp:rw,nosuid,nodev,size=4g \
     --mount "type=bind,src=$PWD/workspace,dst=/workspace" \
     --mount "type=bind,src=$PWD/scratch,dst=/scratch" \
     ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu \
     chemsmart --help

xTB requires the unlimited stack setting. The container entrypoint raises the
soft limit when allowed, but Docker's ``--ulimit stack=-1:-1`` establishes the
host-side hard limit.

Mounted directories must permit UID 1000 to write. On a shared server, use
ACLs or administrator-managed ownership rather than world-writable modes.

Configuration
=============

The image's ``~/.chemsmart/server/local.yaml`` declares local execution, two
CPU cores, 8 GB RAM, zero GPUs, PySCF and xTB from ``/opt/conda/bin``, and
``/scratch``. Resource flags on ``chemsmart run`` override these conservative
defaults.

Mount a site-specific server file if required:

.. code-block:: bash

   --mount type=bind,src=/secure/config/local.yaml,\
   dst=/home/chemsmart/.chemsmart/server/local.yaml,readonly

Do not store provider credentials, private registry tokens, or host secrets in
an image layer, YAML committed to Git, OCI label, or Docker build argument.

PySCF and xTB execution
=======================

Put project YAML and coordinates in the mounted workspace. The following
commands use the same canonical path available to the agent.

PySCF CPU single point:

.. code-block:: bash

   chemsmart run --server local --num-cores 4 --num-gpus 0 \
     --mem-gb 16 --no-scratch pyscf \
     --project /workspace/pyscf.yaml \
     --filename /workspace/water.xyz sp

xTB optimization:

.. code-block:: bash

   chemsmart run --server local --num-cores 4 --num-gpus 0 \
     --mem-gb 16 xtb \
     --project /workspace/xtb.yaml \
     --filename /workspace/water.xyz opt

For Hessian or frequency work, feed the validated optimized PySCF HDF5 or xTB
``xtbopt.xyz`` to the downstream command. Do not silently reuse the initial
geometry.

Live agent observation
======================

The image contains an Alibaba Token Plan profile selecting
``deepseek-v4-flash-0731`` with maximum reasoning effort, thinking
continuation, and no fallback. Supply the key only as a permission-restricted,
read-only runtime file:

.. code-block:: bash

   chmod 600 /secure/provider.env
   chemsmart agent plan \
     --provider alibaba-token-plan \
     --task-file /workspace/task.md \
     --secret-file /run/secrets/provider.env \
     --workspace /workspace

The Docker invocation must bind ``/secure/provider.env`` to
``/run/secrets/provider.env`` read-only. ``agent plan`` may create project YAML
and safe previews, but it cannot execute PySCF or xTB.

Qualification and scientific claims
===================================

The manual candidate workflow builds natively on GitHub's Ubuntu 22.04
x86-64 runner. It checks the read-only/non-root runtime and project previews,
then performs real neutral-singlet water calculations:

- PySCF B3LYP/def2-SVP SP, a separate D3BJ SP, geomeTRIC OPT, and HESS on the
  optimized HDF5 geometry; and
- GFN2-xTB SP, OPT, and HESS on the optimized ``xtbopt.xyz`` geometry.

Acceptance requires exact engine versions, normal termination, charge 0,
multiplicity 1, finite energy, optimization convergence, structured HDF5,
Hessian symmetry, parsed water frequencies, and explicit optimized-geometry
handoff. These checks qualify the installed engines; they do not grade agent
intelligence or reproduce a large scientific benchmark.

The workflow also records one preview-only model answer and visible tool
trajectory. A computational chemist reads it directly and accepts any
scientifically valid ChemSmart YAML/CLI decomposition. No answer-key, exact
DAG comparator, token threshold, or tool-count score decides acceptance.

Only the public calculation summary and visible agent observation are uploaded.
Raw workspaces, HDF5 artifacts, verbose logs, event stores, provider-private
reasoning, and credentials remain private.

Source installation for development
===================================

For source-level development rather than reference benchmark execution:

.. code-block:: bash

   git clone https://github.com/Hongjiseung-ROK/chemsmart.git
   cd chemsmart
   make env
   conda activate chemsmart
   make install-dev
   make configure

A separate PySCF compute environment should use the same direct pins:

.. code-block:: bash

   python -m pip install \
     "pyscf==2.14.0" "numpy==1.26.4" "h5py==3.16.0" \
     "geometric==1.1.1" "pyscf-dispersion==1.5.0"

Point ``PYSCF.EXEFOLDER`` at that environment's ``bin`` directory. This
source path is not byte-equivalent to the pinned container unless all
transitive dependencies and server configuration are also reproduced.

Future server stage
===================

After the image digest is promoted and pulled successfully:

1. provision the Ubuntu VM and durable workspace/scratch volumes;
2. configure private GHCR read credentials outside the image;
3. qualify the image on the actual VM with its CPU, memory, filesystem, and
   network limits;
4. run a small real benchmark through ChemSmart YAML/CLI; and
5. add Gaussian, ORCA, MPI, or a scheduler only as separately qualified
   execution layers.

See :doc:`agent-workflows`, :doc:`pyscf-cli-options`,
:doc:`xtb-cli-options`, and :doc:`configuration-server-settings`.
