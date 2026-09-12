##########################
 GPU Acceleration (PySCF)
##########################

This page covers running PySCF jobs on the GPU4PySCF engine and configuring the interpreter that owns it.

GPU4PySCF is an execution engine of the PySCF program, not a separate program. The same ``sp``, ``opt``, and ``hess``
subcommands and the same project YAML apply; only the engine changes.

**********************
 Selecting the Engine
**********************

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--gpu/--no-gpu``
      -  bool
      -  Run on GPU4PySCF or on CPU

   -  -  ``-g, --num-gpus``
      -  int
      -  Number of GPUs per node for the scheduler

``--gpu`` chooses the engine PySCF runs on. ``-g`` tells the scheduler how many devices to reserve and belongs to
``chemsmart run`` and ``chemsmart sub``, so it comes before the program name. With neither ``--gpu`` nor ``--no-gpu``,
the engine follows ``SERVER.NUM_GPUS`` from the server YAML.

Basic Usage
===========

Single point on the GPU engine:

.. code:: bash

   chemsmart sub -s GPU_SERVER -g 1 pyscf -p test -f molecule.xyz --gpu sp

Optimization on two devices:

.. code:: bash

   chemsmart sub -s GPU_SERVER -g 2 pyscf -p test -f molecule.xyz --gpu opt

Force CPU on a machine that has GPUs available:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz --no-gpu opt

GPU and CPU results are not bit-identical, so the engine that was used is recorded in the results file. ``--no-gpu``
also clears ``CUDA_VISIBLE_DEVICES`` for the calculation, which keeps a CPU run reproducible on a GPU node.

***********************************
 Configuring the PySCF Interpreter
***********************************

PySCF is a Python library, not a binary, so what matters is the interpreter that runs the calculation, not the one
running chemsmart. For GPU work CuPy and a CUDA toolchain usually live in their own environment.

Point chemsmart at the ``bin`` directory whose ``python`` owns PySCF:

.. code:: bash

   chemsmart config pyscf -f ~/miniconda3/envs/pyscf-gpu/bin

The same setting can be written by hand as a ``PYSCF`` block in the server YAML:

.. code:: yaml

   PYSCF:
     EXEFOLDER: ~/miniconda3/envs/pyscf-gpu/bin
     LOCAL_RUN: True
     SCRATCH: False

Omit ``EXEFOLDER`` when PySCF shares the environment chemsmart runs in, and the running interpreter is used. Keeping a
separate CUDA environment is the usual arrangement, since it keeps CuPy away from the NumPy that RDKit and PyMOL pin.

**************************
 When No GPU Is Available
**************************

.. warning::

   A ``--gpu`` request is never downgraded to CPU. If the GPU4PySCF stack or the device is missing, the job fails: the
   error is recorded in the results file with its stage and traceback, ``normal_termination`` is false, and the process
   exits non-zero.

Not every method runs on the GPU engine even with a device present. Double-hybrid functionals, Laplacian-dependent
meta-GGAs, and basis sets with high angular momentum fall outside what GPU4PySCF covers, and a solvent model it does not
implement will fail there while working on CPU.

Use ``--fake`` to render a GPU job on a machine with no device:

.. code:: bash

   chemsmart run --fake pyscf -p test -f molecule.xyz --gpu sp

This writes the driver script without importing PySCF, which is a quick check that a project, basis set, and functional
produce the script you expected before asking a scheduler for GPU time. See :doc:`pyscf-results`.
