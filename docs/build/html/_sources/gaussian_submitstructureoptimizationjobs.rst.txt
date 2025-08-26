Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Structure Optimization Jobs
===================================

ChemSmart provides comprehensive structure optimization capabilities using Gaussian. This section covers basic geometry optimization, constrained optimization with frozen atoms, modified redundant coordinate optimization, and CREST optimization workflows.


Geometry Optimization
----------------------------

Run basic geometry optimization to find the minimum energy structure using the ``opt`` command.

.. code-block:: console

    chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] opt

Run constrained optimization with specific atoms frozen in place using the ``opt`` command with frozen atom options (such as https://www.researchgate.net/post/Freezing-atoms-in-gaussian-how-to-do-it).

.. code-block:: console

    chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] opt -f <atom_indices_to_freeze>

Opt-Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Frozen Atom Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-f, --freeze-atoms``
     - string
     - Indices of atoms to freeze for constrained optimization (1-based indexing)

Basic Usage
^^^^^^^^^^^

**Basic geometry optimization**:

    .. code-block:: console

        chemsmart sub gaussian -p project_name -f input.gjf -c 0 -m 1 opt

**Constrained optimization with frozen atoms**:

*   to submit the geometry optimization job with atoms numbered 1 to 10 frozen, one can do

    .. code-block:: console

        chemsmart sub -s shared gaussian -p frozen_opt -f input.com opt -f 1-10

.. Warning::

    1-indexed numbers are used, instead of 0-indexed numbers in Python language, since most visualization softwares for moleculare are 1-indexed.

Examples
^^^^^^^^

!need to add!


CREST Optimization Jobs
-----------------------

Run CREST-based optimization on multiple conformers using the ``crestopt`` command.

.. code-block:: console

    chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] crestopt [SUBCMD_OPTIONS]

CRESTOPT-Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: CREST Optimization Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-n, --num-confs-to-opt``
     - int
     - Number of conformers to optimize from the CREST ensemble


Basic Usage
^^^^^^^^^^^^^^

**CREST Optimization for all conformers**

    .. code-block:: console

        chemsmart sub gaussian -p crest_optimization -f conformers.xyz crestopt

**with specific number of conformers**

    .. code-block:: console

        chemsmart sub gaussian -p crest_opt -f molecule.xyz crestopt -n 10

.. note::

    If the job terminates before ``<n_conformers_to_opt>`` are all optimized, perhaps due to walltime limit, resubmitting the job will continue crest opt job until all ``<n_conformers_to_opt>`` are optimized. Charge and multiplicity need to be specified.



