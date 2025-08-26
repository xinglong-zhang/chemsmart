Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Conformational Sampling & Dynamics Jobs
========================================

ChemSmart provides powerful tools for conformational sampling and molecular dynamics calculations using Gaussian. This section covers CREST conformational searches and trajectory analysis workflows.

CREST/CREST-link jobs
---------------------

CREST (Conformer-Rotamer Ensemble Sampling Tool) is used for systematic conformational searches to find low-energy conformers.

.. code-block:: console

    chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] crest [SUBCMD_OPTIONS]

CREST-link jobs allow you to run additional Gaussian calculations on the generated conformers.

.. code-block:: console

    chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] crest -j <jobtype> [SUBCMD_OPTIONS]

CREST-Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: CREST Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-N, --num-confs-to-run``
     - int
     - Number of conformers to optimize
   * - ``-j, --jobtype``
     - string
     - Gaussian job type. Options: ["opt", "ts", "modred", "scan", "sp", "grouper"]

Basic Usage
^^^^^^^^^^^

**Basic CREST job**:

    .. code-block:: console

        chemsmart sub gaussian -p project_name -f input.xyz crest

**CREST with specific number of conformers**:

    .. code-block:: console

        chemsmart sub gaussian -p conformers -f molecule.xyz crest -N 10

**CREST-link jobs**:

*   To run opt or modred or ts conformers from CREST run output, do:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test-f <input_file> -c 1 -m 1 crest -j opt

    or

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f <input_file> -c 0 -m 2 crest -j modred -c [1,4]

    or

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f <input_file> -c 0 -m 1 crest -j ts

    respectively

.. note::

     Typically, the ``<input_file>`` is a list of all conformers obtained by CREST program and named ``crest_conformers.xyz``.


Examples
^^^^^^^^^^^

!need to add example here!


Trajectory/Traj-link Analysis
-------------------

Trajectory analysis allows you to process molecular dynamics trajectories and extract specific structures for further analysis.

.. code-block:: console

    chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] traj [SUBCMD_OPTIONS]

.. note::

    Charge and multiplicity need to be specified, as these cannot be obtained from the supplied .traj file.

Trajectory-Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Trajectory Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-N, --num-structures-to-run``
     - int
     - Number of structures from the list of unique structures to run the job on
   * - ``-x, --proportion-structures-to-use``
     - float
     - Proportion of structures from the end of trajectory to use. Values ranges from 0.0 < x <= 1.0. Defaults to 0.1 (last 10% of structures)

Basic Usage
^^^^^^^^^^^

**Basic trajectory analysis**:

    .. code-block:: console

        chemsmart sub gaussian -p trajectory -f trajectory.xyz -c 0 -m 1 traj

**Trajectory analysis with specific number of structures**:

    .. code-block:: console

        chemsmart sub gaussian -p traj_analysis -f md_output.xyz -c 0 -m 1 traj -N 50

**Trajectory analysis with specific proportion of structures and sequential grouping**:

*   to consider the last 20% of the structures in md.traj trajectory file, then uses sequential grouper to group those structures into unique structures and run the 10 lowest energy structures from the list of unique structures found by the grouper:

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f imd.traj traj -x 0.2 -N 10 -g rmsd


Additional grouper option for crest/traj jobs
-------------------
Process the results of the crest/traj task further using multiple molecular similarity-based grouping strategies.

.. code-block:: console

        chemsmart sub gaussian [GAUSSIAN OPTIONS] crest/traj -g <> [SUBCMD_OPTIONS]

Grouper-Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Grouper Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-g, --grouping-strategy``
     - string
     - Grouping strategy to use for grouping. Options: "rmsd", "tanimoto", "formula", "isomorphism", "connectivity" (default = "rmsd")
   * - ``-i, --ignore-hydrogens``
     - bool
     - Ignore H atoms in the grouping (Default = False)
   * - ``-t, --threshold``
     - float
     - Threshold value for grouping (Default = 0.5 for rmsd, 0.9 for tanimoto)
   * - ``-p, --num-procs``
     - int
     - Number of processors to use for grouper (Default=1)


chemsmart sub -s small gaussian -p test -f 1.traj -c 1 -m 1 traj -j opt -x 0.2 -n 10 -g rmsd


Examples
^^^^^^^^^^^

chemsmart run gaussian -p test -f crest_conformers.xyz -l grouped -c 0 -m 1 crest -j opt -g rmsd -t 0.2 -p 4

!need to add example here!
