Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Conformational Sampling & Dynamics Jobs
========================================

ChemSmart provides powerful tools for conformational sampling and molecular dynamics calculations using Gaussian. This section covers CREST conformational searches and trajectory analysis workflows.

CREST/CREST-link jobs
---------------------

CREST (Conformer-Rotamer Ensemble Sampling Tool) is used for systematic conformational searches to find low-energy conformers.

Run "chemsmart sub [GENERAL_OPTIONS] gaussian [GAUSSIAN_OPTIONS] crest [OPTIONS]" to use CREST.

CREST-Specific Options
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

CREST-link jobs allow you to run additional Gaussian calculations on the generated conformers.

.. list-table:: CREST-Link Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-j, --jobtype``
     - string
     - Gaussian job type. Options: ["opt", "ts", "modred", "scan", "sp", "grouper"]

Basic Usage
^^^^^^^^^^^

* **Basic CREST job**:

    .. code-block:: console

        chemsmart sub gaussian -p project_name -f input.xyz crest

* **CREST with specific number of conformers**:

    .. code-block:: console

        chemsmart sub gaussian -p conformers -f molecule.xyz crest -N 10

* **CREST-link with optimization**:

    .. code-block:: console

        chemsmart sub gaussian -p crest_opt -f molecule.xyz crest -j opt

* **CREST-link with single point calculations**:

    .. code-block:: console

        chemsmart sub gaussian -p crest_sp -f molecule.xyz crest -j sp


Trajectory Analysis
-------------------

Trajectory analysis allows you to process molecular dynamics trajectories and extract specific structures for further analysis.

Trajectory-Specific Options
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

Basic Usage
^^^^^^^^^^^

* **Basic trajectory analysis**:

    .. code-block:: console

        chemsmart sub gaussian -p trajectory -f trajectory.xyz traj

* **Trajectory analysis with specific number of structures**:

    .. code-block:: console

        chemsmart sub gaussian -p traj_analysis -f md_output.xyz traj -N 50

Additional grouper option for crest/traj jobs
-------------------

Run "chemsmart sub gaussian [GAUSSIAN OPTIONS] crest/traj -g <> [OPTIONS]" to use

Specific Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Grouper Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-g, --grouping-strategy``
     - string
     - Grouping strategy to use for grouping. Options: "rmsd", "rmsd_simple", "tanimoto", "formula", "isomorphism", "connectivity" (default = "rmsd")
   * - ``-i, --ignore-hydrogens``
     - bool
     - Ignore H atoms in the grouping (Default = False)
   * - ``-p, --num-procs``
     - int
     - Number of processors to use for grouper (default=4)
   * - ``-x, --proportion-structures-to-use``
     - float
     - Proportion of structures from the end of trajectory to use. Values ranges from 0.0 < x <=1.0. Defaults to 0.1 (last 10% of structures)

Basic Usage
^^^^^^^^^^^
