Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Electronic Structure Properties & Analyses Jobs
======================================================

ChemSmart provides comprehensive electronic structure analysis capabilities using Gaussian. This section covers single point calculations, excited state properties, bond analysis, and molecular interaction studies.

Single Point Jobs
-----------------

Run single point energy calculations on optimized geometries.

Basic Usage
^^^^^^^^^^^

* **Basic single point calculation**:

    .. code-block:: console

        chemsmart sub gaussian -p sp_calc -f optimized.log sp

* **Single point with solvent settings**:

    .. code-block:: console

        chemsmart sub gaussian -p sp_solvent -f molecule.log sp -sm smd -si water

DI-AS Jobs
----------

Run Distortion-Interaction/Activation-Strain analysis for reaction mechanisms.

DI-AS Specific Options
^^^^^^^^^^^^^^^^^^^^^

.. list-table:: DI-AS Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-i, --fragment-indices <string>``
     - string
     - Indices of one fragment for DI-AS analysis (required)
   * - ``-n, --every-n-points <int>``
     - int
     - Every nth points along the IRC file to prepare for DI-AS analysis (default=3)
   * - ``-s, --solv/--no-solv``
     - bool
     - Turn on/off solvent for DI-AS job calculations (default=False)
   * - ``-m, --mode``
     - string
     - Mode of DI-AS analysis. Options: irc, ts (default=irc)
   * - ``-c1, --charge-of-fragment1 <int>``
     - int
     - Charge of fragment 1 (default=None)
   * - ``-m1, --multiplicity-of-fragment1 <int>``
     - int
     - Multiplicity of fragment 1 (default=None)
   * - ``-c2, --charge-of-fragment2 <int>``
     - int
     - Charge of fragment 2 (default=None)
   * - ``-m2, --multiplicity-of-fragment2 <int>``
     - int
     - Multiplicity of fragment 2 (default=None)

Basic Usage
^^^^^^^^^^^

* **Basic DI-AS analysis on IRC**:

    .. code-block:: console

        chemsmart sub gaussian -p dias_analysis -f irc_output.log dias -i "1-10"

* **DI-AS analysis with solvent**:

    .. code-block:: console

        chemsmart sub gaussian -p dias_solv -f irc_output.log dias -i "1-15" -s -n 5

* **DI-AS analysis on transition state**:

    .. code-block:: console

        chemsmart sub gaussian -p dias_ts -f ts_optimized.log dias -i "1-12" -m ts

RESP Jobs
---------

Run RESP (Restrained Electrostatic Potential) charge fitting calculations.

Basic Usage
^^^^^^^^^^^

* **Basic RESP calculation**:

    .. code-block:: console

        chemsmart sub gaussian -p resp_calc -f molecule.xyz resp

**Note**: This creates an input file with fixed route for RESP job: ``HF/6-31+G(d) SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=17,6/50=1)``

NCI Jobs
--------

Run Non-Covalent Interaction analysis for intermolecular interactions.

Basic Usage
^^^^^^^^^^^

* **Basic NCI calculation**:

    .. code-block:: console

        chemsmart sub gaussian -p nci_analysis -f complex.xyz nci

TD-DFT Jobs
-----------

Run time-dependent DFT calculations for excited state properties.

TD-DFT Specific Options
^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: TD-DFT Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-s, --states``
     - string
     - States for closed-shell singlet systems. Options: 'singlets', 'triplets', '50-50' (default=singlets)
   * - ``-r, --root <int>``
     - int
     - Specifies the "state of interest". The default is the first excited state (N=1) (default=1)
   * - ``-n, --nstates <int>``
     - int
     - Solve for M states. If 50-50, this gives the number of each type of state to solve (default=3)
   * - ``-e, --eqsolv <string>``
     - string
     - Whether to perform equilibrium or non-equilibrium PCM solvation (default=None)

Basic Usage
^^^^^^^^^^^

* **Basic TD-DFT calculation**:

    .. code-block:: console

        chemsmart sub gaussian -p td_calc -f molecule.xyz td

* **TD-DFT with specific states**:

    .. code-block:: console

        chemsmart sub gaussian -p td_triplets -f molecule.xyz td -s triplets -n 5

* **TD-DFT with 50-50 singlet-triplet mix**:

    .. code-block:: console

        chemsmart sub gaussian -p td_mixed -f molecule.xyz td -s 50-50 -n 4

WBI Jobs
--------

Run Wiberg Bond Index calculations for bond analysis.

Basic Usage
^^^^^^^^^^^

* **Basic WBI calculation**:

    .. code-block:: console

        chemsmart sub gaussian -p wbi_analysis -f molecule.xyz wbi

* **WBI with NBO analysis**:

    .. code-block:: console

        chemsmart sub gaussian -p nbo_wbi -f optimized.log wbi

