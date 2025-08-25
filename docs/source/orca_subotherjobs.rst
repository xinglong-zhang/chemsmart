Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Other ORCA Jobs
======================

ChemSmart provides additional ORCA job types for specialized calculations and direct input file execution. This section covers single point calculations and direct ORCA input file execution.

Single Point Calculations
-------------------------

Single point calculations compute the energy and properties of a molecule at a fixed geometry without optimization.

Run "chemsmart sub [GENERAL_OPTIONS] orca [ORCA_OPTIONS] sp" to perform single point calculations.

Basic Usage
^^^^^^^^^^^

* **Basic single point calculation**:

    .. code-block:: console

        chemsmart sub orca -p project_name -f input.xyz sp

* **Single point with specific method and basis set**:

    .. code-block:: console

        chemsmart sub orca -p sp_calc -f molecule.xyz -x B3LYP -b def2-TZVP sp

* **Single point with high-level method**:

    .. code-block:: console

        chemsmart sub orca -p ccsd_sp -f molecule.xyz -A "CCSD(T)" -b cc-pVTZ sp

* **Single point for molecular properties**:

    .. code-block:: console

        chemsmart sub orca -p properties -f molecule.xyz -x wB97X-D3 -b def2-TZVP sp


Direct ORCA Input File Execution
--------------------------------

This option allows you to run pre-prepared ORCA input files directly without any modifications by ChemSmart.

Run "chemsmart sub orca -f <filename.inp> inp" to execute ORCA input files directly.

Basic Usage
^^^^^^^^^^^

* **Run existing ORCA input file**:

    .. code-block:: console

        chemsmart sub orca -f calculation.inp inp

* **Run input file with project specification**:

    .. code-block:: console

        chemsmart sub orca -p direct_run -f my_calculation.inp inp

* **Submit multiple input files**:

    .. code-block:: console

        chemsmart sub orca -f "*.inp" inp
