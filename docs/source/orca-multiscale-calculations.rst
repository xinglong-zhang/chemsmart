Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

##########################################
 ORCA QM/MM Multiscale Calculations Guide
##########################################

CHEMSMART provides comprehensive tools for multiscale QM/MM calculations using ORCA. This section covers various QM/MM
schemes including additive, subtractive ONIOM, and crystal QM/MM methods for large molecular systems.

**********
 Overview
**********

ORCA QM/MM calculations allow you to treat large systems by combining quantum mechanical (QM) methods for chemically
important regions with molecular mechanical (MM) force fields for the environment. CHEMSMART supports five types of
multiscale calculations:

#. **Additive QM/MM** - Traditional QM/MM with electrostatic embedding
#. **Subtractive QM/QM2** - Two-layer ONIOM scheme with different QM methods
#. **Subtractive QM/QM2/MM** - Three-layer ONIOM with QM, QM2, and MM layers
#. **MOL-CRYSTAL-QMMM** - QM/MM for molecular crystals
#. **IONIC-CRYSTAL-QMMM** - QM/MM for semiconductors and insulators

*********************
 Basic QM/MM Command
*********************

The basic command structure for ORCA QM/MM calculations is:

.. code:: console

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] <JOBTYPE> qmmm [QMMM_OPTIONS]

where ``<JOBTYPE>`` is one of ``opt``, ``ts``, ``sp``, ``scan``, ``modred``, ``qrc``, or ``neb``. Or using the ``run``
command with project settings:

.. code:: console

   chemsmart run [OPTIONS] orca -p qmmm -f <structure_file> opt qmmm [QMMM_OPTIONS]

Example with the ``run`` command:

.. code:: bash

   chemsmart run --no-scratch --fake orca \
     -p qmmm \
     -f tests/data/StructuresTests/xyz/crest_best.xyz \
     opt qmmm \
     -j QM/QM2/MM \
     -lm system.ORCAFF.prms \
     -ha 1-15 \
     -ct 0 \
     -mt 1 \
     -ch 0 \
     -mh 1

**Key Options:**

-  ``--no-scratch``: Disable scratch for this run (overrides YAML / class defaults; see :ref:`scratch-behavior`)
-  ``--fake``: Dry run mode (don't actually submit)
-  ``-p qmmm``: Load settings from ``~/.chemsmart/orca/qmmm.yaml``
-  ``-f``: Input structure file
-  ``opt qmmm``: Parent job type plus nested QM/MM subcommand

************************
 QM/MM-Specific Options
************************

Job Type and Theory Level
=========================

.. list-table:: Job Type and Method Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  Choice
      -  Multiscale calculation type: QMMM, QM/QM2, QM/QM2/MM, MOL-CRYSTAL-QMMM, IONIC-CRYSTAL-QMMM

   -  -  ``-hx, --high-level-functional``
      -  string
      -  DFT functional for high-level (QM) region (e.g., B3LYP, PBE0)

   -  -  ``-hb, --high-level-basis``
      -  string
      -  Basis set for high-level (QM) region (e.g., def2-SVP, def2-TZVP)

   -  -  ``-ix, --intermediate-level-functional``
      -  string
      -  DFT functional for intermediate-level (QM2) region

   -  -  ``-ib, --intermediate-level-basis``
      -  string
      -  Basis set for intermediate-level (QM2) region

   -  -  ``-im, --intermediate-level-method``
      -  string
      -  Built-in method for intermediate-level (QM2) region (XTB, HF-3C, PBEH-3C, R2SCAN-3C, PM3, AM1)

   -  -  ``-lm, --low-level-method``
      -  string
      -  ORCA force-field parameter file written as ``ORCAFFFilename`` (e.g. ``system.ORCAFF.prms``). Alias:
         ``--low-level-force-field``.

Atom Partitioning
=================

.. list-table:: Atom Partition Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-ha, --high-level-atoms``
      -  string
      -  High-level (QM) atom indices (e.g., '1-15,20' or '1:15 20')

   -  -  ``-ia, --intermediate-level-atoms``
      -  string
      -  Intermediate-level (QM2) atom indices (e.g., '16-30')

   -  -  ``-a, --active-atoms``
      -  string
      -  Active atom indices for optimization

Charge and Multiplicity
=======================

High-level (``-ch`` / ``-mh``) charge and multiplicity are **required**. They are written to the ORCA coordinate section
(``* xyz`` line). Total and intermediate values are written to the ``%qmmm`` block and are not used as a fallback for
the coordinate section.

.. list-table:: Charge and Multiplicity Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-ch, --charge-high``
      -  int
      -  Required. High-level (QM) region charge for the ``* xyz`` line

   -  -  ``-mh, --mult-high``
      -  int
      -  Required. High-level (QM) region multiplicity for the ``* xyz`` line

   -  -  ``-ci, --charge-intermediate``
      -  int
      -  Intermediate layer charge for the ``%qmmm`` block (QM/QM2/MM)

   -  -  ``-mi, --mult-intermediate``
      -  int
      -  Intermediate layer multiplicity for the ``%qmmm`` block

   -  -  ``-ct, --charge-total``
      -  int
      -  Total system charge for the ``%qmmm`` block

   -  -  ``-mt, --mult-total``
      -  int
      -  Total system multiplicity for the ``%qmmm`` block

Advanced QM/MM Options
======================

.. list-table:: Advanced QM/MM Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-s, --intermediate-level-solvation``
      -  string
      -  Solvation model for intermediate-level (QM2) region (CPCM, SMD, ALPB)

   -  -  ``-e, --embedding-type``
      -  string
      -  Embedding type: electronic or mechanical

   -  -  ``-h, --high-level-h-bond-length``
      -  string
      -  Custom high-level-H bond lengths, e.g. ``"{'C_H': 1.09, 'N_H': 1.01}"``.

   -  -  ``-d, --delete-la-double-counting``
      -  bool
      -  Remove bend/torsion double counting

   -  -  ``-db, --delete-la-bond-double-counting-atoms``
      -  bool
      -  Remove bond double counting

Optimization Controls
=====================

.. list-table:: Optimization Control Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-o, --optregion-fixed-atoms``
      -  string
      -  Fixed atom indices in optimization

   -  -  ``-ua, --use-active-info-from-pbc``
      -  string
      -  Use active atom info from PDB file

Crystal QM/MM Options
=====================

.. list-table:: Crystal QM/MM Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-cc, --conv-charges``
      -  bool
      -  Use converged charges for crystal QM/MM

   -  -  ``-xn, --conv-charges-max-n-cycles``
      -  int
      -  Max cycles for charge convergence

   -  -  ``-t, --conv-charges-conv-thresh``
      -  float
      -  Charge convergence threshold

   -  -  ``-sc, --scale-formal-charge-mm-atom``
      -  float
      -  MM atom charge scaling factor

   -  -  ``-nc, --n-unit-cell-atoms``
      -  int
      -  Atoms per unit cell (MOL-CRYSTAL-QMMM)

   -  -  ``-ecp, --ecp-layer-ecp``
      -  string
      -  ECP type for boundary region

   -  -  ``-ecpn, --ecp-layer``
      -  int
      -  Number of ECP layers around QM region

   -  -  ``-sc2, --scale-formal-charge-ecp-atom``
      -  float
      -  ECP atom charge scaling factor

****************
 Usage Examples
****************

Additive QM/MM
==============

Basic additive QM/MM calculation with B3LYP for the QM region and an ORCA force-field parameter file for MM:

.. code:: bash

   chemsmart sub orca -p protein_qmmm -f protein.pdb opt qmmm \
       -j QMMM \
       -hx B3LYP -hb def2-SVP \
       -lm protein.ORCAFF.prms \
       -ha 1-20 \
       -ch 0 -mh 1 \
       -ct 0 -mt 1

Subtractive QM/QM2 ONIOM
========================

Two-layer ONIOM calculation with DFT for high level and semi-empirical for intermediate level:

.. code:: console

   chemsmart sub orca -p enzyme_oniom -f enzyme.xyz opt qmmm \
     -j QM/QM2 \
     -hx B3LYP -hb def2-TZVP \
     -ix HF -ib STO-3G \
     -ha 1-15 -ia 16-50 \
     -ch 0 -mh 1 \
     -ci 0 -mi 1

Three-Layer QM/QM2/MM ONIOM
===========================

Three-layer ONIOM with DFT, XTB, and MM:

.. code:: console

   chemsmart sub orca -p complex_system -f system.pdb opt qmmm \
     -j QM/QM2/MM \
     -hx B3LYP -hb def2-SVP \
     -im XTB \
     -lm system.ORCAFF.prms \
     -ha 1-10 -ia 11-30 \
     -ch 0 -mh 1 \
     -ci 0 -mi 1 \
     -ct 0 -mt 1

Crystal QM/MM for Molecular Crystals
====================================

QM/MM calculation for a molecular crystal:

.. code:: console

   chemsmart sub orca -p molecular_crystal -f crystal.cif opt qmmm \
     -j MOL-CRYSTAL-QMMM \
     -hx PBE -hb def2-SVP \
     -ha 1-20 \
     -ch 0 -mh 1 \
     -nc 50 -cc true -xn 30

Advanced QM/MM with Custom Settings
===================================

QM/MM with custom bond lengths and embedding options:

.. code:: console

   chemsmart sub orca -p advanced_qmmm -f system.xyz opt qmmm \
     -j QMMM \
     -hx M06-2X -hb def2-TZVP \
     -lm system.ORCAFF.prms \
     -ha 1-25 \
     -ch -1 -mh 2 \
     -ct -1 -mt 2 \
     -e electronic \
     -h "{'C_H': 1.09, 'N_H': 1.01}" \
     -d true

***********************
 Project Configuration
***********************

You can configure QM/MM settings in a YAML file to avoid repetitive CLI options. This is especially useful for running
multiple jobs with similar settings. Charge and multiplicity are not set in the project YAML; pass them on the command
line (``-ch``, ``-mh``, ``-ct``, ``-mt``, and for QM/QM2 jobs ``-ci``, ``-mi``).

YAML Configuration Location
===========================

Create a YAML file at ``~/.chemsmart/orca/qmmm.yaml``:

Basic 2-Layer QM/MM Configuration
=================================

.. code:: yaml

   # ~/.chemsmart/orca/qmmm.yaml - Basic QM/MM setup
   qmmm:
     jobtype: "QMMM"
     high_level_functional: "B3LYP"
     high_level_basis: "def2-SVP"
     low_level_method: "system.ORCAFF.prms"
     embedding_type: "Electronic"
     freq: false

3-Layer ONIOM Configuration
===========================

.. code:: yaml

   # ~/.chemsmart/orca/qmmm.yaml - 3-Layer ONIOM
   qmmm:
     jobtype: "QM/QM2/MM"

     # High-level (QM) region
     high_level_functional: "PBE0"
     high_level_basis: "def2-TZVP"

     # Intermediate-level (QM2) region
     intermediate_level_functional: "B3LYP"
     intermediate_level_basis: "def2-SVP"
     # Or use a built-in method:
     # intermediate_level_method: "XTB"

     # Low-level (MM) region
     low_level_method: "system.ORCAFF.prms"

     # Additional options
     freq: false
     embedding_type: "Electronic"
     delete_la_double_counting: true

Using Built-in Methods
======================

.. code:: yaml

   # ~/.chemsmart/orca/qmmm.yaml - Using XTB for intermediate level
   qmmm:
     jobtype: "QM/QM2"
     high_level_functional: "PBE0"
     high_level_basis: "def2-TZVP"
     intermediate_level_method: "XTB"  # Built-in method

**Supported Built-in Methods:**

-  ``XTB``, ``XTB0``, ``XTB1`` - xTB semi-empirical methods
-  ``HF-3C`` - Hartree-Fock with 3 corrections
-  ``PBEH-3C`` - PBEh-3c composite method
-  ``R2SCAN-3C`` - r²SCAN-3c composite method
-  ``PM3``, ``AM1`` - Semi-empirical methods

Using YAML Configuration
========================

Once you have a YAML file configured, use it with the ``-p`` flag:

.. code:: console

   # Theory from YAML; charge and multiplicity on the CLI
   chemsmart run orca -p qmmm -f system.pdb opt qmmm \
     -ha 1-20 \
     -ch 0 -mh 1 -ct 0 -mt 1

   # Override YAML settings with CLI options
   chemsmart run orca -p qmmm -f system.pdb opt qmmm \
     -ha 1-20 \
     -hx M06-2X \
     -ch 0 -mh 1 -ct 0 -mt 1

**Settings Priority:**

#. CLI options (highest priority)
#. YAML configuration file
#. Default settings (lowest priority)

Complete YAML Example
=====================

.. code:: yaml

   # ~/.chemsmart/orca/qmmm.yaml - Complete example
   qmmm:
     jobtype: "QM/QM2/MM"

     # Theory levels
     high_level_functional: "PBE0"
     high_level_basis: "def2-TZVP"
     intermediate_level_method: "XTB"
     low_level_method: "system.ORCAFF.prms"

     # Atom partitioning (optional, can be set via CLI)
     # high_level_atoms: [1, 2, 3, 4, 5]
     # intermediate_level_atoms: [6, 7, 8, 9, 10]

     # QM/MM options
     embedding_type: "Electronic"
     intermediate_level_solvation: "CPCM(Water)"
     delete_la_double_counting: true

     # Job control
     freq: false

.. code:: console

   chemsmart run orca -p qmmm -f system.pdb opt qmmm \
     -ha 1-20 \
     -ch 0 -mh 1 \
     -ct 0 -mt 1

************
 Next Steps
************

For more advanced QM/MM workflows, see:

-  **Project Configuration**: Set up custom QM/MM project settings
-  **Server Configuration**: Configure HPC settings for large QM/MM jobs
-  **Force Field Setup**: Prepare MM parameter files
-  **Analysis Tools**: Post-process QM/MM results

.. tip::

   QM/MM calculations can be computationally demanding. Consider using HPC clusters with adequate memory and CPU
   resources, especially for large systems or crystal QM/MM calculations.

.. warning::

   Ensure proper atom partitioning between QM and MM regions. The QM region should include chemically important areas
   like active sites, reaction centers, or defects in crystals.

*****************
 Quick Reference
*****************

CLI Option Summary
==================

``qmmm`` is nested under a parent job type (``opt``, ``ts``, ``sp``, ``scan``, ``modred``, ``qrc``, or ``neb``).

.. list-table:: Key CLI Options
   :header-rows: 1
   :widths: 15 30 55

   -  -  Short
      -  Long
      -  Description

   -  -  ``-hx``
      -  ``--high-level-functional``
      -  High-level (QM) functional

   -  -  ``-hb``
      -  ``--high-level-basis``
      -  High-level (QM) basis set

   -  -  ``-ha``
      -  ``--high-level-atoms``
      -  High-level atom indices

   -  -  ``-ix``
      -  ``--intermediate-level-functional``
      -  Intermediate (QM2) functional

   -  -  ``-ib``
      -  ``--intermediate-level-basis``
      -  Intermediate (QM2) basis set

   -  -  ``-im``
      -  ``--intermediate-level-method``
      -  Intermediate built-in method

   -  -  ``-ia``
      -  ``--intermediate-level-atoms``
      -  Intermediate atom indices

   -  -  ``-ci``
      -  ``--charge-intermediate``
      -  Intermediate charge

   -  -  ``-mi``
      -  ``--mult-intermediate``
      -  Intermediate multiplicity

   -  -  ``-lm``
      -  ``--low-level-method``
      -  ORCA force-field parameter file (``ORCAFFFilename``). Alias: ``--low-level-force-field``.

Intermediate Layer Options
==========================

Use ``-ix``/``-ib``/``-im``/``-ia``/``-ci``/``-mi`` for the intermediate (QM2) functional, basis, built-in method,
atoms, charge, and multiplicity. YAML keys are ``intermediate_level_*``.

*******************************
 ORCAQMMMJobSettings Reference
*******************************

The ``ORCAQMMMJobSettings`` class provides comprehensive configuration options for ORCA multiscale calculations with
enhanced equality comparison support. This class now includes proper ``__eq__`` method implementation for accurate
settings comparison during job validation and configuration management.

Core Configuration Parameters
=============================

Job Type and Methods
--------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``jobtype``
      -  str
      -  Calculation type (QMMM, QM/QM2, QM/QM2/MM, MOL-CRYSTAL-QMMM, IONIC-CRYSTAL-QMMM)

   -  -  ``high_level_functional``
      -  str
      -  DFT functional for high-level (QM) region (B3LYP, PBE0, etc.)

   -  -  ``high_level_basis``
      -  str
      -  Basis set for high-level (QM) region (def2-SVP, def2-TZVP, etc.)

   -  -  ``intermediate_level_functional``
      -  str
      -  DFT functional for intermediate-level (QM2) layer

   -  -  ``intermediate_level_basis``
      -  str
      -  Basis set for intermediate-level (QM2) layer

   -  -  ``intermediate_level_method``
      -  str
      -  Built-in method for intermediate-level (QM2) (XTB, HF-3C, PBEH-3C, R2SCAN-3C, PM3, AM1)

   -  -  ``low_level_method``
      -  str
      -  ORCA force-field parameter file for ``ORCAFFFilename`` (e.g. ``system.ORCAFF.prms``). YAML also accepts
         ``low_level_force_field``.

Settings Comparison and Validation
----------------------------------

The ``ORCAQMMMJobSettings`` class now supports robust equality comparison including:

-  **Complete attribute comparison**: All QMMM-specific parameters are included in equality checks
-  **Configuration validation**: Enables accurate detection of setting changes during job execution
-  **Type safety**: Proper type checking prevents invalid comparisons
-  **Inheritance support**: Maintains compatibility with parent ORCA settings class

This enhancement improves reliability when: - Comparing job configurations across multiple runs - Validating project
settings consistency - Detecting configuration changes in automated workflows - Debugging multiscale calculation setup
issues

Settings Atom Partitioning
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``high_level_atoms``
      -  list/str
      -  High-level (QM) region atom indices (e.g., [1,2,3] or "1-10,15")

   -  -  ``intermediate_level_atoms``
      -  list/str
      -  Intermediate-level (QM2) region atom indices

   -  -  ``active_atoms``
      -  list/str
      -  Active atoms for optimization (default: all atoms)

Settings Charge and Multiplicity
--------------------------------

``charge_high`` and ``mult_high`` are required for input generation: they are written to the coordinate section.
``charge_total`` / ``charge_intermediate`` (and their multiplicities) belong in the ``%qmmm`` block and do not
substitute for the high-level values.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``charge_high``
      -  int
      -  Required. High-level (QM) region charge (``* xyz`` line)

   -  -  ``mult_high``
      -  int
      -  Required. High-level (QM) region spin multiplicity (``* xyz`` line)

   -  -  ``charge_intermediate``
      -  int
      -  Intermediate layer charge in the ``%qmmm`` block (QM/QM2/MM)

   -  -  ``mult_intermediate``
      -  int
      -  Intermediate layer multiplicity in the ``%qmmm`` block

   -  -  ``charge_total``
      -  int
      -  Total system charge in the ``%qmmm`` block

   -  -  ``mult_total``
      -  int
      -  Total system multiplicity in the ``%qmmm`` block

QM/MM Advanced Options
======================

Link Atom and Interface Control
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``qm_h_bond_length``
      -  dict
      -  Custom QM-H bond distances {(atom1, atom2): distance}

   -  -  ``delete_la_double_counting``
      -  bool
      -  Remove bend/torsion double counting for link atoms

   -  -  ``delete_la_bond_double_counting_atoms``
      -  bool
      -  Remove bond double counting for link atoms

   -  -  ``embedding_type``
      -  str
      -  Electronic or mechanical embedding ("electronic"/"mechanical")

Solvation and Environment
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``qm2_solvation``
      -  str
      -  Solvation model for QM2 region (CPCM, SMD, ALPB(Water))

Crystal QM/MM Parameters
========================

Charge Convergence Control
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``conv_charges``
      -  bool
      -  Use converged charges for crystal QM/MM

   -  -  ``conv_charges_max_n_cycles``
      -  int
      -  Maximum cycles for charge convergence

   -  -  ``conv_charges_conv_thresh``
      -  float
      -  Convergence threshold for charges

   -  -  ``scale_formal_charge_mm_atom``
      -  float
      -  Scaling factor for MM atomic charges

Crystal Structure Parameters
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``n_unit_cell_atoms``
      -  int
      -  Number of atoms per unit cell (MOL-CRYSTAL-QMMM)

   -  -  ``ecp_layer_ecp``
      -  str
      -  ECP type for boundary region (IONIC-CRYSTAL-QMMM)

   -  -  ``ecp_layer``
      -  int
      -  Number of ECP layers around QM region

   -  -  ``scale_formal_charge_ecp_atom``
      -  float
      -  Scaling factor for ECP atomic charges

Optimization Controls
=====================

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``optregion_fixed_atoms``
      -  list/str
      -  Fixed atom indices in geometry optimization

   -  -  ``use_active_info_from_pbc``
      -  bool
      -  Use active atom information from PDB file

*****************
 Troubleshooting
*****************

Common Issues and Solutions
===========================

**Settings Comparison Errors**
   With the enhanced equality comparison, you can now reliably detect when settings have changed between job runs. The
   simplified comparison (without detailed difference logging) ensures faster execution.

**QM/MM Interface Optimization**
   If the QM/MM interface is problematic, adjust custom bond lengths or reconsider atom partitioning. The improved
   settings validation helps ensure consistent interface parameter application.

**Force Field and Method Compatibility**
   Ensure your system has required force field files and that QM methods are compatible. The enhanced error checking
   provides clearer messages about missing dependencies.

**Crystal QM/MM Setup**
   For crystal calculations, ensure proper unit cell definition and charge convergence settings. The validation methods
   help detect incompatible parameter combinations.

**Performance Optimization**
   -  Use appropriate basis sets for each layer (larger for QM, smaller for QM2)
   -  Consider semi-empirical methods for large QM2 regions
   -  Balance QM region size with available computational resources
   -  Use mechanical embedding for very large systems to reduce computational cost

For additional support and examples, see the CHEMSMART community forums and documentation repository.
