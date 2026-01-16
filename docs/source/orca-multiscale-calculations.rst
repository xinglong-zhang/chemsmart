Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

##########################################
 ORCA QM/MM Multiscale Calculations Guide
##########################################

ChemSmart provides comprehensive tools for multiscale QM/MM calculations using ORCA. This section covers various QM/MM
schemes including additive, subtractive ONIOM, and crystal QM/MM methods for large molecular systems.

**********
 Overview
**********

ORCA QM/MM calculations allow you to treat large systems by combining quantum mechanical (QM) methods for chemically
important regions with molecular mechanical (MM) force fields for the environment. ChemSmart supports five types of
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

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] qmmm [QMMM_OPTIONS]

Or using the ``run`` command with project settings:

.. code:: console

   chemsmart run [OPTIONS] orca -p qmmm -f <structure_file> qmmm [QMMM_OPTIONS]

Example with the ``run`` command:

.. code:: bash

   chemsmart run --no-scratch --fake orca \
     -p qmmm \
     -f tests/data/StructuresTests/xyz/crest_best.xyz \
     qmmm \
     -j QM/QM2/MM \
     -lf amber \
     -ct 0 \
     -mt 1 \
     -ha 1-15 \
     -R

**Key Options:**

-  ``--no-scratch``: Don't use scratch directory
-  ``--fake``: Dry run mode (don't actually submit)
-  ``-p qmmm``: Load settings from ``~/.chemsmart/orca/qmmm.yaml``
-  ``-f``: Input structure file
-  ``-R``: Run the job after setup

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

   -  -  ``-j, --job-type``
      -  Choice
      -  Multiscale calculation type: QMMM, QM/QM2, QM/QM2/MM, MOL-CRYSTAL-QMMM, IONIC-CRYSTAL-QMMM

   -  -  ``-hx, --high-level-functional``
      -  string
      -  DFT functional for high-level (QM) region (e.g., B3LYP, PBE0)

   -  -  ``-hb, --high-level-basis``
      -  string
      -  Basis set for high-level (QM) region (e.g., def2-SVP, def2-TZVP)

   -  -  ``-intx, --intermediate-level-functional``
      -  string
      -  DFT functional for intermediate-level (QM2) region

   -  -  ``-intb, --intermediate-level-basis``
      -  string
      -  Basis set for intermediate-level (QM2) region

   -  -  ``-im, --intermediate-level-method``
      -  string
      -  Built-in method for intermediate-level (QM2) region (XTB, HF-3C, PBEH-3C, R2SCAN-3C, PM3, AM1)

   -  -  ``-lf, --low-level-force-field``
      -  string
      -  Force field for low-level (MM) region (MMFF, AMBER, CHARMM)

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

   -  -  ``-inta, --intermediate-level-atoms``
      -  string
      -  Intermediate-level (QM2) atom indices (e.g., '16-30')

   -  -  ``-a, --active-atoms``
      -  string
      -  Active atom indices for optimization

Charge and Multiplicity
=======================

.. list-table:: Charge and Multiplicity Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-ch, --charge-high``
      -  int
      -  High-level (QM) region charge

   -  -  ``-mh, --mult-high``
      -  string
      -  High-level (QM) region multiplicity

   -  -  ``-mi, --charge-intermediate``
      -  int
      -  Intermediate layer charge (for QM/QM2/MM)

   -  -  ``-minti, --mult-intermediate``
      -  string
      -  Intermediate layer multiplicity

   -  -  ``-ct, --charge-total``
      -  int
      -  Total system charge

   -  -  ``-mt, --mult-total``
      -  string
      -  Total system multiplicity

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

   -  -  ``-h, --qm-h-bond-length``
      -  dict
      -  Custom QM-H bond lengths

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

Basic additive QM/MM calculation with B3LYP for QM region and AMBER force field:

.. code:: bash

   chemsmart sub orca -p protein_qmmm -f protein.pdb qmmm \
       -j QMMM \
       -hx B3LYP -hb def2-SVP \
       -lf amber99 \
       -ha 1-20 \
       -ch 0 -mh 1 \
       -ct 0 -mt 1

Subtractive QM/QM2 ONIOM
========================

Two-layer ONIOM calculation with DFT for high level and semi-empirical for intermediate level:

.. code:: console

   chemsmart sub orca -p enzyme_oniom -f enzyme.xyz qmmm \
     -j QM/QM2 \
     -hx B3LYP -hb def2-TZVP \
     -intx HF -intb STO-3G \
     -ha 1-15 -inta 16-50 \
     -ch 0 -mh 1 \
     -mi 0 -minti 1

Three-Layer QM/QM2/MM ONIOM
===========================

Three-layer ONIOM with DFT, XTB, and MM:

.. code:: console

   chemsmart sub orca -p complex_system -f system.pdb qmmm \
     -j QM/QM2/MM \
     -hx B3LYP -hb def2-SVP \
     -im XTB \
     -lf amber99 \
     -ha 1-10 -inta 11-30 \
     -ch 0 -mh 1 \
     -mi 0 -minti 1 \
     -ct 0 -mt 1

Crystal QM/MM for Molecular Crystals
====================================

QM/MM calculation for a molecular crystal:

.. code:: console

   chemsmart sub orca -p molecular_crystal -f crystal.cif qmmm \
     -j MOL-CRYSTAL-QMMM \
     -hx PBE -hb def2-SVP \
     -ha 1-20 \
     -ch 0 -mh 1 \
     -nc 50 -cc true -xn 30

Advanced QM/MM with Custom Settings
===================================

QM/MM with custom bond lengths and embedding options:

.. code:: console

   chemsmart sub orca -p advanced_qmmm -f system.xyz qmmm \
     -j QMMM \
     -hx M06-2X -hb def2-TZVP \
     -lf charmm36 \
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
multiple jobs with similar settings.

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
     low_level_force_field: "amber"
     charge_total: 0
     mult_total: 1
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
     charge_high: 0
     mult_high: 1

     # Intermediate-level (QM2) region
     intermediate_level_functional: "B3LYP"
     intermediate_level_basis: "def2-SVP"
     # Or use a built-in method:
     # intermediate_level_method: "XTB"
     charge_intermediate: 0
     mult_intermediate: 1

     # Low-level (MM) region
     low_level_force_field: "amber"

     # Total system
     charge_total: 0
     mult_total: 1

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
     charge_total: 0
     mult_total: 1

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

   # Minimal command - all settings from YAML
   chemsmart run orca -p qmmm -f system.pdb qmmm \
     -ha 1-20 \
     -R

   # Override YAML settings with CLI options
   chemsmart run orca -p qmmm -f system.pdb qmmm \
     -ha 1-20 \
     -hx M06-2X \
     -R

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
     low_level_force_field: "amber"

     # Charge and multiplicity
     charge_high: 0
     mult_high: 1
     charge_intermediate: 0
     mult_intermediate: 1
     charge_total: 0
     mult_total: 1

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

   chemsmart run orca -p qmmm -f system.pdb qmmm \
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

Updated CLI options after refactoring from "medium" to "intermediate" naming:

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

   -  -  ``-intx``
      -  ``--intermediate-level-functional``
      -  Intermediate (QM2) functional

   -  -  ``-intb``
      -  ``--intermediate-level-basis``
      -  Intermediate (QM2) basis set

   -  -  ``-im``
      -  ``--intermediate-level-method``
      -  Intermediate built-in method

   -  -  ``-inta``
      -  ``--intermediate-level-atoms``
      -  Intermediate atom indices

   -  -  ``-mi``
      -  ``--charge-intermediate``
      -  Intermediate charge

   -  -  ``-minti``
      -  ``--mult-intermediate``
      -  Intermediate multiplicity

   -  -  ``-lf``
      -  ``--low-level-force-field``
      -  Low-level (MM) force field

Naming Update Summary
=====================

**Old → New:** The "medium" terminology has been updated to "intermediate" throughout:

-  CLI options: ``-mx`` → ``-intx``, ``-mm`` → ``-im``, ``-cm`` → ``-mi``
-  YAML parameters: ``medium_level_*`` → ``intermediate_level_*``
-  Settings attributes: ``charge_medium`` → ``charge_intermediate``

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

   -  -  ``low_level_force_field``
      -  str
      -  Force field for low-level (MM) region (MMFF, AMBER, CHARMM)

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

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``charge_high``
      -  int
      -  High-level (QM) region charge

   -  -  ``mult_high``
      -  int
      -  High-level (QM) region spin multiplicity

   -  -  ``charge_intermediate``
      -  int
      -  Intermediate layer charge (QM2 layer in 3-layer ONIOM)

   -  -  ``mult_intermediate``
      -  int
      -  Intermediate layer multiplicity

   -  -  ``charge_total``
      -  int
      -  Total system charge

   -  -  ``mult_total``
      -  int
      -  Total system multiplicity

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

For additional support and examples, see the ChemSmart community forums and documentation repository.
