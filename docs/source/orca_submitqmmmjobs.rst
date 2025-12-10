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

   -  -  ``-qx, --qm-functional``
      -  string
      -  DFT functional for QM region (e.g., B3LYP, PBE0)

   -  -  ``-qb, --qm-basis``
      -  string
      -  Basis set for QM region (e.g., def2-SVP, def2-TZVP)

   -  -  ``-qx2, --qm2-functional``
      -  string
      -  DFT functional for QM2 region

   -  -  ``-qb2, --qm2-basis``
      -  string
      -  Basis set for QM2 region

   -  -  ``-qm2m, --qm2-methods``
      -  string
      -  Built-in method for QM2 region (XTB, HF-3C, PBEH-3C)

   -  -  ``-mf, --mm-force-field``
      -  string
      -  Force field for MM region (MMFF, AMBER, CHARMM)

Atom Partitioning
=================

.. list-table:: Atom Partition Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-qa, --qm-atoms``
      -  string
      -  QM atom indices (e.g., '1-15,20' or '1:15 20')

   -  -  ``-qa2, --qm2-atoms``
      -  string
      -  QM2 atom indices (e.g., '16-30')

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

   -  -  ``-cq, --charge-qm``
      -  int
      -  QM region charge

   -  -  ``-mq, --mult-qm``
      -  string
      -  QM region multiplicity

   -  -  ``-cm, --charge-medium``
      -  int
      -  Medium layer charge (for QM/QM2/MM)

   -  -  ``-mm, --mult-medium``
      -  string
      -  Medium layer multiplicity

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

   -  -  ``-s, --qm2-solvation``
      -  string
      -  Solvation model for QM2 region (CPCM, SMD)

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

.. code:: console

   chemsmart sub orca -p protein_qmmm -f protein.pdb qmmm -j QMMM -qx B3LYP -qb def2-SVP -mf amber99 -qa 1-20 -cq 0 -mq 1 -ct 0 -mt 1

Subtractive QM/QM2 ONIOM
========================

Two-layer ONIOM calculation with DFT for high level and semi-empirical for low level:

.. code:: console

   chemsmart sub orca -p enzyme_oniom -f enzyme.xyz qmmm -j QM/QM2 -qx B3LYP -qb def2-TZVP -qx2 HF -qb2 STO-3G -qa 1-15 -qa2 16-50 -cq 0 -mq 1 -cm 0 -mm 1

Three-Layer QM/QM2/MM ONIOM
===========================

Three-layer ONIOM with DFT, semi-empirical, and MM:

.. code:: console

   chemsmart sub orca -p complex_system -f system.pdb qmmm -j QM/QM2/MM -qx B3LYP -qb def2-SVP -qm2m HF-3C -mf amber99 -qa 1-10 -qa2 11-30 -cq 0 -mq 1 -cm 0 -mm 1 -ct 0 -mt 1

Crystal QM/MM for Molecular Crystals
====================================

QM/MM calculation for a molecular crystal:

.. code:: console

   chemsmart sub orca -p molecular_crystal -f crystal.cif qmmm -j MOL-CRYSTAL-QMMM -qx PBE -qb def2-SVP -qa 1-20 -cq 0 -mq 1 -nc 50 -cc true -xn 30

Advanced QM/MM with Custom Settings
===================================

QM/MM with custom bond lengths and embedding options:

.. code:: console

   chemsmart sub orca -p advanced_qmmm -f system.xyz qmmm -j QMMM -qx M06-2X -qb def2-TZVP -mf charmm36 -qa 1-25 -cq -1 -mq 2 -ct -1 -mt 2 -e electronic -h "{'C_H': 1.09, 'N_H': 1.01}" -d true

***********************
 Project Configuration
***********************

You can also configure QM/MM settings in your project YAML file. Create a ``qmmm.yaml`` file in your project directory:

.. code:: yaml

   # ~/.chemsmart/orca/qmmm.yaml
   qm_functional: B3LYP
   qm_basis: def2-SVP
   mm_force_field: amber99
   job_type: opt
   embedding_type: electronic
   delete_la_double_counting: true

Then use it with the project flag:

.. code:: console

   chemsmart sub orca -p qmmm -f system.pdb qmmm -qa 1-20 -cq 0 -mq 1 -ct 0 -mt 1

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

*******************************
 ORCAQMMMJobSettings Reference
*******************************

The ``ORCAQMMMJobSettings`` class provides comprehensive configuration options for ORCA multiscale calculations. This
section provides detailed documentation for programmatic configuration.

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

   -  -  ``qm_functional``
      -  str
      -  DFT functional for QM region (B3LYP, PBE0, etc.)

   -  -  ``qm_basis``
      -  str
      -  Basis set for QM region (def2-SVP, def2-TZVP, etc.)

   -  -  ``qm2_functional``
      -  str
      -  DFT functional for QM2 intermediate layer

   -  -  ``qm2_basis``
      -  str
      -  Basis set for QM2 intermediate layer

   -  -  ``qm2_method``
      -  str
      -  Built-in method for QM2 (XTB, HF-3C, PBEH-3C)

   -  -  ``mm_force_field``
      -  str
      -  Force field for MM region (MMFF, AMBER, CHARMM)

Settings Atom Partitioning
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``qm_atoms``
      -  list/str
      -  QM region atom indices (e.g., [1,2,3] or "1-10,15")

   -  -  ``qm2_atoms``
      -  list/str
      -  QM2 region atom indices

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

   -  -  ``charge_qm``
      -  int
      -  Charge of QM region

   -  -  ``mult_qm``
      -  int
      -  Multiplicity of QM region

   -  -  ``charge_medium``
      -  int
      -  Charge of medium system (QM2 layer)

   -  -  ``mult_medium``
      -  int
      -  Multiplicity of medium system

   -  -  ``charge_total``
      -  int
      -  Total system charge

   -  -  ``mult_total``
      -  int
      -  Total system multiplicity

Advanced Settings Options
=========================

Embedding and Interactions
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``embedding_type``
      -  str
      -  Electronic (default) or mechanical embedding

   -  -  ``qm2_solvation``
      -  str
      -  Solvation model for QM2 (CPCM, SMD, etc.)

   -  -  ``delete_la_double_counting``
      -  bool
      -  Remove bend/torsion double counting

   -  -  ``delete_la_bond_double_counting_atoms``
      -  bool
      -  Remove bond double counting

Custom Bond Parameters
----------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``qm_h_bond_length``
      -  dict
      -  Custom QM-H bond lengths {(atom1, atom2): length}

Settings Optimization Controls
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``optregion_fixed_atoms``
      -  list/str
      -  Fixed atom indices during optimization

   -  -  ``use_active_info_from_pbc``
      -  bool
      -  Use active atom info from PDB file

Crystal QM/MM Parameters
========================

For MOL-CRYSTAL-QMMM and IONIC-CRYSTAL-QMMM calculations:

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``conv_charges``
      -  bool
      -  Use converged charges (default: True)

   -  -  ``conv_charges_max_n_cycles``
      -  int
      -  Max charge convergence cycles

   -  -  ``conv_charges_conv_thresh``
      -  float
      -  Charge convergence threshold

   -  -  ``scale_formal_charge_mm_atom``
      -  float
      -  MM atom charge scaling factor

   -  -  ``n_unit_cell_atoms``
      -  int
      -  Atoms per unit cell (required for MOL-CRYSTAL-QMMM)

   -  -  ``ecp_layer_ecp``
      -  str
      -  ECP type for boundary region

   -  -  ``ecp_layer``
      -  int
      -  Number of ECP layers around QM region

   -  -  ``scale_formal_charge_ecp_atom``
      -  float
      -  ECP atom charge scaling factor

***********************
 Configuration Methods
***********************

YAML Configuration Files
========================

Create project-specific QM/MM settings in YAML format:

**Basic QM/MM Configuration** (``~/.chemsmart/orca/qmmm.yaml``):

.. code:: yaml

   # Basic additive QM/MM settings
   jobtype: "QMMM"
   qm_functional: "B3LYP"
   qm_basis: "def2-SVP"
   mm_force_field: "amber99"
   embedding_type: "electronic"
   delete_la_double_counting: true

**ONIOM Configuration** (``~/.chemsmart/orca/oniom.yaml``):

.. code:: yaml

   # Three-layer ONIOM settings
   jobtype: "QM/QM2/MM"
   qm_functional: "B3LYP"
   qm_basis: "def2-TZVP"
   qm2_functional: "HF"
   qm2_basis: "STO-3G"
   mm_force_field: "amber99"
   embedding_type: "electronic"

**Crystal QM/MM Configuration** (``~/.chemsmart/orca/crystal.yaml``):

.. code:: yaml

   # Molecular crystal QM/MM settings
   jobtype: "MOL-CRYSTAL-QMMM"
   qm_functional: "PBE"
   qm_basis: "def2-SVP"
   n_unit_cell_atoms: 50
   conv_charges: true
   conv_charges_max_n_cycles: 30
   conv_charges_conv_thresh: 0.01

**Using Project Configuration**

Use the YAML configuration with the project flag:

.. code:: console

   chemsmart sub orca -p qmmm -f system.pdb qmmm -qa 1-20 -cq 0 -mq 1 -ct 0 -mt 1

Python Configuration
====================

Programmatic configuration using the settings class:

.. code:: python

   from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

   # Create basic QM/MM settings
   qmmm_settings = ORCAQMMMJobSettings(
       jobtype="QMMM",
       qm_functional="B3LYP",
       qm_basis="def2-SVP",
       mm_force_field="amber99",
       qm_atoms="1-20",
       charge_qm=0,
       mult_qm=1,
       charge_total=0,
       mult_total=1,
       embedding_type="electronic",
   )

   # ONIOM configuration
   oniom_settings = ORCAQMMMJobSettings(
       jobtype="QM/QM2/MM",
       qm_functional="B3LYP",
       qm_basis="def2-TZVP",
       qm2_method="HF-3C",
       mm_force_field="charmm36",
       qm_atoms=[1, 2, 3, 4, 5],
       qm2_atoms=list(range(6, 21)),
       charge_qm=0,
       mult_qm=1,
       charge_medium=0,
       mult_medium=1,
   )

*******************************
 Best Practices and Validation
*******************************

Required Parameters
===================

-  **QM region**: Must specify ``qm_functional`` and ``qm_basis`` (unless using built-in methods)
-  **Atom partitioning**: Must define ``qm_atoms`` for all calculations
-  **Charges**: Must specify appropriate charge/multiplicity for each layer
-  **Force fields**: Required for calculations involving MM regions

Validation Rules
================

#. **Method conflicts**: Cannot specify both functional/basis and built-in methods for the same layer
#. **Crystal QM/MM**: Multiplicity should not be specified for crystal calculations
#. **Force field requirements**: MM calculations require valid force field specification
#. **Charge/multiplicity consistency**: All layers must have compatible charge/multiplicity values

Performance Considerations
==========================

-  **QM region size**: Keep QM regions manageable (typically < 100 atoms for routine calculations)
-  **Basis set selection**: Balance accuracy vs. computational cost
-  **Crystal calculations**: Can be very demanding - ensure adequate computational resources
-  **Convergence settings**: Adjust charge convergence parameters for crystal QM/MM

Error Prevention
================

-  **Atom indexing**: Ensure atom indices are valid and don't overlap inappropriately
-  **File dependencies**: Verify force field parameter files exist and are accessible
-  **Resource allocation**: Ensure sufficient memory/CPU for large QM/MM systems
-  **Boundary effects**: Carefully choose QM/MM boundaries to avoid artifacts

*******************
 Advanced Features
*******************

**Expanded Next Steps**

For more advanced QM/MM workflows, see:

-  **Project Configuration**: Set up custom QM/MM project settings
-  **Server Configuration**: Configure HPC settings for large QM/MM jobs
-  **Force Field Setup**: Prepare MM parameter files
-  **Analysis Tools**: Post-process QM/MM results

***************
 API Reference
***************

.. autoclass:: chemsmart.jobs.orca.settings.ORCAQMMMJobSettings
   :members:
   :undoc-members:
   :show-inheritance:

**********
 See Also
**********

-  :doc:`orca_generalcliforallorcajobs` - General ORCA CLI options
-  :doc:`configuration_setuptheprojectsettings` - Project configuration guide
