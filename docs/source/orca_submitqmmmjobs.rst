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

where ``<JOBTYPE>`` is one of ``opt``, ``ts``, ``sp``, ``scan``, ``modred``, ``qrc``, or ``neb``.

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
      -  Built-in method for intermediate-level (QM2) region (XTB, HF-3C, PBEH-3C)

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
      -  High-level atom indices (e.g., '1-15,20' or '1:15 20')

   -  -  ``-ia, --intermediate-level-atoms``
      -  string
      -  Intermediate-level atom indices (e.g., '16-30')

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
      -  Solvation model for intermediate-level region (CPCM, SMD)

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

.. code:: console

   chemsmart sub orca -p protein_qmmm -f protein.pdb opt qmmm -j QMMM -hx B3LYP -hb def2-SVP -lm protein.ORCAFF.prms -ha 1-20 -ch 0 -mh 1 -ct 0 -mt 1

Subtractive QM/QM2 ONIOM
========================

Two-layer ONIOM calculation with DFT for high level and semi-empirical for low level:

.. code:: console

   chemsmart sub orca -p enzyme_oniom -f enzyme.xyz opt qmmm -j QM/QM2 -hx B3LYP -hb def2-TZVP -ix HF -ib STO-3G -ha 1-15 -ia 16-50 -ch 0 -mh 1 -ci 0 -mi 1

Three-Layer QM/QM2/MM ONIOM
===========================

Three-layer ONIOM with DFT, semi-empirical, and MM:

.. code:: console

   chemsmart sub orca -p complex_system -f system.pdb opt qmmm -j QM/QM2/MM -hx B3LYP -hb def2-SVP -im HF-3C -lm system.ORCAFF.prms -ha 1-10 -ia 11-30 -ch 0 -mh 1 -ci 0 -mi 1 -ct 0 -mt 1

Crystal QM/MM for Molecular Crystals
====================================

QM/MM calculation for a molecular crystal:

.. code:: console

   chemsmart sub orca -p molecular_crystal -f crystal.cif opt qmmm -j MOL-CRYSTAL-QMMM -hx PBE -hb def2-SVP -ha 1-20 -ch 0 -mh 1 -nc 50 -cc true -xn 30

Advanced QM/MM with Custom Settings
===================================

QM/MM with custom bond lengths and embedding options:

.. code:: console

   chemsmart sub orca -p advanced_qmmm -f system.xyz opt qmmm -j QMMM -hx M06-2X -hb def2-TZVP -lm system.ORCAFF.prms -ha 1-25 -ch -1 -mh 2 -ct -1 -mt 2 -e electronic -h "{'C_H': 1.09, 'N_H': 1.01}" -d true

***********************
 Project Configuration
***********************

You can also configure QM/MM settings in your project YAML file. Create a ``qmmm.yaml`` file in your project directory.
Charge and multiplicity are not set in the project YAML; pass them on the command line.

.. code:: yaml

   # ~/.chemsmart/orca/qmmm.yaml
   qmmm:
     jobtype: "QMMM"
     high_level_functional: B3LYP
     high_level_basis: def2-SVP
     low_level_method: system.ORCAFF.prms
     embedding_type: electronic
     delete_la_double_counting: true

Then use it with the project flag:

.. code:: console

   chemsmart sub orca -p qmmm -f system.pdb opt qmmm -ha 1-20 -ch 0 -mh 1 -ct 0 -mt 1

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
      -  Built-in method for intermediate-level (XTB, HF-3C, PBEH-3C)

   -  -  ``low_level_method``
      -  str
      -  ORCA force-field parameter file for ``ORCAFFFilename`` (e.g. ``system.ORCAFF.prms``). YAML also accepts
         ``low_level_force_field``.

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
      -  High-level region atom indices (e.g., [1,2,3] or "1-10,15")

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
      -  Required. Charge of high-level (QM) region (``* xyz`` line)

   -  -  ``mult_high``
      -  int
      -  Required. Multiplicity of high-level (QM) region (``* xyz`` line)

   -  -  ``charge_intermediate``
      -  int
      -  Intermediate (QM2) layer charge in the ``%qmmm`` block

   -  -  ``mult_intermediate``
      -  int
      -  Intermediate (QM2) layer multiplicity in the ``%qmmm`` block

   -  -  ``charge_total``
      -  int
      -  Total system charge in the ``%qmmm`` block

   -  -  ``mult_total``
      -  int
      -  Total system multiplicity in the ``%qmmm`` block

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

   -  -  ``intermediate_level_solvation``
      -  str
      -  Solvation model for intermediate-level (CPCM, SMD, etc.)

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

   -  -  ``high_level_h_bond_length``
      -  dict
      -  Custom high-level-H bond lengths {(atom1, atom2): length}

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

Create project-specific QM/MM settings in YAML format. Charge and multiplicity are not set in the project YAML; pass
them on the command line (``-ch``, ``-mh``, ``-ct``, ``-mt``).

**Basic QM/MM Configuration** (``~/.chemsmart/orca/qmmm.yaml``):

.. code:: yaml

   # Basic additive QM/MM settings
   qmmm:
     jobtype: "QMMM"
     high_level_functional: "B3LYP"
     high_level_basis: "def2-SVP"
     low_level_method: "system.ORCAFF.prms"
     embedding_type: "electronic"
     delete_la_double_counting: true

**ONIOM Configuration** (``~/.chemsmart/orca/oniom.yaml``):

.. code:: yaml

   # Three-layer ONIOM settings
   qmmm:
     jobtype: "QM/QM2/MM"
     high_level_functional: "B3LYP"
     high_level_basis: "def2-TZVP"
     intermediate_level_functional: "HF"
     intermediate_level_basis: "STO-3G"
     low_level_method: "system.ORCAFF.prms"
     embedding_type: "electronic"

**Crystal QM/MM Configuration** (``~/.chemsmart/orca/crystal.yaml``):

.. code:: yaml

   # Molecular crystal QM/MM settings
   qmmm:
     jobtype: "MOL-CRYSTAL-QMMM"
     high_level_functional: "PBE"
     high_level_basis: "def2-SVP"
     n_unit_cell_atoms: 50
     conv_charges: true
     conv_charges_max_n_cycles: 30
     conv_charges_conv_thresh: 0.01

**Using Project Configuration**

Use the YAML configuration with the project flag:

.. code:: console

   chemsmart sub orca -p qmmm -f system.pdb opt qmmm -ha 1-20 -ch 0 -mh 1 -ct 0 -mt 1

Python Configuration
====================

Programmatic configuration using the settings class:

.. code:: python

   from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

   # Create basic QM/MM settings
   qmmm_settings = ORCAQMMMJobSettings(
       jobtype="QMMM",
       high_level_functional="B3LYP",
       high_level_basis="def2-SVP",
       low_level_method="system.ORCAFF.prms",
       high_level_atoms="1-20",
       charge_high=0,
       mult_high=1,
       charge_total=0,
       mult_total=1,
       embedding_type="electronic",
   )

   # ONIOM configuration
   oniom_settings = ORCAQMMMJobSettings(
       jobtype="QM/QM2/MM",
       high_level_functional="B3LYP",
       high_level_basis="def2-TZVP",
       intermediate_level_method="HF-3C",
       low_level_method="system.ORCAFF.prms",
       high_level_atoms=[1, 2, 3, 4, 5],
       intermediate_level_atoms=list(range(6, 21)),
       charge_high=0,
       mult_high=1,
       charge_intermediate=0,
       mult_intermediate=1,
   )

*******************************
 Best Practices and Validation
*******************************

Required Parameters
===================

-  **QM region**: Must specify ``high_level_functional`` and ``high_level_basis`` (unless using built-in methods)
-  **Atom partitioning**: Must define ``high_level_atoms`` for all calculations
-  **High-level charge/mult**: Must specify ``charge_high`` and ``mult_high`` (written to ``* xyz``)
-  **Total/intermediate charges**: Specify ``charge_total`` / ``charge_intermediate`` for the ``%qmmm`` block as needed
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

-  :doc:`orca-cli-options` - General ORCA CLI options
-  :doc:`configuration-project-settings` - Project configuration guide
-  :doc:`orca-multiscale-calculations` - ORCA multiscale / QM/MM calculations
