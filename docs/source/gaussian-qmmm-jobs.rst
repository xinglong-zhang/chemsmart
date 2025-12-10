Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

#########################################
 Gaussian QM/MM ONIOM Calculations Guide
#########################################

ChemSmart provides comprehensive tools for multiscale QM/MM calculations using Gaussian's ONIOM (Our own N-layered
Integrated molecular Orbital and molecular Mechanics) methodology. This section covers various QM/MM schemes for
accurate treatment of large molecular systems including enzymes, organometallic complexes, and materials.

**********
 Overview
**********

Gaussian ONIOM calculations enable efficient treatment of large molecular systems by partitioning them into multiple
layers with different levels of theory. ChemSmart supports flexible ONIOM calculation setups:

**2-Layer ONIOM**: High(QM):Low(MM) - **High level**: Quantum mechanics (DFT, ab initio, MP2, etc.) - **Low level**:
Molecular mechanics force fields (AMBER, UFF, DREIDING)

**3-Layer ONIOM**: High(QM):Medium(QM):Low(MM) - **High level**: High-accuracy quantum mechanics - **Medium level**:
Lower-cost quantum mechanics - **Low level**: Molecular mechanics force fields

Key advantages of ONIOM methodology: - Automatic link atom handling for covalent boundaries - Customizable scale factors
for accurate QM/MM interfaces - Support for multiple charge/multiplicity specifications - Integration with all major
force fields available in Gaussian - Compatibility with various job types (optimization, frequencies, transition states)

*********************
 Basic QM/MM Command
*********************

The basic command structure for Gaussian ONIOM QM/MM calculations is:

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] qmmm [QMMM_OPTIONS]

************************
 QM/MM-Specific Options
************************

Job Type and Theory Levels
==========================

.. list-table:: Job Type and Method Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  Choice
      -  ONIOM job type: sp, opt, freq, ts, irc

   -  -  ``-hx, --high-level-functional``
      -  string
      -  DFT functional for high layer (e.g., B3LYP, M06-2X, wB97X-D)

   -  -  ``-hb, --high-level-basis``
      -  string
      -  Basis set for high layer (e.g., 6-31G*, def2-TZVP, cc-pVTZ)

   -  -  ``-hf, --high-level-force-field``
      -  string
      -  Force field for high layer (if MM, rare)

   -  -  ``-mx, --medium-level-functional``
      -  string
      -  DFT functional for medium layer

   -  -  ``-mb, --medium-level-basis``
      -  string
      -  Basis set for medium layer

   -  -  ``-mf, --medium-level-force-field``
      -  string
      -  Force field for medium layer

   -  -  ``-lx, --low-level-functional``
      -  string
      -  DFT functional for low layer (if QM, uncommon)

   -  -  ``-lb, --low-level-basis``
      -  string
      -  Basis set for low layer (if QM, uncommon)

   -  -  ``-lf, --low-level-force-field``
      -  string
      -  Force field for low layer (AMBER=HardFirst, UFF, DREIDING)

Atom Partitioning
=================

.. list-table:: Layer Assignment Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-ha, --high-level-atoms``
      -  string
      -  Atom indices for high layer (e.g., '1-15,20' or '1:15 20')

   -  -  ``-ma, --medium-level-atoms``
      -  string
      -  Atom indices for medium layer (3-layer only)

   -  -  ``-la, --low-level-atoms``
      -  string
      -  Atom indices for low layer (usually auto-assigned)

   -  -  ``-b, --bonded-atoms``
      -  string
      -  Bonds crossing layer boundaries: e.g., '(1,2),(5,6)'

   -  -  ``-s, --scale-factors``
      -  dict
      -  Custom link atom scale factors: {(atom1,atom2): [low,med,high]}

Charge and Multiplicity
=======================

.. list-table:: Charge and Multiplicity Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-cr, --real-charge``
      -  int
      -  Charge of complete molecular system

   -  -  ``-mr, --real-multiplicity``
      -  int
      -  Spin multiplicity of complete system (2S+1)

   -  -  ``-ci, --int-charge``
      -  int
      -  Charge of high+medium layers (3-layer only)

   -  -  ``-mi, --int-multiplicity``
      -  int
      -  Multiplicity of high+medium layers

   -  -  ``-cm, --model-charge``
      -  int
      -  Charge of high layer only

   -  -  ``-mm, --model-multiplicity``
      -  int
      -  Multiplicity of high layer only

****************
 Usage Examples
****************

2-Layer Enzyme Active Site
==========================

Basic enzyme QM/MM calculation with DFT for active site and AMBER for protein:

.. code:: console

   chemsmart sub gaussian -p enzyme_qmmm -f protein.pdb qmmm -j opt -hx B3LYP -hb 6-31G* -lf AMBER=HardFirst -ha 1-25 -cr 0 -mr 1 -cm 0 -mm 1 -b "(25,26)"

3-Layer Organometallic Catalyst
===============================

Multi-layer calculation with high-accuracy DFT for metal center:

.. code:: console

   chemsmart sub gaussian -p catalyst_oniom -f complex.xyz qmmm -j freq -hx M06-2X -hb def2-TZVP -mx B3LYP -mb 6-31G* -lf UFF -ha 1-10 -ma 11-50 -cr -1 -mr 2 -ci 0 -mi 1 -cm 0 -mm 1 -b "(10,11),(50,51)"

Transition State Search
=======================

ONIOM transition state optimization for enzyme catalysis:

.. code:: console

   chemsmart sub gaussian -p ts_qmmm -f reactant.com qmmm -j ts -hx wB97X-D -hb 6-311++G(d,p) -lf AMBER=HardFirst -ha 1-30 -cr 0 -mr 1 -cm 0 -mm 1 -b "(30,31)" -s "{(30,31): [0.709, 0.709, 0.709]}"

Frequency Analysis
==================

Vibrational analysis of QM/MM optimized structure:

.. code:: console

   chemsmart sub gaussian -p freq_qmmm -f optimized.com qmmm -j freq -hx B3LYP -hb 6-31G* -lf AMBER=HardFirst -ha 1-20 -cr 0 -mr 1 -cm 0 -mm 1

Single Point Energy
===================

High-accuracy single point calculation on QM/MM geometry:

.. code:: console

   chemsmart sub gaussian -p sp_qmmm -f geometry.xyz qmmm -j sp -hx CCSD(T) -hb aug-cc-pVDZ -lf AMBER=HardFirst -ha 1-15 -cr 0 -mr 1 -cm 0 -mm 1

***********************************
 GaussianQMMMJobSettings Reference
***********************************

The ``GaussianQMMMJobSettings`` class provides comprehensive configuration options for Gaussian ONIOM calculations with
enhanced equality comparison support. This class now includes proper ``__eq__`` method implementation for accurate
settings comparison during job validation and configuration management.

Core Configuration Parameters
=============================

Job Type and Theory Levels
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``jobtype``
      -  str
      -  ONIOM calculation type (sp, opt, freq, ts, irc)

   -  -  ``high_level_functional``
      -  str
      -  DFT functional for high layer (B3LYP, M06-2X, wB97X-D, etc.)

   -  -  ``high_level_basis``
      -  str
      -  Basis set for high layer (6-31G*, def2-TZVP, cc-pVTZ, etc.)

   -  -  ``high_level_force_field``
      -  str
      -  Force field for high layer (if MM, uncommon)

   -  -  ``medium_level_functional``
      -  str
      -  DFT functional for medium layer

   -  -  ``medium_level_basis``
      -  str
      -  Basis set for medium layer

   -  -  ``medium_level_force_field``
      -  str
      -  Force field for medium layer

   -  -  ``low_level_functional``
      -  str
      -  DFT functional for low layer (if QM, rare)

   -  -  ``low_level_basis``
      -  str
      -  Basis set for low layer (if QM, rare)

   -  -  ``low_level_force_field``
      -  str
      -  Force field for low layer (AMBER=HardFirst, UFF, DREIDING, etc.)

Settings Comparison and Validation
----------------------------------

The ``GaussianQMMMJobSettings`` class now supports robust equality comparison including:

-  **Complete attribute comparison**: All QMMM-specific parameters are included in equality checks
-  **Configuration validation**: Enables accurate detection of setting changes during job execution
-  **Debugging support**: Detailed difference logging when settings don't match
-  **Type safety**: Proper type checking prevents invalid comparisons

This enhancement improves reliability when: - Comparing job configurations across multiple runs - Validating project
settings consistency - Detecting configuration changes in automated workflows - Debugging job setup issues

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
      -  High layer atom indices (e.g., [1,2,3] or "1-10,15")

   -  -  ``medium_level_atoms``
      -  list/str
      -  Medium layer atom indices (3-layer only)

   -  -  ``low_level_atoms``
      -  list/str
      -  Low layer atom indices (usually auto-assigned)

   -  -  ``bonded_atoms``
      -  list
      -  Bonds crossing boundaries: [(atom1, atom2), (atom3, atom4)]

   -  -  ``scale_factors``
      -  dict
      -  Link atom scale factors: {(atom1, atom2): [low, med, high]}

Settings Charge and Multiplicity
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   -  -  Parameter
      -  Type
      -  Description

   -  -  ``real_charge``
      -  int
      -  Total charge of complete molecular system

   -  -  ``real_multiplicity``
      -  int
      -  Spin multiplicity of complete system (2S+1)

   -  -  ``int_charge``
      -  int
      -  Charge of high+medium layers (3-layer only)

   -  -  ``int_multiplicity``
      -  int
      -  Multiplicity of high+medium layers

   -  -  ``model_charge``
      -  int
      -  Charge of high layer only

   -  -  ``model_multiplicity``
      -  int
      -  Multiplicity of high layer only

Advanced ONIOM Features
=======================

Link Atom Handling
------------------

ONIOM automatically handles covalent bonds crossing layer boundaries by:

-  **Link atom placement**: Hydrogen atoms placed along cut bonds
-  **Scale factor control**: Distance scaling for optimal QM/MM interface
-  **Charge redistribution**: Prevents charge buildup at boundaries
-  **Geometric constraints**: Maintains chemical sensibility

Common scale factors: - **0.709**: Standard C-H bond ratio (recommended for most cases) - **0.724**: C-N bond ratio -
**0.659**: C-O bond ratio - **Custom values**: Can be optimized for specific systems

Force Field Integration
-----------------------

Gaussian supports multiple MM force fields:

.. list-table:: Available Force Fields
   :header-rows: 1
   :widths: 25 75

   -  -  Force Field
      -  Description and Applications
   -  -  ``AMBER=HardFirst``
      -  General-purpose for proteins, nucleic acids, standard organic molecules
   -  -  ``UFF``
      -  Universal force field for any element combination
   -  -  ``DREIDING``
      -  Good for organic and main-group compounds
   -  -  ``MM3``
      -  Accurate for small organic molecules
   -  -  ``AMBER=SoftFirst``
      -  Alternative AMBER parameterization

***********************
 Configuration Methods
***********************

YAML Configuration Files
========================

Create project-specific ONIOM settings in YAML format:

**Basic 2-Layer QM/MM** (``~/.chemsmart/gaussian/qmmm.yaml``):

.. code:: yaml

   # Basic enzyme active site calculation
   jobtype: "opt"
   high_level_functional: "B3LYP"
   high_level_basis: "6-31G*"
   low_level_force_field: "AMBER=HardFirst"
   real_charge: 0
   real_multiplicity: 1
   model_charge: 0
   model_multiplicity: 1

**3-Layer ONIOM Configuration** (``~/.chemsmart/gaussian/oniom3.yaml``):

.. code:: yaml

   # Advanced organometallic catalyst
   jobtype: "freq"
   high_level_functional: "M06-2X"
   high_level_basis: "def2-TZVP"
   medium_level_functional: "B3LYP"
   medium_level_basis: "6-31G*"
   low_level_force_field: "UFF"
   real_charge: -1
   real_multiplicity: 2
   int_charge: 0
   int_multiplicity: 1
   model_charge: 0
   model_multiplicity: 1
   bonded_atoms: [(10, 11), (50, 51)]
   scale_factors:
     (10, 11): [0.709, 0.709, 0.709]
     (50, 51): [0.709, 0.709, 0.709]

Usage with project settings:

.. code:: console

   chemsmart sub gaussian -p qmmm -f complex.xyz qmmm -ha "1-10" -ma "11-50"

*****************
 Troubleshooting
*****************

Common Issues and Solutions
===========================

**Settings Comparison Errors**
   With the enhanced equality comparison, you can now reliably detect when settings have changed between job runs. Check
   the logging output for detailed differences when settings don't match expectations.

**Link Atom Optimization**
   If QM/MM interface is problematic, adjust scale factors or reconsider layer boundaries. The improved settings
   validation helps ensure consistent scale factor application.

**Force Field Compatibility**
   Ensure your system has the required force field files. The enhanced error checking provides clearer messages about
   missing dependencies.

**Performance Optimization**
   -  Use smaller basis sets for medium layers in 3-layer calculations
   -  Consider semi-empirical methods for large QM2 regions
   -  Balance QM region size with computational resources

For additional support and examples, see the ChemSmart community forums and documentation repository.
