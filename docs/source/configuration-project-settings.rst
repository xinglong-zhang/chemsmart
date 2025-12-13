##################
 Project Settings
##################

Configure project-specific settings for Gaussian and ORCA calculations.

***************************
 Gaussian Project Settings
***************************

The ``~/.chemsmart/gaussian/`` directory contains project settings files specifying DFT functionals, basis sets, and
other calculation parameters.

Job Type Settings
=================

Project settings are organized by job type, allowing you to specify different computational methods for different phases
of your calculations:

-  ``gas``: Used for geometry optimization, transition state searches, and other gas-phase structural calculations.
-  ``solv``: Used for single point energy calculations in solution.
-  ``td``: Used for TD-DFT (time-dependent DFT) calculations for excited states and UV-Vis spectra.

Example 1: Basic Multi-Phase Settings
=====================================

This example (``~/.chemsmart/gaussian/test.yaml``) demonstrates how to configure different settings for gas-phase
optimizations, solution-phase single points, and TD-DFT calculations:

.. code:: yaml

   gas:
     ab_initio: Null
     functional: m062x
     basis: def2svp
     semiempirical: Null
     solvent_model: smd
     solvent_id: dichloroethane
   solv:
     functional: m062x
     basis: def2tzvp
     freq: False
     solvent_model: smd
     solvent_id: dichloroethane
   td:
     functional: cam-b3lyp
     basis: genecp
     heavy_elements: ['I']
     heavy_elements_basis: def2-SVPD
     light_elements_basis: def2SVP
     freq: False

In this configuration:

-  **Gas phase** optimizations use M062X/def2-SVP with SMD(dichloroethane) implicit solvation
-  **Solution phase** single points use M062X/def2-TZVP with higher basis set for better energies
-  **TD-DFT** calculations use CAM-B3LYP with mixed basis sets (GENECP) for systems containing iodine

Example 2: Mixed Element Basis Sets
===================================

For systems with transition metals or heavy elements, you can specify different basis sets for different elements
(``~/.chemsmart/gaussian/test2.yaml``):

.. code:: yaml

   gas:
     functional: b3lyp
     basis: genecp
     additional_route_parameters: empiricaldispersion=gd3bj
     heavy_elements: ['Pd', 'Ag', 'Br', 'Cu', 'Mn']
     heavy_elements_basis: def2-TZVPPD
     light_elements_basis: def2-SVP
   solv:
     functional: b3lyp
     freq: False
     basis: def2qzvp
     additional_route_parameters: empiricaldispersion=gd3bj
     solvent_model: smd
     solvent_id: TetraHydroFuran

This configuration:

-  Uses B3LYP-D3(BJ) functional with Grimme's D3 dispersion correction
-  Assigns def2-TZVPPD basis to heavy elements (Pd, Ag, Br, Cu, Mn)
-  Assigns def2-SVP basis to light elements (H, C, N, O, etc.)
-  Solution phase calculations use uniform def2-QZVP basis for high accuracy

Example 3: Custom Solvent Parameters
====================================

For non-standard solvents, you can define custom solvent parameters (``~/.chemsmart/gaussian/test3.yaml``):

.. code:: yaml

   gas:
     functional: mn15
     basis: genecp
     heavy_elements: ['Rh','Ag','Br','Cu', 'Pd']
     heavy_elements_basis: def2tzvpd
     light_elements_basis: def2svp
   solv:
     functional: mn15
     basis: def2qzvp
     freq: False
     solvent_model: smd
     solvent_id: generic
     custom_solvent : |
       solventname=1,1,1,3,3,3-hexafluoropropan-2-ol
       eps=16.7
       epsinf=1.625625 !(not really needed for gs properties)
       HBondAcidity=0.77
       HBondBasicity=0.10
       SurfaceTensionAtInterface=23.23
       CarbonAromaticity=0.0
       ElectronegativeHalogenicity=0.60

This example shows how to:

-  Define custom solvent properties (dielectric constant, H-bonding parameters, etc.)
-  Use MN15 functional with mixed basis sets for transition metal systems
-  Set ``solvent_id: generic`` and provide detailed parameters via ``custom_solvent``

Special Settings
================

To run all calculations with solvent (skip gas phase), set ``gas: Null``:

.. code:: yaml

   gas: Null

***********************
 ORCA Project Settings
***********************

The ``~/.chemsmart/orca/`` directory contains ORCA project settings files.

Job Type Settings
=================

Similar to Gaussian, ORCA project settings support multiple job types:

-  ``gas``: Used for geometry optimization and transition state searches in gas phase.
-  ``solv``: Used for single point energy calculations, often with higher-level methods.

Example: DFT Optimization with High-Level Single Points
=======================================================

This example (``~/.chemsmart/orca/test.yaml``) demonstrates a common workflow: geometry optimization with DFT followed
by high-accuracy single point calculations using coupled cluster methods:

.. code:: yaml

   gas:
     functional: M062X
     basis: def2-SVP
   solv:
     ab_initio: DLPNO-CCSD(T)
     functional: Null
     basis: Extrapolate(2/3,cc)
     aux_basis: AutoAux
     defgrid: DEFGRID3
     freq: False
     scf_tol: TightSCF
     scf_algorithm: KDIIS
     scf_maxiter: 500
     mdci_cutoff: Normal
     mdci_density: None
     dipole: False
     solvent_model: SMD
     solvent_id: "toluene"

This configuration:

-  **Gas phase**: Performs geometry optimizations using M062X/def2-SVP DFT method

-  **Solution phase**: Performs high-accuracy single point energies with DLPNO-CCSD(T) using:

   -  Complete basis set (CBS) extrapolation from cc-pVDZ and cc-pVTZ (``Extrapolate(2/3,cc)``)
   -  Automatic auxiliary basis set selection (``AutoAux``)
   -  Tight SCF convergence with KDIIS algorithm
   -  SMD implicit solvation model with toluene as solvent
   -  MDCI (Modified Davidson Configuration Interaction) cutoff parameters

This workflow is efficient for obtaining highly accurate energies on DFT-optimized geometries, commonly used for
thermochemistry and reaction energetics.

Key ORCA-Specific Parameters
============================

-  ``ab_initio``: Specifies post-HF methods (e.g., DLPNO-CCSD(T), CCSD, MP2)
-  ``functional``: DFT functional (set to ``Null`` when using pure ab initio methods)
-  ``basis``: Basis set specification, supports extrapolation schemes
-  ``aux_basis``: Auxiliary basis for RI approximations (``AutoAux`` for automatic selection)
-  ``defgrid``: Integration grid quality (DEFGRID1-3)
-  ``scf_tol``: SCF convergence threshold (``TightSCF``, ``VeryTightSCF``)
-  ``scf_algorithm``: SCF convergence algorithm (``KDIIS``, ``SOSCF``)
-  ``mdci_cutoff``: MDCI method cutoff settings (``Loose``, ``Normal``, ``Tight``)
-  ``mdci_density``: Density treatment in MDCI (must be the string ``"None"``, not YAML null value)

*******************
 Scratch Directory
*******************

Set up scratch directories for Gaussian and ORCA jobs:

.. code:: bash

   ln -s /path/to/scratch/ ~/scratch

.. note::

   If ``freq: False`` is not set, frequency calculations are performed by default for all geometry optimization jobs.
