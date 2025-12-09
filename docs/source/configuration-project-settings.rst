##################
 Project Settings
##################

Configure project-specific settings for Gaussian and ORCA calculations.

***************************
 Gaussian Project Settings
***************************

The ``~/.chemsmart/gaussian/`` directory contains project settings files specifying DFT functionals, basis sets, and
other calculation parameters.

Example project file (``~/.chemsmart/gaussian/test.yaml``):

.. code:: yaml

   gas:
     functional: m062x
     basis: def2svp
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

Settings behavior:

-  ``gas``: Used for geometry optimization, transition state searches, etc.
-  ``solv``: Used for single point calculations.
-  ``td``: Used for TD-DFT calculations.

To run all calculations with solvent, set ``gas: Null``:

.. code:: yaml

   gas: Null

***********************
 ORCA Project Settings
***********************

The ``~/.chemsmart/orca/`` directory contains ORCA project settings files.

Example project file (``~/.chemsmart/orca/test.yaml``):

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

This runs gas-phase jobs with M062X/def2-SVP and single point calculations with DLPNO-CCSD(T)/CBS using cc-pVDZ/cc-pVTZ
extrapolation in SMD(toluene).

*******************
 Scratch Directory
*******************

Set up scratch directories for Gaussian and ORCA jobs:

.. code:: bash

   ln -s /path/to/scratch/ ~/scratch

.. note::

   If ``freq: False`` is not set, frequency calculations are performed by default for all geometry optimization jobs.
