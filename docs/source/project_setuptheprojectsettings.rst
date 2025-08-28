Welcome to the tutorials! We’re thrilled to have you here. Please go through the code examples, and don’t hesitate to
contact our team if you have questions or feedback.

#############################
 Set up the Project Settings
#############################

Chemsmart allows user to execute projects either locally or on a cluster. Users can refer to the following examples to
configure projects as required.

-  The ``~/.chemsmart/gaussian/`` directory contains files related to gaussian project settings, which contain DFT
   functional and basis set etc, that is required to write the input file for running a gaussian job. For example, we
   can specify a test project settings in ``~/.chemsmart/gaussian/test.yaml`` with the following information:

      .. code:: console

         gas:
           functional: m062x  # quotes required for string with spaces
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
           ##solvent_model: smd
           ##solvent_id: DiethylEther

   By default, the ``gas`` phase settings are used for all jobs such as geometry optimization, transition state search
   etc, and the ``solv`` settings are used for single point calculations; the ``td`` settings are used to run TD-DFT
   calculations. One can specify additional project settings in ``~/.chemsmart/gaussian/`` in a similar way to adapt to
   each project that one wishes to run. If setting

      .. code:: console

         gas: Null

   Then all jobs will use settings specified in ``solv``, i.e., all calculations will be run in implicit solvation
   model.

-  The ``~/.chemsmart/orca/`` directory contains files related to ORCA project settings, which contain DFT functional
   and basis set etc, that is required to write the input file for running an ORCA job. For example, we can specify a
   test project settings in ``~/.chemsmart/orca/test.yaml`` with the following information:

      .. code:: console

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

   This will run jobs in the gas phase (geometry and TS opt etc) using M062X/def2-SVP method and run single point with
   solvent correction using DLPNO-CCSD(T)/CBS with cc-pVDZ/cc-pVTZ extrapolation in SMD(toluene), for example. Again,
   users can customize different settings in different ``~/.chemsmart/orca/*project_settings*.yaml`` files to adapt to
   different project requirements.

-  One also need to set up scratch directories where scratch jobs may be run (for Gaussian and ORCA jobs, by default,
   these are run in scratch folder), one may do ``ls -s /path/to/scratch/ ~/scratch``.

.. note::

   If ``freq: False`` is not set in the project settings, frequency calculation will be performed by default for all
   geometry optimization jobs.
