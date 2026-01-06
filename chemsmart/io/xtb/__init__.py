class XTBRefs:
    """
    Reference data for xTB (extended tight-binding) semi-empirical quantum chemistry software.

    This class provides comprehensive reference lists of xTB input keywords, groups,
    methods, optimization levels, and solvation models. Used for validating xTB input
    files and providing autocomplete functionality.
    """

    XTB_GROUPS = [
        "chrg",  # Set the charge of the molecule
        "spin",  # Set Nalpha-Nbeta of the molecule
        "cma",  # Shift molecule to center of mass
        "constrain",  # Constrain the gradient by appling potentials
        "cube",
        "external",
        "fix",  # Exact fixing
        "gbsa",  # Generalized born (GB) model with solvent accessable surface area
        "gfn",  # GFN Hamiltonian
        "hess",
        "metadyn",
        "md",
        "modef",
        "oniom",
        "opt",
        "path",
        "scan",
        "scc",
        "split",
        "stm",
        "symmetry",
        "thermo",
        "wall",
        "write",
    ]

    XTB_METHODS = ["GFN0", "GFN1", "GFN2", "GFNFF"]

    XTB_OPT_LEVEL = [
        "crude",
        "sloppy",
        "loose",
        "lax",
        "normal",
        "tight",
        "vtight",
        "extreme",
    ]

    XTB_OPT_ENGINE = [
        "rf",  # Approximate Normal Coordinate Rational Function optimizer
        "lbfgs",  # L-BFGS Approximate Normal Coordinate optimizer
        "inertial",  # Fast Inertial Relaxation Engine
    ]

    XTB_SOLVENT_MODELS = [
        "gbsa",  # Generalized born model with solvent accessable surface area contributions
        "alpb",  # analytical linearized Poisson-Boltzmann model
    ]

    @property
    def xtb_groups(self):
        """
        Get xTB input groups/blocks.

        Returns:
            list: Available xTB input groups in lowercase
        """
        return [group.lower() for group in self.XTB_GROUPS]

    @property
    def xtb_methods(self):
        """
        Get xTB calculation methods.

        Returns:
            list: Available xTB methods (GFN0, GFN1, GFN2, GFNFF) in lowercase
        """
        return [method.lower() for method in self.XTB_METHODS]

    @property
    def xtb_opt_levels(self):
        """
        Get xTB optimization convergence levels.

        Returns:
            list: Available optimization levels in lowercase
        """
        return [level.lower() for level in self.XTB_OPT_LEVEL]

    @property
    def xtb_opt_engines(self):
        """
        Get xTB optimization engines/algorithms.

        Returns:
            list: Available optimization engines in lowercase
        """
        return [engine.lower() for engine in self.XTB_OPT_ENGINE]

    @property
    def xtb_solvent_models(self):
        """
        Get xTB implicit solvation models.

        Returns:
            list: Available solvation models in lowercase
        """
        return [model.lower() for model in self.XTB_SOLVENT_MODELS]


# Global constants for convenient access to xTB reference data
xtb_ref = XTBRefs()
XTB_ALL_GROUPS = xtb_ref.xtb_groups
XTB_ALL_METHODS = xtb_ref.xtb_methods
XTB_ALL_OPT_LEVELS = xtb_ref.xtb_opt_levels
XTB_ALL_OPT_ENGINES = xtb_ref.xtb_opt_engines
XTB_ALL_SOLVENT_MODELS = xtb_ref.xtb_solvent_models
