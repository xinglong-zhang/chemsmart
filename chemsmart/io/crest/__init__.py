class CRESTRefs:
    """
    Reference data for CREST (Conformer–Rotamer Ensemble Sampling Tool) software.

    This class provides comprehensive reference lists of CREST input keywords, methods,
    optimization levels, and solvation models. Used for validating CREST input
    files and providing autocomplete functionality.
    """

    CREST_METHODS = ["gfn1", "gfn2", "gfnff", "gfn2//gfnff"]

    # xTB ANCOPT levels accepted by CREST --optlev
    CREST_OPT_LEVELS = [
        "crude",
        "sloppy",
        "loose",
        "lax",
        "normal",
        "tight",
        "vtight",
        "extreme",
    ]

    CREST_JOB_TYPES = ["conformers"]

    CREST_SOLVENT_MODELS = [
        "gbsa",  # generalized born (GB) model with solvent accessible surface (SASA) model
        "alpb",  # analytical linearized Poisson-Boltzmann (ALPB) model
    ]

    CREST_SOLVENT_IDS = [
        "acetone",
        "acetonitrile",
        "aniline",
        "benzaldehyde",
        "benzene",
        "ch2cl2",
        "chcl3",
        "cs2",
        "dioxane",
        "dmf",
        "dmso",
        "ethanol",
        "ether",
        "ethylacetate",
        "furane",
        "hexadecane",
        "hexane",
        "h2o",  # Alias for water in GBSA
        "methanol",
        "n-hexane",  # GFN2-xTB only in GBSA
        "nitromethane",
        "octanol",
        "woctanol",
        "phenol",
        "thf",
        "toluene",
        "water",
    ]

    CREST_QUICK_MODES = ["quick", "squick", "mquick"]

    @property
    def crest_methods(self):
        """
        Get CREST calculation methods.

        Returns:
            list: Available CREST methods in lowercase
        """
        return [method.lower() for method in self.CREST_METHODS]

    @property
    def crest_opt_levels(self):
        """
        Get CREST optimization convergence levels.

        Returns:
            list: Available optimization levels in lowercase
        """
        return [level.lower() for level in self.CREST_OPT_LEVELS]

    @property
    def crest_jobtypes(self):
        """
        Get CREST job types.

        Returns:
            list: Available job types in lowercase
        """
        return [jobtype.lower() for jobtype in self.CREST_JOB_TYPES]

    @property
    def crest_solvent_models(self):
        """
        Get CREST implicit solvation models.

        Returns:
            list: Available solvation models in lowercase
        """
        return [model.lower() for model in self.CREST_SOLVENT_MODELS]

    @property
    def crest_solvent_ids(self):
        """
        Get CREST solvent identifiers for implicit solvation models.

        Returns:
            list: Available solvent IDs in lowercase
        """
        return [solvent_id.lower() for solvent_id in self.CREST_SOLVENT_IDS]

    @property
    def crest_quick_modes(self):
        """
        Get CREST quick modes.

        Returns:
            list: Available quick modes in lowercase
        """
        return [quick_mode.lower() for quick_mode in self.CREST_QUICK_MODES]


# Global constants for convenient access to CREST reference data
crest_ref = CRESTRefs()
CREST_ALL_METHODS = crest_ref.crest_methods
CREST_ALL_OPT_LEVELS = crest_ref.crest_opt_levels
CREST_ALL_JOB_TYPES = crest_ref.crest_jobtypes
CREST_ALL_SOLVENT_MODELS = crest_ref.crest_solvent_models
CREST_ALL_SOLVENT_IDS = crest_ref.crest_solvent_ids
CREST_ALL_QUICK_MODES = crest_ref.crest_quick_modes
