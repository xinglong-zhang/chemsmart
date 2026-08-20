class GaussianRefs:
    """Reference data for Gaussian quantum chemistry calculations.

    This class contains lists of supported methods, basis sets, solvation
    models, and other parameters available in Gaussian software.
    """

    # Ab initio methods supported by Gaussian
    g_ab_initio = [
        "hf",
        "uhf",
        "rhf",
        "mp2",
        "mp4",
        "cis",
        "cisdt",
        "ccd",
        "ccsd",
        "ccsd(t)",
    ]

    # Density functional theory (DFT) functionals
    g_functionals = [
        "pbe",
        "tpss",
        "vwn",
        "lyp",
        "p86",
        "b86",
        "b97",
        "pw91",
        "m05",
        "m06",
        "mn15",
    ]

    # Basis set families and prefixes
    g_bases = [
        "3-",
        "6-",
        "def",
        "def2",
        "lan",
        "cc",
        "aug",
        "gen",
        "genecp",
        "sv",
        "tzv",
    ]

    # Semi-empirical methods
    g_semiempirical = [
        "AM1",
        "PM3",
        "PM6",
        "PM7",
        "CNDO",
        "INDO",
        "MNDO",
        "MINDO3",
    ]

    # Complete standalone functional keywords accepted on a Gaussian 16
    # route line. Curated from the Gaussian 16 Rev C.01 keyword
    # documentation ("Density Functional (DFT) Methods"); provenance is
    # manual -- the Gaussian job-submission hold on this host forbids the
    # binary probe the ORCA vocabulary received -- so this list holds only
    # keywords the documentation names explicitly. It is a vocabulary of
    # standalone route tokens, unlike g_functionals above, which is a
    # substring-fragment list used for route classification and must not be
    # mistaken for one.
    g_all_functionals = [
        # local / pure
        "lsda",
        "svwn",
        "svwn5",
        "blyp",
        "bp86",
        "bpw91",
        "bpbe",
        "bvwn",
        "bvp86",
        "pbepbe",
        "hcth",
        "hcth93",
        "hcth147",
        "hcth407",
        "thcth",
        "tpsstpss",
        "revtpssrevtpss",
        "m06l",
        "m11l",
        "mn12l",
        "mn15l",
        "n12",
        "sogga11",
        "b97d",
        "b97d3",
        "vsxc",
        # global hybrids
        "b3lyp",
        "b3p86",
        "b3pw91",
        "o3lyp",
        "x3lyp",
        "b1b95",
        "b1lyp",
        "b98",
        "b971",
        "b972",
        "bmk",
        "bhandh",
        "bhandhlyp",
        "mpw1pw91",
        "mpw1lyp",
        "mpw1pbe",
        "mpw3pbe",
        "pbe1pbe",
        "pbeh1pbe",
        "hseh1pbe",
        "ohse2pbe",
        "ohse1pbe",
        "tpssh",
        "thcthhyb",
        "apf",
        "apfd",
        "m05",
        "m052x",
        "m06",
        "m062x",
        "m06hf",
        "m08hx",
        "mn15",
        "pw6b95",
        "pw6b95d3",
        "sogga11x",
        "hissbpbe",
        "x3lyp",
        # range-separated
        "cam-b3lyp",
        "lc-whpbe",
        "wb97",
        "wb97x",
        "wb97xd",
        "m11",
        "mn12sx",
        "n12sx",
        # double hybrids
        "b2plyp",
        "b2plypd3",
        "mpw2plyp",
    ]

    # Solvation models available in Gaussian
    g_solvation_models = [
        "smd",
        "cpcm",
        "iefpcm",
        "pcm",
        "scipcm",
        "ipcm",
        "dipole",
    ]

    # Additional route line parameters
    g_additional_route_parameters = [
        "force",
        "nosymm",
        "guess",
        "stable",
        "scf",
        "geom",  # Used for geom=check keyword when reading checkpoint files
    ]

    # Additional optimization options
    g_additional_opt_options = [
        "maxstep",
        "maxcycles",
        "restart",
        "calcall",
        "readfc",
        "recalcfc",
        "vcd",
        "noraman",
        "tight",
        "verytight",
        "loose",
        "expert",
        "z-matrix",
        "cartesian",
        "gic",
        "micro",
        "quadmacro",
    ]

    # Gaussian route line verbosity tags
    g_dieze_tags = ["#n", "#p", "#t"]

    @property
    def gaussian_ab_initio(self):
        """Get list of ab initio methods."""
        return self.g_ab_initio

    @property
    def gaussian_dft_fuctionals(self):
        """Get list of DFT functionals."""
        return self.g_functionals

    @property
    def gaussian_basis_sets(self):
        """Get list of basis set families."""
        return self.g_bases

    @property
    def gaussian_semiempirical_methods(self):
        """Get list of semi-empirical methods."""
        return self.g_semiempirical

    @property
    def gaussian_all_functionals(self):
        """Standalone G16 functional keywords, lower-case and deduplicated."""

        return sorted(
            {name.strip().lower() for name in self.g_all_functionals}
        )

    @property
    def gaussian_solvation_models(self):
        """Get list of solvation models."""
        return self.g_solvation_models

    @property
    def gaussian_additional_route_parameters(self):
        """Get list of additional route parameters."""
        return self.g_additional_route_parameters

    @property
    def gaussian_additional_opt_options(self):
        """Get list of additional optimization options."""
        return self.g_additional_opt_options

    @property
    def gaussian_dieze_tags(self):
        """Get list of route line verbosity tags."""
        return self.g_dieze_tags


class BSEMetadata:
    """
    Basis Set Exchange metadata and interface.
    """

    def __init__(self):
        try:
            import basis_set_exchange as bse

            self.bse = bse
        except ImportError as e:
            raise ImportError(
                "basis_set_exchange module needed.\n"
                "See https://github.com/MolSSI-BSE/basis_set_exchange "
                "for installation."
            ) from e

    def all_bases_names(self):
        """Get all available basis set names in lowercase.

        Returns:
            list: List of all basis set names converted to lowercase
        """
        all_bases = self.bse.get_all_basis_names()
        # Convert all to lowercase for consistency
        return [s.lower() for s in all_bases]

    def all_formats(self):
        """Get all available output formats.

        Returns:
            dict: Dictionary mapping format keys to format names
        """
        return self.bse.get_formats()
        # Available formats include:
        # {
        #     'nwchem': 'NWChem',
        #     'gaussian94': 'Gaussian',
        #     'psi4': 'Psi4',
        #     'molcas': 'Molcas',
        #     'qchem': 'Q-Chem',
        #     'orca': 'ORCA',
        #     'dalton': 'Dalton',
        #     'qcschema': 'QCSchema',
        #     'cp2k': 'CP2K',
        #     'pqs': 'PQS',
        #     'demon2k': 'deMon2K',
        #     'gamess_us': 'GAMESS US',
        #     'turbomole': 'Turbomole',
        #     'gamess_uk': 'GAMESS UK',
        #     'molpro': 'Molpro',
        #     'cfour': 'CFOUR',
        #     'acesii': 'ACES II',
        #     'xtron': 'xTron',
        #     'bsedebug': 'BSE Debug',
        #     'json': 'JSON',
        #     'bdf': 'BDF'
        # }

    def list_of_available_formats(self):
        """Get list of available output format keys.

        Available formats include: nwchem, gaussian94, psi4, molcas, qchem,
        orca, dalton, qcschema, cp2k, pqs, demon2k, gamess_us, turbomole,
        gamess_uk, molpro, cfour, acesii, xtron, bsedebug, json, bdf.

        Returns:
            list: List of format keys
        """
        return list(self.all_formats().keys())


# Initialize Gaussian reference data
gaussian_ref = GaussianRefs()

# Export commonly used reference lists for easy access
GAUSSIAN_AB_INITIO = gaussian_ref.gaussian_ab_initio
GAUSSIAN_FUNCTIONALS = gaussian_ref.gaussian_dft_fuctionals
GAUSSIAN_BASES = gaussian_ref.gaussian_basis_sets
GAUSSIAN_SEMIEMPIRICAL = gaussian_ref.gaussian_semiempirical_methods
GAUSSIAN_SOLVATION_MODELS = gaussian_ref.gaussian_solvation_models
GAUSSIAN_ALL_FUNCTIONALS = gaussian_ref.gaussian_all_functionals
GAUSSIAN_ADDITIONAL_OPT_OPTIONS = gaussian_ref.gaussian_additional_opt_options
GAUSSIAN_DIEZE_TAGS = gaussian_ref.gaussian_dieze_tags
GAUSSIAN_ADDITIONAL_ROUTE_PARAMETERS = (
    gaussian_ref.gaussian_additional_route_parameters
)
