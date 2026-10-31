import functools


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
    g_bases = ["3-", "6-", "def", "def2", "lan", "cc", "aug", "gen", "genecp"]

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

    def basis_names_matching(self, element, shells):
        """Find the basis sets whose functions for an element are identical
        to the given shells (inverse of retrieving a basis set by name).

        The BSE basis sets are laid out the way their Gaussian (gaussian94)
        format is written, i.e. as the shells Gaussian reads from a gen/genecp
        section, and compared independently of shell order.

        Args:
            element (str): Element symbol, e.g. 'Ag'.
            shells (list): Shell dicts as given by
                Gaussian16Output.heavy_elements_basis_functions[element].

        Returns:
            list: Sorted lowercase BSE basis set names; empty if none match.
        """
        return list(
            _bse_basis_names_matching(element, _basis_shells_key(shells))
        )


def _significant(value, digits=8):
    """Round to significant digits so that values printed by Gaussian
    compare equal to the BSE values they were read from."""
    return float(f"{value:.{digits - 1}e}")


def _basis_shells_key(shells):
    """Order-independent, hashable key for a list of shell dicts."""
    return tuple(
        sorted(
            (
                shell["shell"],
                tuple(
                    tuple(_significant(value) for value in primitive)
                    for primitive in shell["primitives"]
                ),
            )
            for shell in shells
        )
    )


@functools.lru_cache(maxsize=None)
def _bse_basis_names_matching(element, shells_key):
    """Scan all BSE basis sets containing the element for ones whose
    shells give `shells_key`.  Cached because a scan takes seconds."""
    import basis_set_exchange as bse

    z = str(bse.lut.element_Z_from_sym(element))
    matches = []
    for name, metadata in bse.get_metadata().items():
        latest = metadata["versions"][metadata["latest_version"]]
        if z not in latest["elements"]:
            continue
        basis = bse.get_basis(name, elements=[element])
        # Same shell layout as BSE's gaussian94 writer
        basis = bse.manip.uncontract_general(basis, False)
        basis = bse.manip.uncontract_spdf(basis, 1, False)
        bse_shells = []
        for shell in basis["elements"][z].get("electron_shells", []):
            label = bse.lut.amint_to_char(
                shell["angular_momentum"], hij=True
            ).upper()
            primitives = [
                (
                    float(exponent),
                    *(float(c[i]) for c in shell["coefficients"]),
                )
                for i, exponent in enumerate(shell["exponents"])
            ]
            bse_shells.append({"shell": label, "primitives": primitives})
        if bse_shells and _basis_shells_key(bse_shells) == shells_key:
            # Display name, as accepted by bse.get_basis and listed by
            # all_bases_names (metadata keys are mangled, e.g. 6-31g_st_)
            matches.append(metadata["display_name"].lower())
    return tuple(sorted(matches))


# Initialize Gaussian reference data
gaussian_ref = GaussianRefs()

# Export commonly used reference lists for easy access
GAUSSIAN_AB_INITIO = gaussian_ref.gaussian_ab_initio
GAUSSIAN_FUNCTIONALS = gaussian_ref.gaussian_dft_fuctionals
GAUSSIAN_BASES = gaussian_ref.gaussian_basis_sets
GAUSSIAN_SEMIEMPIRICAL = gaussian_ref.gaussian_semiempirical_methods
GAUSSIAN_SOLVATION_MODELS = gaussian_ref.gaussian_solvation_models
GAUSSIAN_ADDITIONAL_OPT_OPTIONS = gaussian_ref.gaussian_additional_opt_options
GAUSSIAN_DIEZE_TAGS = gaussian_ref.gaussian_dieze_tags
GAUSSIAN_ADDITIONAL_ROUTE_PARAMETERS = (
    gaussian_ref.gaussian_additional_route_parameters
)
