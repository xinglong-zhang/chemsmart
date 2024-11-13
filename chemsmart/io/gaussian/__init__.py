class GaussianRefs:
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
    g_bases = ["3-", "6-", "def", "def2", "lan", "cc", "aug", "gen"]
    g_solvation_models = [
        "smd",
        "cpcm",
        "iefpcm",
        "pcm",
        "scipcm",
        "ipcm",
        "dipole",
    ]
    g_additional_route_parameters = ["force", "nosymm", "guess"]
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
    g_dieze_tags = ["#n", "#p", "#t"]

    @property
    def gaussian_ab_initio(self):
        return self.g_ab_initio

    @property
    def gaussian_dft_fuctionals(self):
        return self.g_functionals

    @property
    def gaussian_basis_sets(self):
        return self.g_bases

    @property
    def gaussian_solvation_models(self):
        return self.g_solvation_models

    @property
    def gaussian_additional_route_parameters(self):
        return self.g_additional_route_parameters

    @property
    def gaussian_additional_opt_options(self):
        return self.g_additional_opt_options

    @property
    def gaussian_dieze_tags(self):
        return self.g_dieze_tags


class BSEMetadata:
    def __init__(self):
        try:
            import basis_set_exchange as bse

            self.bse = bse
        except ImportError as e:
            raise ImportError(
                "basis_set_exchange module needed.\n"
                "see https://github.com/MolSSI-BSE/basis_set_exchange for installation."
            ) from e

    def all_bases_names(self):
        all_bases = self.bse.get_all_basis_names()
        # convert all to lower case
        return [s.lower() for s in all_bases]

    def all_formats(self):
        return self.bse.get_formats()
        # all formats = {
        # 'nwchem': 'NWChem',
        # 'gaussian94': 'Gaussian',
        # 'psi4': 'Psi4',
        # 'molcas': 'Molcas',
        # 'qchem': 'Q-Chem',
        # 'orca': 'ORCA',
        # 'dalton': 'Dalton',
        # 'qcschema': 'QCSchema',
        # 'cp2k': 'CP2K',
        # 'pqs': 'PQS',
        # 'demon2k': 'deMon2K',
        # 'gamess_us': 'GAMESS US',
        # 'turbomole': 'Turbomole',
        # 'gamess_uk': 'GAMESS UK',
        # 'molpro': 'Molpro',
        # 'cfour': 'CFOUR',
        # 'acesii': 'ACES II',
        # 'xtron': 'xTron',
        # 'bsedebug': 'BSE Debug',
        # 'json': 'JSON',
        # 'bdf': 'BDF'
        # }

    def list_of_available_formats(self):
        """Available formats.

        list_of_formats = ['nwchem', 'gaussian94', 'psi4', 'molcas', 'qchem', 'orca', 'dalton', 'qcschema', 'cp2k',
        'pqs', 'demon2k', 'gamess_us', 'turbomole', 'gamess_uk', 'molpro', 'cfour', 'acesii',
        'xtron', 'bsedebug', 'json', 'bdf']
        """
        return list(self.all_formats().keys())
