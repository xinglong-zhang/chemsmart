import logging

from chemsmart.io.crest import (
    CREST_ALL_METHODS,
    CREST_ALL_OPT_LEVELS,
    CREST_ALL_SOLVENT_IDS,
    CREST_ALL_SOLVENT_MODELS,
)

logger = logging.getLogger(__name__)


class CRESTRoute:
    """
    Parser for CREST route (program call) specifications.

    This class parses CREST command line strings to extract job types and
    calculation parameters from the program call, enabling programmatic
    analysis of CREST calculation settings.

    Args:
        route_string (str): CREST command line string
    """

    def __init__(self, route_string):
        """
        Initialize CREST route parser.

        Args:
            route_string (str): CREST program call string to parse
                (e.g., "crest mol.xyz --gfn2 --chrg 0")
        """
        self.route_string = route_string.lower()
        self.route_inputs = self.route_string.split()

    @property
    def method(self):
        """Extract the computational method (GFN method)."""
        return self.gfn_version

    @property
    def basis(self):
        """
        Return default xTB basis set specification.

        CREST uses xTB Hamiltonians with a built-in default basis
        and does not accept user-specified basis sets.
        """
        return "default"

    @property
    def gfn_version(self):
        """
        Extract GFN version from the route.

        CREST supported xTB GFN versions:
        - 'gfn1': GFN1-xTB
        - 'gfn2': GFN2-xTB (default)
        - 'gfnff': GFN-FF
        - 'gfn2//gfnff': GFN2 energetics with GFN-FF sampling

        Command line flags: --gfn0, --gfn1, --gfn2, --gfnff (or --gff), --gfn2//gfnff,
        --gfn INT

        Returns:
            str: GFN version identifier or None
        """
        for method in CREST_ALL_METHODS:
            if f"--{method}" in self.route_inputs:
                return method

        # Check for --gfnff alias --gff
        if "--gff" in self.route_inputs:
            return "gfnff"

        # Check for --gfn INT format
        if "--gfn" in self.route_inputs:
            gfn_index = self.route_inputs.index("--gfn") + 1
            if gfn_index < len(self.route_inputs):
                try:
                    gfn_num = int(self.route_inputs[gfn_index])
                    if gfn_num == 0:
                        return "gfn0"
                    elif gfn_num == 1:
                        return "gfn1"
                    elif gfn_num == 2:
                        return "gfn2"
                except (ValueError, IndexError):
                    pass

        return None

    @property
    def optimization_level(self):
        """
        Extract optimization level from the route.

        Command line flags: --optlev LEVEL, -optlev LEVEL,
        --opt [LEVEL], -o [LEVEL], --optlevel [LEVEL], --ancopt [LEVEL]

        Returns:
            str: Optimization level identifier or None if not specified
        """
        for flag in (
            "--optlev",
            "-optlev",
            "--opt",
            "-o",
            "--optlevel",
            "--ancopt",
        ):
            if flag in self.route_inputs:
                opt_index = self.route_inputs.index(flag)
                if opt_index + 1 < len(self.route_inputs):
                    opt_level = self.route_inputs[opt_index + 1]
                    if opt_level == "verytight":
                        return "vtight"
                    if opt_level in CREST_ALL_OPT_LEVELS:
                        return opt_level
        return None

    @property
    def charge(self):
        """
        Extract molecular charge specification from the route.

        Command line flags: --chrg INT, -chrg INT

        Returns:
            int: Molecular charge value or 0 if not specified
        """
        for flag in ("--chrg", "-chrg"):
            if flag in self.route_inputs:
                chrg_index = self.route_inputs.index(flag)
                if chrg_index + 1 < len(self.route_inputs):
                    try:
                        return int(self.route_inputs[chrg_index + 1])
                    except (ValueError, IndexError):
                        pass
        return 0

    @property
    def uhf(self):
        """
        Extract unpaired electron (Nalpha - Nbeta) specification from the route.

        Command line flags: --uhf INT, -uhf INT

        Returns:
            int: Number of unpaired electrons or 0 if not specified
        """
        for flag in ("--uhf", "-uhf"):
            if flag in self.route_inputs:
                uhf_index = self.route_inputs.index(flag)
                if uhf_index + 1 < len(self.route_inputs):
                    try:
                        return int(self.route_inputs[uhf_index + 1])
                    except (ValueError, IndexError):
                        pass
        return 0

    @property
    def multiplicity(self):
        """Spin multiplicity (uhf + 1)."""
        return self.uhf + 1

    @property
    def jobtype(self):
        """Extract the primary job type from the route."""
        return self.get_jobtype()

    def get_jobtype(self):
        """
        Extract job type from route specification.

        CREST's primary purpose is conformer–rotamer ensemble sampling,
        so the job type is always conformational search.

        Returns:
            str: Job type ('conformers').
        """
        return "conformers"

    @property
    def solvent_model(self):
        """
        Extract solvation model specification from the route.

        Available solvation models:
        - 'alpb': Analytical Linearized Poisson-Boltzmann (ALPB) model
        - 'gbsa': Generalized Born (GB) model with SASA

        Command line flags: --alpb, --gbsa (or -g)

        Returns:
            str: Solvent model identifier or None if not specified
        """
        return self.get_solvent_model()

    def get_solvent_model(self):
        """
        Extract solvent model from route specification.

        Returns:
            str: Solvent model identifier ('alpb', 'gbsa') or None
        """
        for model in CREST_ALL_SOLVENT_MODELS:
            if f"--{model}" in self.route_inputs:
                return model
        if "-g" in self.route_inputs:
            return "gbsa"
        return None

    @property
    def solvent_id(self):
        """
        Extract solvent identity specification from the route.

        Returns:
            str: Solvent identity (e.g., 'water', 'toluene') or None
        """
        if self.solvent_model is None:
            return None
        for solvent in self.route_inputs:
            if solvent in CREST_ALL_SOLVENT_IDS:
                return solvent
        return None

    @property
    def energy_window(self):
        """
        Extract energy window in kcal/mol from the route.

        Command line flags: --ewin REAL, -ewin REAL

        Returns:
            float: Energy window or None if not specified
        """
        for flag in ("--ewin", "-ewin"):
            if flag in self.route_inputs:
                ewin_index = self.route_inputs.index(flag)
                if ewin_index + 1 < len(self.route_inputs):
                    try:
                        return float(self.route_inputs[ewin_index + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def nci(self):
        """
        Check if non-covalent interaction mode is requested.

        Command line flags: --nci, -nci

        Returns:
            bool: True if --nci mode was used
        """
        return "--nci" in self.route_inputs or "-nci" in self.route_inputs

    @property
    def constrained(self):
        """
        Check if a constraint input file was specified.

        Command line flags: --cinp, -cinp

        Returns:
            bool: True if --cinp was specified
        """
        return "--cinp" in self.route_inputs or "-cinp" in self.route_inputs

    @property
    def no_reference_topology_check(self):
        """
        Check if reference topology check is disabled.

        Command line flags: --noreftopo, -noreftopo

        Returns:
            bool: True if --noreftopo was used
        """
        return (
            "--noreftopo" in self.route_inputs
            or "-noreftopo" in self.route_inputs
        )

    @property
    def threads(self):
        """
        Extract number of threads from the route.

        Command line flags: -T INT, --T INT (stored lowercase as -t / --t)

        Returns:
            int: Number of threads or None if not specified
        """
        for flag in ("-t", "--t"):
            if flag in self.route_inputs:
                t_index = self.route_inputs.index(flag)
                if t_index + 1 < len(self.route_inputs):
                    try:
                        return int(self.route_inputs[t_index + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def quick_mode(self):
        """
        Extract quick mode from the route.

        Command line flags: --quick, --squick, --mquick

        Returns:
            str: Quick mode identifier or None if not specified
        """
        for mode in ("mquick", "squick", "quick"):
            if (
                f"--{mode}" in self.route_inputs
                or f"-{mode}" in self.route_inputs
            ):
                return mode
        return None

    @property
    def rmsd_threshold(self):
        """
        Extract RMSD threshold in Angstrom from the route.

        Command line flags: --rthr REAL, -rthr REAL
        Default in CREST: 0.125 Ang.

        Returns:
            float: RMSD threshold or None if not specified
        """
        for flag in ("--rthr", "-rthr"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return float(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def energy_threshold(self):
        """
        Extract energy threshold in kcal/mol from the route.

        Command line flags: --ethr REAL, -ethr REAL
        Default in CREST: 0.05 kcal/mol.

        Returns:
            float: Energy threshold or None if not specified
        """
        for flag in ("--ethr", "-ethr"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return float(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def bconst_threshold(self):
        """
        Extract rotational constant threshold from the route.

        Command line flags: --bthr REAL, -bthr REAL
        Default in CREST: 0.01 (= 1%).

        Returns:
            float: Rotational constant threshold or None if not specified
        """
        for flag in ("--bthr", "-bthr"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return float(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def population_threshold(self):
        """
        Extract Boltzmann population threshold from the route.

        Command line flags: --pthr REAL, -pthr REAL
        Default in CREST: 0.05 (= 5%).

        Returns:
            float: Population threshold or None if not specified
        """
        for flag in ("--pthr", "-pthr"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return float(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def temperature(self):
        """
        Extract CREGEN temperature in Kelvin from the route.

        Command line flags: --temp REAL, -temp REAL
        Default in CREST: 298.15 K.

        Returns:
            float: Temperature in Kelvin or None if not specified
        """
        for flag in ("--temp", "-temp"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return float(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def no_topology_check(self):
        """
        Check if topology checks in CREGEN are disabled.

        Command line flags: --notopo, -notopo

        Returns:
            bool: True if --notopo was used
        """
        return (
            "--notopo" in self.route_inputs or "-notopo" in self.route_inputs
        )

    @property
    def shake(self):
        """
        Extract SHAKE mode for MD from the route.

        Command line flags: --shake INT, -shake INT
        Values: 0=off, 1=H-only, 2=all bonds. Default in CREST: 2.

        Returns:
            int: SHAKE mode or None if not specified
        """
        for flag in ("--shake", "-shake"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return int(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def md_timestep(self):
        """
        Extract MD time step in fs from the route.

        Command line flags: --tstep INT, -tstep INT
        Default in CREST: 5 fs.

        Returns:
            int: MD time step in fs or None if not specified
        """
        for flag in ("--tstep", "-tstep"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return int(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def md_length(self):
        """
        Extract MD length from the route.

        Command line flags: --mdlen, -mdlen, --len, -len

        Value may be a length in ps (e.g. 50) or a multiplicative
        factor for the default length (e.g. x2.0).

        Returns:
            float or str: MD length in ps, multiplicative factor string
                (e.g. 'x2.0'), or None if not specified
        """
        for flag in ("--mdlen", "-mdlen", "--len", "-len"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    value = self.route_inputs[idx + 1]
                    if value.startswith("x"):
                        return value
                    try:
                        return float(value)
                    except (ValueError, IndexError):
                        return value
        return None

    @property
    def md_dump_step(self):
        """
        Extract xyz dump step to trajectory in fs from the route.

        Command line flags: --mddump INT, -mddump INT
        Default in CREST: 100 fs.

        Returns:
            int: Dump step in fs or None if not specified
        """
        for flag in ("--mddump", "-mddump"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return int(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def vbias_dump_interval(self):
        """
        Extract Vbias dump frequency in ps from the route.

        Command line flags: --vbdump REAL, -vbdump REAL
        Default in CREST: 1.0 ps.

        Returns:
            float: Vbias dump frequency in ps or None if not specified
        """
        for flag in ("--vbdump", "-vbdump"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return float(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def additional_md_temperature(self):
        """
        Extract temperature for additional normal MDs from the route.

        Command line flags: --tnmd REAL, -tnmd REAL
        Default in CREST: 400 K.

        Returns:
            float: Temperature in Kelvin or None if not specified
        """
        for flag in ("--tnmd", "-tnmd"):
            if flag in self.route_inputs:
                idx = self.route_inputs.index(flag)
                if idx + 1 < len(self.route_inputs):
                    try:
                        return float(self.route_inputs[idx + 1])
                    except (ValueError, IndexError):
                        pass
        return None
