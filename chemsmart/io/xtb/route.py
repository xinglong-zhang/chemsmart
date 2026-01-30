import logging

from chemsmart.io.xtb import (
    XTB_ALL_METHODS,
    XTB_ALL_OPT_LEVELS,
    XTB_ALL_SOLVENT_IDS,
    XTB_ALL_SOLVENT_MODELS,
)

logger = logging.getLogger(__name__)


class XTBRoute:
    """
    Parser for xTB route (program call) specifications.

    This class parses xTB command line strings to extract job types and
    calculation parameters from the program call, enabling programmatic
    analysis of xTB calculation settings.

    The parser recognizes various xTB command line flags including:
    - GFN methods (GFN0, GFN1, GFN2, GFN-FF)
    - Job types (single point, optimization, Hessian, MD, etc.)
    - Molecular properties (charge, unpaired electrons)
    - Optimization levels (crude to extreme)
    - Solvation models (ALPB, GBSA, COSMO, etc.)
    - Property calculations (dipole, Mulliken, Wiberg, etc.)

    Args:
        route_string (str): xTB program call string
    """

    def __init__(self, route_string):
        """
        Initialize xTB route parser.

        Args:
            route_string (str): xTB program call string to parse (e.g., "xtb coord --opt")
        """
        self.route_string = route_string.lower()
        self.route_inputs = self.route_string.split()

    @property
    def gfn_version(self):
        """
        Extract GFN version from the route.

        xTB GFN (Geometry, Frequency, Non-covalent interactions) versions:
        - 'gfn0': GFN0-xTB - Minimal basis parametrization
        - 'gfn1': GFN1-xTB - GFN-xTB version 1 (faster, lower accuracy)
        - 'gfn2': GFN2-xTB - GFN-xTB version 2 (default, balanced accuracy/speed)
        - 'gfnff': GFN-FF - GFN force field (fastest, suitable for large systems)

        Command line flags: --gfn 0|1|2, --gfnff (or --gff)

        Returns:
            str: GFN version identifier ('gfn0', 'gfn1', 'gfn2', 'gfnff') or None
        """
        # Check for explicit --gfnX flags
        for method in XTB_ALL_METHODS:
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

        xTB optimization levels (ancopt convergence criteria):
        - 'crude'     : Very rough optimization, fastest, lowest accuracy
        - 'sloppy'    : Sloppy optimization for quick pre-relaxation
        - 'loose'     : Loose convergence criteria
        - 'lax'       : Intermediate between loose and normal
        - 'normal'    : Default optimization level in xTB
        - 'tight'     : Tight convergence criteria
        - 'vtight'    : Very tight convergence criteria (verytight in CLI)
        - 'extreme'   : Extremely tight convergence, highest accuracy

        Command line flags: --opt [LEVEL], -o [LEVEL], --optlevel [LEVEL], --metaopt [LEVEL]

        Returns:
            str: Optimization level identifier or None if not specified
        """
        for flag in ("--opt", "-o", "--optlevel", "--metaopt"):
            if flag in self.route_inputs:
                opt_index = self.route_inputs.index(flag)
                if opt_index + 1 < len(self.route_inputs):
                    opt_level = self.route_inputs[opt_index + 1]
                    # Normalize 'verytight' to 'vtight'
                    if opt_level == "verytight":
                        return "vtight"
                    if opt_level in XTB_ALL_OPT_LEVELS:
                        return opt_level
        return None

    @property
    def charge(self):
        """
        Extract molecular charge specification from the route.

        Command line flags: --chrg INT, -c INT

        This overrides the .CHRG file and xcontrol option.

        Returns:
            int: Molecular charge value or None if not specified
        """
        for flag in ("--chrg", "-c"):
            if flag in self.route_inputs:
                chrg_index = self.route_inputs.index(flag)
                if chrg_index + 1 < len(self.route_inputs):
                    try:
                        return int(self.route_inputs[chrg_index + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def uhf(self):
        """
        Extract unpaired electron (Nalpha - Nbeta) specification from the route.

        Command line flags: --uhf INT, -u INT

        This specifies the number of unpaired electrons and overrides the .UHF file
        and xcontrol option.

        Returns:
            int: Number of unpaired electrons or None if not specified
        """
        for flag in ("--uhf", "-u"):
            if flag in self.route_inputs:
                uhf_index = self.route_inputs.index(flag)
                if uhf_index + 1 < len(self.route_inputs):
                    try:
                        return int(self.route_inputs[uhf_index + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def jobtype(self):
        """
        Extract the primary job type from the route.

        xTB job types:
        - 'sp'    : Single point calculation (--scc, --sp, or default)
        - 'opt'   : Geometry optimization (--opt, --omd)
        - 'hess'  : Hessian/frequency calculation (--hess, --ohess, --bhess)
        - 'md'    : Molecular dynamics (--md, --metadyn)
        - 'path'  : Reaction path calculation (--path)
        - 'modef' : Mode following (--modef)
        """
        return self.get_jobtype()

    def get_jobtype(self):
        """
        Extract job type from route specification.

        Checks for various job type keywords in order of precedence:
        1. Optimization: --opt, --omd, --metaopt, --ohess, --bhess
        2. Hessian: --hess
        3. MD: --md, --metadyn
        4. Path: --path
        5. Mode following: --modef
        6. Default: sp (single point)

        Returns:
            str: Job type ('opt', 'hess', 'md', 'path', 'modef', or 'sp')
        """
        # Check for optimization keywords (including ohess and bhess which do opt first)
        if any(
            flag in self.route_inputs
            for flag in (
                "--opt",
                "-o",
                "--omd",
                "--metaopt",
                "--ohess",
                "--bhess",
            )
        ):
            return "opt"
        # Check for Hessian calculation
        if "--hess" in self.route_inputs:
            return "hess"
        # Check for molecular dynamics
        if "--md" in self.route_inputs or "--metadyn" in self.route_inputs:
            return "md"
        # Check for path calculation
        if "--path" in self.route_inputs:
            return "path"
        # Check for mode following
        if "--modef" in self.route_inputs:
            return "modef"
        # Default is single point
        return "sp"

    @property
    def freq(self):
        """
        Check if frequency/Hessian calculation is requested.

        xTB frequency keywords: --hess, --ohess, --bhess

        - --hess  : Numerical Hessian on input geometry
        - --ohess : Numerical Hessian on optimized geometry
        - --bhess : Biased numerical Hessian on optimized geometry

        Returns:
            bool: True if frequency calculation is requested
        """
        return any(
            flag in self.route_inputs
            for flag in ("--hess", "--ohess", "--bhess")
        )

    @property
    def grad(self):
        """
        Check if gradient calculation is requested.

        Command line flag: --grad

        Returns:
            bool: True if gradient calculation is requested
        """
        return "--grad" in self.route_inputs

    @property
    def solvent_model(self):
        """
        Extract solvation model specification from the route.

        Available solvation models:
        - 'alpb'    : Analytical Linearized Poisson-Boltzmann (ALPB) model
        - 'gbsa'    : Generalized Born (GB) model with Solvent Accessible Surface Area (SASA)
        - 'cosmo'   : Domain Decomposition Conductor-like Screening Model (ddCOSMO)
        - 'tmcosmo' : COSMO with TM (Turbomole) convention
        - 'cpcmx'   : Extended Conduction-like Polarizable Continuum Model (CPCM-X)

        Command line flags: --alpb, --gbsa (or -g), --cosmo, --tmcosmo, --cpcmx

        Returns:
            str: Solvent model identifier or None if not specified
        """
        return self.get_solvent_model()

    def get_solvent_model(self):
        """
        Extract solvent model from route specification.

        Returns:
            str: Solvent model identifier ('alpb', 'gbsa', 'cosmo', 'tmcosmo', 'cpcmx')
                 or None if not specified
        """
        for model in XTB_ALL_SOLVENT_MODELS:
            if f"--{model}" in self.route_inputs:
                return model
        # Short flag for GBSA
        if "-g" in self.route_inputs:
            return "gbsa"
        return None

    @property
    def solvent_id(self):
        """
        Extract solvent identity specification from the route.

        Available solvents (case-insensitive):
        - ALPB: acetone, acetonitrile, aniline, benzaldehyde, benzene, ch2cl2, chcl3,
                cs2, dioxane, dmf, dmso, ether, ethylacetate, furane, hexandecane, hexane,
                methanol, nitromethane, octanol, woctanol, phenol, toluene, thf, water
        - GBSA: acetone, acetonitrile, benzene (GFN1 only), ch2cl2, chcl3, cs2,
                dmf (GFN2 only), dmso, ether, h2o/water, methanol, n-hexane (GFN2 only),
                thf, toluene
        - COSMO: All ALPB solvents, or dielectric constant (EPSILON)
        - CPCMX: All solvents in Minnesota Solvation Database

        Returns:
            str: Solvent identity (e.g., 'water', 'ethanol') or None if not specified
        """
        if self.solvent_model is None:
            return None
        for solvent in self.route_inputs:
            if solvent in XTB_ALL_SOLVENT_IDS:
                return solvent
        return None

    @property
    def accuracy(self):
        """
        Extract accuracy for SCC (Self-Consistent Charge) calculation.

        Command line flags: --acc REAL, -a REAL

        Lower values mean higher accuracy (default = 1.0).

        Returns:
            float: Accuracy value or None if not specified
        """
        for flag in ("--acc", "-a"):
            if flag in self.route_inputs:
                acc_index = self.route_inputs.index(flag)
                if acc_index + 1 < len(self.route_inputs):
                    try:
                        return float(self.route_inputs[acc_index + 1])
                    except (ValueError, IndexError):
                        pass
        return None

    @property
    def electronic_temperature(self):
        """
        Extract electronic temperature specification.

        Command line flag: --etemp REAL

        Temperature in Kelvin (default = 300K).

        Returns:
            float: Electronic temperature in Kelvin or None if not specified
        """
        if "--etemp" in self.route_inputs:
            etemp_index = self.route_inputs.index("--etemp")
            if etemp_index + 1 < len(self.route_inputs):
                try:
                    return float(self.route_inputs[etemp_index + 1])
                except (ValueError, IndexError):
                    pass
        return None

    @property
    def ptb(self):
        """
        Check if PTB (Density Tight-Binding) single-point calculation is requested.

        Command line flag: --ptb

        PTB provides electronic structure and properties (charges, bond orders, dipole)
        but not energy-related properties (total energy, gradients, frequencies).

        Returns:
            bool: True if PTB calculation is requested
        """
        return "--ptb" in self.route_inputs

    @property
    def spinpol(self):
        """
        Check if spin-polarization is enabled for xTB methods.

        Command line flag: --spinpol

        Note: Requires tblite library.

        Returns:
            bool: True if spin-polarization is enabled
        """
        return "--spinpol" in self.route_inputs

    @property
    def ceh(self):
        """
        Check if CEH (Charge-Extended Hückel) charges calculation is requested.

        Command line flag: --ceh

        Writes charges to ceh.charges file.

        Returns:
            bool: True if CEH calculation is requested
        """
        return "--ceh" in self.route_inputs

    @property
    def pop(self):
        """
        Check if Mulliken population analysis printout is requested.

        Command line flag: --pop

        Returns:
            bool: True if population analysis is requested
        """
        return "--pop" in self.route_inputs

    @property
    def wbo(self):
        """
        Check if Wiberg bond order printout is requested.

        Command line flag: --wbo

        Returns:
            bool: True if Wiberg bond order calculation is requested
        """
        return "--wbo" in self.route_inputs

    @property
    def dipole(self):
        """
        Check if dipole moment printout is requested.

        Command line flag: --dipole

        Returns:
            bool: True if dipole calculation is requested
        """
        return "--dipole" in self.route_inputs

    @property
    def molden(self):
        """
        Check if Molden file output is requested.

        Command line flag: --molden

        Returns:
            bool: True if Molden file output is requested
        """
        return "--molden" in self.route_inputs

    @property
    def lmo(self):
        """
        Check if localized molecular orbitals (LMO) calculation is requested.

        Command line flag: --lmo

        Returns:
            bool: True if LMO calculation is requested
        """
        return "--lmo" in self.route_inputs

    @property
    def fod(self):
        """
        Check if FOD (Fractional Occupation Density) calculation is requested.

        Command line flag: --fod

        Returns:
            bool: True if FOD calculation is requested
        """
        return "--fod" in self.route_inputs

    @property
    def esp(self):
        """
        Check if electrostatic potential on VdW grid calculation is requested.

        Command line flag: --esp

        Returns:
            bool: True if ESP calculation is requested
        """
        return "--esp" in self.route_inputs

    @property
    def stm(self):
        """
        Check if STM (Scanning Tunneling Microscopy) image calculation is requested.

        Command line flag: --stm

        Returns:
            bool: True if STM calculation is requested
        """
        return "--stm" in self.route_inputs

    @property
    def vip(self):
        """
        Check if vertical ionization potential calculation is requested.

        Command line flag: --vip

        Requires param_ipea-xtb.txt parameters and GFN1 Hamiltonian.

        Returns:
            bool: True if VIP calculation is requested
        """
        return "--vip" in self.route_inputs

    @property
    def vea(self):
        """
        Check if vertical electron affinity calculation is requested.

        Command line flag: --vea

        Requires param_ipea-xtb.txt parameters and GFN1 Hamiltonian.

        Returns:
            bool: True if VEA calculation is requested
        """
        return "--vea" in self.route_inputs

    @property
    def vipea(self):
        """
        Check if both VIP and VEA calculations are requested.

        Command line flag: --vipea

        Requires param_ipea-xtb.txt parameters and GFN1 Hamiltonian.

        Returns:
            bool: True if both VIP and VEA calculations are requested
        """
        return "--vipea" in self.route_inputs

    @property
    def vfukui(self):
        """
        Check if Fukui indices calculation is requested.

        Command line flag: --vfukui

        Returns:
            bool: True if Fukui indices calculation is requested
        """
        return "--vfukui" in self.route_inputs

    @property
    def vomega(self):
        """
        Check if electrophilicity index calculation is requested.

        Command line flag: --vomega

        Requires param_ipea-xtb.txt parameters and GFN1 Hamiltonian.

        Returns:
            bool: True if electrophilicity index calculation is requested
        """
        return "--vomega" in self.route_inputs

    @property
    def alpha(self):
        """
        Check if static molecular dipole polarizability calculation is requested.

        Command line flag: --alpha

        Returns:
            bool: True if polarizability calculation is requested
        """
        return "--alpha" in self.route_inputs

    @property
    def cma(self):
        """
        Check if center of mass shift and principal axis transformation is requested.

        Command line flag: --cma

        Shifts molecule to center of mass and transforms Cartesian coordinates
        into the coordinate system of the principal axis.

        Returns:
            bool: True if CMA transformation is requested
        """
        return "--cma" in self.route_inputs
