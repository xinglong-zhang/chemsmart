import logging

from chemsmart.io.orca import (
    ORCA_ALL_AB_INITIO,
    ORCA_ALL_AUXILIARY_BASIS_SETS,
    ORCA_ALL_BASIS_SETS,
    ORCA_ALL_DISPERSION_CORRECTIONS,
    ORCA_ALL_EXTRAPOLATION_BASIS_SETS,
    ORCA_ALL_FUNCTIONALS,
    ORCA_ALL_JOB_TYPES,
    ORCA_ALL_SCF_ALGORITHMS,
    ORCA_SCF_CONVERGENCE,
)

logger = logging.getLogger(__name__)


class ORCARoute:
    """
    Parser for ORCA route (calculation method) specifications.

    This class parses ORCA route strings to extract computational methods,
    basis sets, job types, and other calculation parameters. It validates
    keywords against known ORCA options and provides convenient access
    to different route components.

    Args:
        route_string (str): ORCA route string starting with '!'
    """

    def __init__(self, route_string):
        """
        Initialize ORCA route parser.

        Args:
            route_string (str): ORCA route string to parse (e.g., '! B3LYP def2-TZVP')
        """
        self.route_string = route_string.lower()
        self.route_inputs = self.route_string.split()

    @property
    def route_keywords(self):
        """
        Extract and clean route keywords from route string.

        Returns:
            list: List of individual route keywords with '!' removed
        """
        route_keywords = []
        for raw_route_input in self.route_inputs:
            if len(raw_route_input) != 0:
                route_input = raw_route_input.replace("!", "").strip()
                if len(route_input) != 0:
                    route_keywords.append(route_input)
        return route_keywords

    @property
    def functional(self):
        """
        Extract DFT functional from route keywords.

        Returns:
            str: DFT functional name or None if not found
        """
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_FUNCTIONALS:
                return route_keyword
        return None

    @property
    def ab_initio(self):
        """
        Extract ab initio method from route keywords.

        Returns:
            str: Ab initio method name or None if not found
        """
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_AB_INITIO:
                return route_keyword
        return None

    @property
    def dispersion(self):
        """
        Extract dispersion correction from route keywords.

        Returns:
            str: Dispersion correction type or None if not specified
        """
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_DISPERSION_CORRECTIONS:
                return route_keyword
        return None

    @property
    def basis(self):
        """
        Extract basis set from route keywords.

        Returns:
            str: Basis set name or None if not found
        """
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_BASIS_SETS:
                return route_keyword
        return None

    @property
    def auxiliary_basis(self):
        """Extract auxiliary basis set from route keywords."""
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_AUXILIARY_BASIS_SETS:
                return route_keyword
        return None

    @property
    def extrapolation_basis(self):
        """Extract basis set extrapolation scheme from route keywords."""
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_EXTRAPOLATION_BASIS_SETS:
                return route_keyword
        return None

    @property
    def defgrid(self):
        """Extract integration grid specification from route keywords."""
        for route_keyword in self.route_keywords:
            if route_keyword.startswith("defgrid"):
                return route_keyword
        return None

    @property
    def scf_tol(self):
        """Extract SCF convergence tolerance from route keywords."""
        for route_input in self.route_inputs:
            if any(conv in route_input for conv in ORCA_SCF_CONVERGENCE):
                if route_input.endswith("scf"):
                    return route_input[:-3]
                return route_input
        return None

    @property
    def scf_algorithm(self):
        """Extract SCF algorithm specification from route keywords."""
        for route_input in self.route_inputs:
            if any(alg in route_input for alg in ORCA_ALL_SCF_ALGORITHMS):
                return route_input
        return None

    @property
    def jobtype(self):
        """
        Extract job type from route keywords.

        Returns:
            str: Job type (e.g., 'opt', 'freq') or 'sp' as default
        """
        for route_input in self.route_inputs:
            if route_input in ORCA_ALL_JOB_TYPES:
                return route_input
        return "sp"

    @property
    def freq(self):
        """Check if frequency calculation is requested."""
        for route_input in self.route_inputs:
            if "freq" in route_input:
                return route_input == "freq"
        return False

    @property
    def numfreq(self):
        """Check if numerical frequency calculation is requested."""
        for route_input in self.route_inputs:
            if "freq" in route_input:
                return route_input == "numfreq"
        return False
