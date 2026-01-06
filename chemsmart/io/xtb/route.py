import logging

from chemsmart.io.xtb import XTB_ALL_METHODS

logger = logging.getLogger(__name__)


class XTBRoute:
    """
    Parser for XTB route (program call) specifications.

    This class parses XTB program call strings to extract job types and
    calculation parameters.

    Args:
        route_string (str): XTB program call string
    """

    def __init__(self, route_string):
        """
        Initialize XTB route parser.

        Args:
            route_string (str): XTB program call string to parse (e.g., "xtb coord --opt")
        """
        self.route_string = route_string.lower()
        self.route_inputs = self.route_string.split()

    @property
    def gfn_version(self):
        """
        Extract GFN version from the route.

        XTB GFN versions:
        - 'gfn0': GFN0-xTB
        - 'gfn1': GFN1-xTB
        - 'gfn2': GFN2-xTB
        - 'gfnff': GFN-FF

        Returns:
            str: GFN version identifier ('gfn0', 'gfn1', 'gfn2', 'gfnff') or None
        """
        for gfn_version in XTB_ALL_METHODS:
            if gfn_version in self.route_string:
                return gfn_version
        return None

    @property
    def job_type(self):
        """
        Extract the primary job type from the route.

        XTB job types:
        - 'opt': Geometry optimization (--opt, --ohess)
        - 'sp': Single point calculation (default)

        Returns:
            str: Job type identifier ('opt' or 'sp')
        """
        return self.get_job_type()

    def get_job_type(self):
        """
        Extract job type from route specification.

        Returns:
            str: Job type ('opt' or 'sp')
        """
        # Check for optimization keywords
        if "--opt" in self.route_string or "--ohess" in self.route_string:
            return "opt"
        return "sp"

    @property
    def freq(self):
        """
        Check if frequency calculation is requested.

        Returns:
            bool: True if frequency calculation is requested
        """
        return self.get_frequency()

    def get_frequency(self):
        """
        Check for frequency calculation in route.

        XTB frequency keywords: --hess, --ohess

        Returns:
            bool: True if 'hess' appears in route string
        """
        return "hess" in self.route_string
