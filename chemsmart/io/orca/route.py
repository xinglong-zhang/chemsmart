import logging

from chemsmart.io.orca import (
    ORCA_ALL_AB_INITIO,
    ORCA_ALL_AUXILIARY_BASIS_SETS,
    ORCA_ALL_BASIS_SETS,
    ORCA_ALL_DISPERSION_CORRECTIONS,
    ORCA_ALL_EXTRAPOLATION_BASIS_SETS,
    ORCA_ALL_FUNCTIONALS,
)

logger = logging.getLogger(__name__)


class ORCARoute:
    def __init__(self, route_string):
        self.route_string = route_string.lower()
        self.route_inputs = self.route_string.split()

    @property
    def route_keywords(self):
        route_keywords = []
        for raw_route_input in self.route_inputs:
            if len(raw_route_input) != 0:
                route_input = raw_route_input.replace("!", "").strip()
                if len(route_input) != 0:
                    route_keywords.append(route_input)
        return route_keywords

    @property
    def functional(self):
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_FUNCTIONALS:
                return route_keyword
        return None

    @property
    def ab_initio(self):
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_AB_INITIO:
                return route_keyword
        return None

    @property
    def dispersion(self):
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_DISPERSION_CORRECTIONS:
                return route_keyword
        return None

    @property
    def basis(self):
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_BASIS_SETS:
                return route_keyword
        return None

    @property
    def auxiliary_basis(self):
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_AUXILIARY_BASIS_SETS:
                return route_keyword
        return None

    @property
    def extrapolation_basis(self):
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_EXTRAPOLATION_BASIS_SETS:
                return route_keyword
        return None

    @property
    def defgrid(self):
        for route_keyword in self.route_keywords:
            if route_keyword.startswith("defgrid"):
                return route_keyword
        return None

    @property
    def scf_tol(self):
        from chemsmart.io.orca import ORCA_SCF_CONVERGENCE as convs

        for route_input in self.route_inputs:
            if any(conv in route_input for conv in convs):
                if route_input.endswith("scf"):
                    return route_input[:-3]
                return route_input
        return None

    @property
    def scf_algorithm(self):
        from chemsmart.io.orca import ORCA_ALL_SCF_ALGORITHMS as algs

        for route_input in self.route_inputs:
            if any(alg in route_input for alg in algs):
                return route_input
        return None

    @property
    def job_type(self):
        from chemsmart.io.orca import ORCA_ALL_JOB_TYPES as ORCA_JOB_TYPES

        for route_input in self.route_inputs:
            if route_input in ORCA_JOB_TYPES:
                return route_input
        return "sp"

    @property
    def freq(self):
        for route_input in self.route_inputs:
            if "freq" in route_input:
                return route_input == "freq"
        return False

    @property
    def numfreq(self):
        for route_input in self.route_inputs:
            if "freq" in route_input:
                return route_input == "numfreq"
        return False
