import logging

from chemsmart.io.gaussian import (
    GAUSSIAN_AB_INITIO,
    GAUSSIAN_ADDITIONAL_OPT_OPTIONS,
    GAUSSIAN_ADDITIONAL_ROUTE_PARAMETERS,
    GAUSSIAN_BASES,
    GAUSSIAN_DIEZE_TAGS,
    GAUSSIAN_FUNCTIONALS,
    GAUSSIAN_SEMIEMPIRICAL,
)
from chemsmart.io.gaussian import (
    GAUSSIAN_SOLVATION_MODELS as gaussian_solvation_models,
)

logger = logging.getLogger(__name__)


class GaussianRoute:
    def __init__(self, route_string):
        if "\n" in route_string:
            route_string = route_string.replace("\n", " ")
        self.route_string = route_string.lower()
        self.route_inputs = self.route_string.split()

    @property
    def dieze_tag(self):
        return self.get_dieze_tag()

    @property
    def job_type(self):
        return self.get_job_type()

    @job_type.setter
    def job_type(self, value):
        self._job_type = value

    @property
    def freq(self):
        return self.get_freqeuncy()

    @property
    def numfreq(self):
        return self.get_numfreq()

    @property
    def force(self):
        return "force" in self.route_string

    @property
    def ab_initio(self):
        return self.get_ab_initio()

    @property
    def functional(self):
        functional, _ = self.get_functional_and_basis()
        return functional

    @property
    def basis(self):
        _, basis = self.get_functional_and_basis()
        return basis

    @property
    def semiempirical(self):
        for each_input in self.route_inputs:
            if any(
                semiemp.lower() in each_input
                for semiemp in GAUSSIAN_SEMIEMPIRICAL
            ):
                return each_input.upper()
        return None

    @property
    def solvent_model(self):
        return self.get_solvent_model()

    @property
    def solvent_id(self):
        return self.get_solvent_id()

    @property
    def additional_solvent_options(self):
        return self.get_additional_solvent_options()

    @property
    def solv(self):
        return self.solvent_model is not None and self.solvent_id is not None

    @property
    def additional_opt_options_in_route(self):
        return self.get_additional_opt_options()

    @property
    def additional_route_parameters(self):
        return self.get_additional_route_parameters()

    def get_dieze_tag(self):
        dieze_tag = None
        # dieze_tag '# ', '#N', '#P' '#T'
        if "#" in self.route_string and any(
            self.route_string.startswith(tag) for tag in GAUSSIAN_DIEZE_TAGS
        ):
            dieze_tag = self.route_string[0:2]
        return dieze_tag

    def get_job_type(self):
        # get job type: opt/ts/sp/ircf/ircr
        if "ts" in self.route_string:
            job_type = "ts"
        elif (
            "opt" in self.route_string
            and "ts" not in self.route_string
            and "modred" not in self.route_string
            and "stable=opt" not in self.route_string
        ):
            job_type = "opt"
        elif "opt=modred" in self.route_string:
            job_type = "modred"  # would include scan jobs too
        elif "irc" in self.route_string and "forward" in self.route_string:
            job_type = "ircf"
        elif "irc" in self.route_string and "reverse" in self.route_string:
            job_type = "ircr"
        elif "output=wfn" in self.route_string:
            job_type = "nci"
        elif (
            "Pop=MK IOp(6/33=2,6/41=10,6/42=17,6/50=1)".lower()
            in self.route_string
        ):
            job_type = "resp"
        elif "stable=opt" in self.route_string:
            job_type = (
                "link"  # so far only using stable=opt to determine link job
            )
        else:
            job_type = "sp"
        return job_type

    def get_freqeuncy(self):
        # get freq: T/F
        return "freq " in self.route_string

    def get_numfreq(self):
        return "freq=numer" in self.route_string

    def get_ab_initio(self):
        # get ab initio method by looking through the route string
        ab_initio = None
        for each_input in self.route_inputs:
            if any(ab in each_input for ab in GAUSSIAN_AB_INITIO):
                ab_initio = each_input
        return ab_initio

    def get_functional_and_basis(self):
        functional = None
        basis = None
        for each_input in self.route_inputs:
            # obtain functional and basis
            if "/" in each_input:
                func_basis = each_input.split(
                    "/"
                )  # TODO # not necessarily for non-standard route e.g.
                # pbepbe 6-31g(d,p)/auto force scrf=(dipole,solvent=water) pbc=gammaonly'
                if len(func_basis) == 2:
                    functional = func_basis[0]
                    basis = func_basis[1]
                elif len(func_basis) == 3:  # e.g., tpsstpss/def2tzvp/fit
                    functional = func_basis[0]
                    basis = f"{func_basis[1]}/{func_basis[2]}"  # note if the basis set for density fitting is written
                    # as 'def2tzvp fit', then the job fails to run
            else:  # '/' not in route
                if any(
                    functional in each_input
                    for functional in GAUSSIAN_FUNCTIONALS
                ):
                    functional = each_input
                if "empiricaldispersion" in each_input:
                    dispersion = each_input
                    functional_with_dispersion = functional + " " + dispersion
                    functional = functional_with_dispersion
                if (
                    any(basisset in each_input for basisset in GAUSSIAN_BASES)
                    and "generic" not in each_input
                ):
                    basis = each_input

        # non standard input by user e.g., `#wb897xd` without space
        if functional is not None and functional.startswith("#"):
            functional = functional[1:]

        return functional, basis

    def get_additional_route_parameters(self):
        additional_route = [
            each_input
            for each_input in self.route_inputs
            if any(
                route_parameter in each_input
                for route_parameter in GAUSSIAN_ADDITIONAL_ROUTE_PARAMETERS
            )
        ]

        return (
            " ".join(additional_route) if len(additional_route) != 0 else None
        )

    def get_additional_opt_options(self):
        additional_opt_options = []
        for each_input in self.route_inputs:
            if "opt" in each_input:
                opt_route = (
                    each_input.replace("opt=", "opt")
                    .split("opt")[-1]
                    .replace("(", "")
                    .replace(")", "")
                )
                if len(opt_route) != 0:
                    opt_options = opt_route.split(",")
                    for opt_option in opt_options:
                        # if 'eigentest' in opt_option and 'ts' in route_input:
                        #     additional_opt_options.append(opt_option)   # add `no/eigentest` only for ts jobs
                        # <-- `eigentest` already included in writing in GaussianSettings.write_gaussian_input()
                        if any(
                            option in opt_option
                            for option in GAUSSIAN_ADDITIONAL_OPT_OPTIONS
                        ):
                            additional_opt_options.append(
                                opt_option
                            )  # noqa: PERF401

        return (
            ",".join(additional_opt_options)
            if len(additional_opt_options) != 0
            else None
        )

    def get_solvent_model(self):
        if "scrf" in self.route_string:
            scrf_string = ""
            for each_input in self.route_inputs:
                if "scrf" in each_input:
                    scrf_string = each_input

            # get solvation model e.g., scrf=(cpcm,solvent=toluene)
            # if none of the models are present, set default model to pcm as in Gaussian
            if all(
                model not in scrf_string for model in gaussian_solvation_models
            ):
                solvent_model = "pcm"
            else:
                # some model is present, then parse route_input
                solvent_model = (
                    scrf_string.strip()
                    .split("(")[-1]
                    .strip()
                    .split(",")[0]
                    .strip()
                )
            return solvent_model
        return None

    def get_solvent_id(self):
        if "scrf" in self.route_string:
            scrf_string = ""
            for each_input in self.route_inputs:
                if "scrf" in each_input:
                    scrf_string = each_input

                    # get solvent identity
                    if "solvent" in scrf_string:
                        solvent_id = (
                            scrf_string.strip()
                            .split("solvent=")[-1]
                            .split(",")[0]
                            .split(")")[0]
                        )
                        if "read" in each_input:
                            # include read in solvent_id
                            solvent_id = f"{solvent_id},read"
                        return solvent_id
            return None
        return None

    def get_additional_solvent_options(self):
        if "scrf" in self.route_string:
            scrf_string = ""
            for each_input in self.route_inputs:
                if "scrf" in each_input:
                    scrf_string = each_input
                    scrf_line = scrf_string.split("(")[-1].split(")")[0]
                    scrf_line_elements = [
                        e.strip() for e in scrf_line.split(",")
                    ]
                    if len(scrf_line_elements) <= 2:
                        return None

                    # Identify elements to remove
                    elements_to_remove = set()

                    # Remove solvent model and solvent id
                    for element in scrf_line_elements:
                        if element in gaussian_solvation_models:
                            elements_to_remove.add(element)
                        if element.startswith(
                            "solvent="
                        ):  # Check if it's the solvent id
                            elements_to_remove.add(element)
                        if element == "read":
                            # included as part of solvent id
                            elements_to_remove.add(element)

                    # Remove identified elements
                    filtered_elements = [
                        e
                        for e in scrf_line_elements
                        if e not in elements_to_remove
                    ]
                    return (
                        ",".join(filtered_elements)
                        if filtered_elements
                        else None
                    )
        return None
