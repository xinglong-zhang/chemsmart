import logging
import re

from chemsmart.io.orca import (
    ORCA_ALL_AB_INITIO,
    ORCA_ALL_AUXILIARY_BASIS_SETS,
    ORCA_ALL_BASIS_SETS_LOWER,
    ORCA_ALL_DISPERSION_CORRECTIONS,
    ORCA_ALL_EXTRAPOLATION_BASIS_SETS,
    ORCA_ALL_FUNCTIONALS,
    ORCA_ALL_JOB_TYPES,
    ORCA_ALL_QM2_BUILT_IN_METHODS,
    ORCA_ALL_SCF_ALGORITHMS,
    ORCA_ALL_SOLVENT_MODELS,
    ORCA_SCF_CONVERGENCE,
    is_orca_neb_joboption,
    normalize_orca_neb_joboption,
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
            route_string (str): ORCA route string
            to parse (e.g., '! B3LYP def2-TZVP')
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
            density_fitted_method = (
                route_keyword[3:] if route_keyword.startswith("ri-") else None
            )
            if route_keyword in {"rimp2", "mp2ri"}:
                return "mp2"
            if (
                route_keyword in ORCA_ALL_AB_INITIO
                or density_fitted_method in ORCA_ALL_AB_INITIO
                or route_keyword.startswith(("dlpno-cc", "ro-dlpno-cc"))
            ):
                return density_fitted_method or route_keyword
        return None

    @property
    def semiempirical(self):
        """Extract a built-in semiempirical or xTB method.

        The ORCA writer emits the project-facing GFN aliases as native
        ``XTB0``/``XTB1``/``XTB2`` keywords.  Parse both vocabularies back to
        the same canonical setting so a generated input retains its method
        meaning when ChemSmart independently reads it for preview validation.
        """

        from chemsmart.jobs.orca.settings import (
            ORCA_SEMIEMPIRICAL_ALIASES,
            _normalize_orca_semiempirical,
        )

        recognized = {
            *ORCA_SEMIEMPIRICAL_ALIASES,
            *(
                value.casefold()
                for value in ORCA_SEMIEMPIRICAL_ALIASES.values()
            ),
        }
        for route_keyword in self.route_keywords:
            if route_keyword.casefold() in recognized:
                return _normalize_orca_semiempirical(route_keyword)
        return None

    @property
    def method(self):
        """Extract the computational method from the simple-input route."""
        return self.functional or self.ab_initio or self.semiempirical

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
    def mdci_cutoff(self):
        """Return the project-facing DLPNO accuracy preset."""

        by_keyword = {
            "loosepno": "loose",
            "normalpno": "normal",
            "tightpno": "tight",
        }
        for route_keyword in self.route_keywords:
            value = by_keyword.get(route_keyword)
            if value is not None:
                return value
        return None

    @property
    def basis(self):
        """
        Extract basis set from route keywords.

        Returns:
            str: Basis set name or None if not found
        """
        for route_keyword in self.route_keywords:
            if route_keyword in ORCA_ALL_BASIS_SETS_LOWER:
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
    def aux_basis(self):
        """The declared-parameter spelling of :attr:`auxiliary_basis`.

        Project YAML and the capability domains call this ``aux_basis``; the
        round-trip contract looks a parameter up on the route by its declared
        name, and the historic property name silently failed that lookup.
        """

        return self.auxiliary_basis

    @property
    def solvent_model(self):
        """Implicit-solvation model named on the route, if any.

        ChemSmart writes every ORCA solvation request onto the route as
        ``MODEL`` or ``MODEL(solvent)`` -- CPCM(water), SMD(dmso), CPCMC,
        COSMORS. No parser read it back, so the preview validator compared a
        requested solvent model against nothing: the exact shape that once
        argued a session out of a correct ``ri_approximation: none``.
        """

        for route_keyword in self.route_keywords:
            head = route_keyword.split("(", 1)[0]
            if head in ORCA_ALL_SOLVENT_MODELS:
                return head
        return None

    @property
    def solvent_id(self):
        """Solvent named inside the solvation model's parentheses, if any.

        Read from the raw lower-cased route rather than token-wise, because
        23 ORCA solvent names contain spaces ("carbon disulfide") and the
        tokenizer splits them.
        """

        match = re.search(
            r"\b(?:cpcmc|cpcm|smd|cosmors)\(([^)]+)\)", self.route_string
        )
        if match is None:
            return None
        return match.group(1).strip()

    @property
    def auxiliary_basis_role(self):
        """Return the scientific fitting role of the auxiliary basis."""

        from chemsmart.jobs.orca.settings import _orca_auxiliary_basis_role

        return _orca_auxiliary_basis_role(self.auxiliary_basis)

    @property
    def correlation_auxiliary_basis(self):
        """Return an AuxC/AutoAux token, never a Coulomb-only fitting set."""

        if self.auxiliary_basis_role in {"correlation", "autoaux"}:
            return self.auxiliary_basis
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
            if route_keyword.startswith("defgrid") or (
                route_keyword.startswith("grid")
                and route_keyword[4:].isdigit()
            ):
                return route_keyword
        return None

    @property
    def ri_approximation(self):
        """Extract the resolution-of-the-identity choice from the route.

        The writer emits ``NoRI``, ``RI``, ``RIJCOSX`` or ``RIJK`` from the
        project's ``ri_approximation``.  Nothing read them back, so a preview
        validator that parses the generated input and compares it with the
        requested settings reported a critical finding for a correct project.

        Observed live: an n-butane session reproducing a protocol that calls
        for conventional four-index integrals wrote ``ri_approximation: none``,
        was told the generated input was invalid, and cleared the finding by
        deleting the key -- which is the one edit that lets ORCA's own default
        density fitting back in.  The harness pushed the model away from the
        scientifically correct input.
        """

        from chemsmart.jobs.orca.settings import ORCA_RI_KEYWORDS

        by_keyword = {
            keyword.lower(): name for name, keyword in ORCA_RI_KEYWORDS.items()
        }
        for route_keyword in self.route_keywords:
            match = by_keyword.get(route_keyword.lower())
            if match is not None:
                return match
        # ``RI-MP2`` owns correlation fitting and does not require the bare
        # RI-J keyword.  Preserve ChemSmart's typed ``ri`` setting when a
        # generated native input is read back during fake preview.
        if any(
            keyword in {"ri-mp2", "rimp2", "mp2ri"}
            for keyword in self.route_keywords
        ):
            return "ri"
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
        for route_input in self.route_keywords:
            if is_orca_neb_joboption(route_input):
                return "neb"
            if route_input in {"optts", "scants"}:
                return "ts"
            # IRC is a native ORCA simple-input job keyword, but it is not
            # listed in the small legacy ORCA_JOB_TYPES reference table.
            # Treating it as the default SP made a correctly materialized
            # ``! IRC`` input fail ChemSmart's own preview round trip.
            if route_input == "irc":
                return "irc"
            if route_input in ORCA_ALL_JOB_TYPES:
                return route_input
        return "sp"

    @property
    def neb_joboption(self):
        """Return the native NEB variant from the simple-input route."""

        for route_input in self.route_keywords:
            if is_orca_neb_joboption(route_input):
                return normalize_orca_neb_joboption(route_input)
        return None

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

    @property
    def vpt2(self):
        """Check if ORCA's anharmonic VPT2 analysis is requested."""

        return "vpt2" in self.route_keywords

    @property
    def qmmm_jobtype(self):
        for route_keyword in self.route_keywords:
            if "qm/" in route_keyword or "qmmm" in route_keyword:
                return route_keyword

    @property
    def qm_functional(self):
        return self.functional

    @property
    def qm_basis(self):
        return self.basis

    @property
    def qm2_method(self):
        # only available when QM2 methods are ORCA built-in methods
        for route_keyword in self.route_keywords:
            if "qm" in route_keyword:
                qmmm_methods = route_keyword.split("/")
                for method in qmmm_methods:
                    if method in ORCA_ALL_QM2_BUILT_IN_METHODS:
                        return method
