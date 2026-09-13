"""Guides: the leaves of the stem-and-leaf surface.

The stem is what every session reads: the core tools, the core operation
vocabulary, and the universal rules. A guide is a family unit -- extra
tools and operations that only that family needs, a few hundred words of
guidance, and the rules placed on it -- opened by the host from four
signals (the task text, the workspace, the planned DAG, the previous
run's endings) or pulled by the model with ``open_guide``. Opening a
guide changes what the model can express and how much it reads; it never
changes what the host approves or verifies.

Every activation is an event carrying the signal and the new tool-schema
digest, so a reading can count them.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, Mapping

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.terminal_states import REPAIRABLE_NODE_STATES

#: The typed terminal states that open the recovery guide.
_RECOVERY_STATES = REPAIRABLE_NODE_STATES


@dataclass(frozen=True)
class GuideV1:
    guide_id: str
    title: str
    tier: str
    body: str
    tools: tuple[str, ...] = ()
    operations: tuple[str, ...] = ()
    activation_terms: tuple[str, ...] = ()
    jobtypes: tuple[str, ...] = ()
    analysis_kinds: tuple[str, ...] = ()
    terminal_states: tuple[str, ...] = ()
    workspace_kinds: tuple[str, ...] = ()
    #: Programs whose presence in the planned DAG opens this guide. A
    #: program leaf cannot key on a jobtype (``hess`` is xTB's too) and
    #: the crossprogram guide keys on the *count* of programs, so the
    #: fourth signal a program family needs is the program itself.
    programs: tuple[str, ...] = ()
    rule_placement: str = field(default="")

    def __post_init__(self) -> None:
        if not self.guide_id or " " in self.guide_id:
            raise ContractError(
                f"guide id must be one token: {self.guide_id!r}"
            )
        object.__setattr__(
            self,
            "rule_placement",
            self.rule_placement or f"leaf:{self.guide_id}",
        )

    def index_entry(self) -> str:
        return f"{self.guide_id}: {self.title}"


GUIDES: tuple[GuideV1, ...] = (
    GuideV1(
        guide_id="structure",
        title="building a starting structure: compose, derive, edit, append, displace",
        tier="T3",
        tools=(
            "compose_molecular_arrangement",
            "derive_molecular_species",
            "edit_molecular_geometry",
            "append_molecular_atom",
            "displace_along_vibrational_mode",
            "break_symmetry",
        ),
        activation_terms=(
            "complex",
            "fragment",
            "protonat",
            "deprotonat",
            "homoly",
            "radical",
            "dihedral",
            "torsion",
            "rotamer",
            "conformer",
            "transition state",
            "saddle",
            "build",
            "arrange",
            "dimer",
            "adduct",
        ),
        jobtypes=("ts",),
        body=(
            "Every host-built geometry is a starting structure, never a "
            "relaxed one: compose places two identity-bound fragments at one "
            "or two explicit atomic contacts; derive keeps an ordered subset "
            "of one parent's atoms (homolysis, deprotonation, fragment "
            "extraction are that one operation); edit sets one internal "
            "coordinate as a rigid motion of a side you name; append places "
            "one atom by three internal coordinates against three anchors; "
            "displace steps a frequency-bearing result along one of its own "
            "printed modes. None of them infers an electronic state -- "
            "removing a hydrogen gives a radical or an anion depending on "
            "where its electron went -- so bind charge and multiplicity "
            "explicitly afterwards, and the consuming stage is a new "
            "workflow for review. A requested value is never refused on "
            "scientific merit: the optimisation that consumes the structure "
            "grades it, and requested-versus-relaxed is the delivered "
            "observable. For a transition-state guess, place the forming "
            "bonds with compose's second contact rather than editing a "
            "distance between separate fragments, which is not an edit at "
            "all. Read which atoms move in a mode before you name it."
        ),
    ),
    GuideV1(
        guide_id="scan",
        title="relaxed coordinate scans and what to do with the surface",
        tier="T2",
        tools=("bind_scan_point_geometry",),
        # coordinate_at_minimum and coordinate_at_maximum used to be this
        # guide's operations. "Which of these is lowest" is every
        # comparison question, not a scan question: a spin-state task
        # declared the ground state's multiplicity as a required output,
        # could not open this guide by any of its terms, was refused for
        # want of a producer, and built the label itself (NOVEL-2 ino1,
        # 2026-09-04). They are stem operations now.
        operations=(),
        activation_terms=(
            "scan",
            "profile",
            "surface",
            "rotational barrier",
            "torsional",
        ),
        jobtypes=("scan",),
        body=(
            "A relaxed scan's driven coordinate is a fact about this molecule "
            "in this calculation: it lives on the workflow node "
            "(internal_coordinates on compile_command), not in project YAML. "
            "A scan ends at a surface, and which point travels is a "
            "scientific judgement. The validated minimum-energy sampled "
            "point may feed a downstream optimisation inside one approval "
            "through the declared producer edge; any other point is an "
            "explicit scan-point binding whose consuming stage is a new "
            "workflow. coordinate_at_minimum and coordinate_at_maximum read "
            "the surface's own ordered vectors -- a barrier position cannot "
            "come from max alone. A step that failed to converge leaves the "
            "surface so far readable as it stands; say how far it reached. A "
            "scan grid does not locate a stationary point: optimise the "
            "minima you find and characterise them by frequencies before "
            "calling any of them a minimum, and count the distinct "
            "stationary points with their orders. scan_coordinate_values is "
            "dimensionless -- positional numbers in the scan's own unit -- "
            "so declare it with unit '1'; a physical distance or angle is "
            "measured from a delivered geometry with the "
            "distance/angle/dihedral operations."
        ),
    ),
    GuideV1(
        guide_id="constants",
        title="literature constants, conventions, pKa and redox potentials",
        tier="T3",
        operations=("constant", "gibbs_to_pka", "gibbs_to_redox_potential"),
        activation_terms=(
            "pka",
            "acid",
            "base",
            "redox",
            "reduction potential",
            "oxidation",
            "electrode",
            "standard state",
            "proton",
            "pcet",
            "electrochem",
        ),
        body=(
            "A value the record supplies rather than the calculation -- the "
            "aqueous proton free energy, a standard-state correction, a "
            "reference acid's measured pKa, an electrode's absolute "
            "potential -- is selected by registered name through the "
            "constant operation; the host owns the value, unit and "
            "standard-state convention, and a literal is recorded as "
            "model-authored. Constants that look independent are often "
            "matched pairs: read the convention family and the purpose "
            "phrase before combining two, and prefer a registered composed "
            "value where one exists. A family says nothing about standard "
            "state, so two entries on one scale can still need the term "
            "that bridges them. gibbs_to_pka owns pKa = dG/(RT ln 10); "
            "gibbs_to_redox_potential owns E = -dG/(nF) with the IUPAC sign, "
            "so a favourable reduction has a negative free energy and a "
            "positive potential, and referencing an electrode stays ordinary "
            "subtraction so the electrode you chose stays visible. Its n is "
            "the electron count, and you can derive it instead of typing "
            "it: subtract the two states' own charge selectors. A typed n "
            "is a number of yours, so nothing downstream of it can serve as "
            "measured evidence for an uncertainty; a derived one keeps the "
            "whole chain the host's. Continuum "
            "solvation of a small localised anion carries a documented "
            "systematic of roughly ten kcal/mol; state it beside the number "
            "and license no accuracy claim."
        ),
    ),
    GuideV1(
        guide_id="cbs",
        title="complete-basis-set extrapolation",
        tier="T3",
        operations=(
            "exponential_cbs_limit",
            "scf_exponential_cbs_limit",
            "scf_inverse_power_cbs_limit",
            "correlation_inverse_power_cbs_limit",
        ),
        activation_terms=(
            "basis-set limit",
            "basis set limit",
            "cbs",
            "extrapolat",
            "complete basis",
        ),
        body=(
            "The basis-set limit is one named operation, not fifteen "
            "arithmetic nodes: a session that rebuilt the three-point "
            "exponential form from multiply, subtract and divide nodes was "
            "the reason these operations exist. SCF and correlation energies "
            "converge by different laws -- exponential and inverse-power "
            "respectively -- so extrapolate them separately and add, never "
            "the total energy under one law. The cardinal numbers must be "
            "consecutive and the exponent, where one is required, comes from "
            "the method's own protocol and is recorded as such -- supply "
            "extrapolation_exponent only when the protocol you are "
            "reproducing states it. When the protocol just says the "
            "energy was extrapolated exponentially and you have three "
            "successive cardinal numbers, prefer exponential_cbs_limit: "
            "it fits the decay from the data and introduces no constant "
            "of your own."
        ),
    ),
    GuideV1(
        guide_id="ensemble",
        title="conformer ensembles and Boltzmann averaging",
        tier="T2",
        operations=("boltzmann_populations", "boltzmann_average"),
        activation_terms=(
            "boltzmann",
            "ensemble",
            "population",
            "conformer",
            "averag",
        ),
        body=(
            "A conformer set is a sample, not the ensemble; say what was "
            "sampled. Populations come from free energies at the stated "
            "temperature with the degeneracy of each multiply-realisable "
            "state (an enantiomeric pair counts twice). A Boltzmann average "
            "of a vector magnitude is linear in the property unless the "
            "observable is a mean square; say which you took. A 0.0000 "
            "energy tie between mirror-image minima is correct physics, not "
            "a defect."
        ),
    ),
    GuideV1(
        guide_id="spectroscopy",
        title="rotational constants, moments of inertia, excitations",
        tier="T1",
        operations=(
            "principal_moments_of_inertia",
            "rigid_rotor_constants",
            "linear_rotor_constant",
            "center_of_mass",
            "photon_wavelength",
        ),
        activation_terms=(
            "rotational constant",
            "moment of inertia",
            "microwave",
            "excitation",
            "absorption",
            "wavelength",
            "oscillator strength",
            "td-dft",
            "tddft",
            "uv",
        ),
        jobtypes=("td",),
        body=(
            "Rotational constants follow from the principal moments of the "
            "optimised geometry; a linear molecule has one constant and its "
            "own operation. Excited-state selectors answer per manifold "
            "root, singlet and triplet apart, with oscillator strengths "
            "beside energies; a wavelength is the photon operation on an "
            "excitation energy, never a hand conversion. PySCF stores "
            "excitation energies in hartree where the log-parsing programs "
            "print electronvolts; the reader states its unit and the "
            "arithmetic is canonical."
        ),
    ),
    GuideV1(
        guide_id="database",
        title="workspace databases and batches",
        tier="T3",
        tools=("inspect_database_records", "extract_database_record_geometry"),
        activation_terms=(
            "database",
            ".db",
            "records",
            "batch",
            "each record",
            "every entry",
        ),
        workspace_kinds=("chemsmart_db",),
        body=(
            "A workspace database is an inspectable artifact whose stored "
            "fields -- charge, multiplicity, energy, optimised flags -- are "
            "observations from the records' own provenance, never bindings. "
            "Enumerate the records, extract one record's exact coordinates "
            "into a lineage-carrying geometry artifact, and bind identity "
            "and electronic state explicitly per record; a stored state that "
            "contradicts the electron count is flagged loudly in the review, "
            "not silently corrected and not silently copied. N records plan "
            "as N disconnected sub-DAGs in one workflow under one decision; "
            "execution is record-major, one record's failure settles that "
            "record while the others deliver, and there is deliberately no "
            "aggregate quantity: a batch of N is N observations."
        ),
    ),
    GuideV1(
        guide_id="crossprogram",
        title="geometries that cross programs, numbers that must not",
        tier="T3",
        activation_terms=(
            "xtb geometry",
            "cheaper",
            "semi-empirical",
            "single point on",
            "single points on",
            "on those geometries",
            "geometries from",
            "mixed level",
            "composite",
            "gfn2 then",
            "two programs",
        ),
        body=(
            "The optimised-geometry handoff is keyed on the producing "
            "program and refuses any change of atom identity or order, so an "
            "xTB optimisation may feed an ORCA or PySCF single point with "
            "parent atom i as child atom i. A typed value carries its unit "
            "and dimension, not the method that produced it, so a "
            "tight-binding energy and a hybrid-DFT energy subtract without "
            "complaint: mixing levels is a method when it is deliberate and "
            "a mistake when it is not, and the displayed chain names the "
            "level behind every input so the reviewer can tell. Naming the "
            "level is necessary and not sufficient: ORCA's B3LYP and PySCF's "
            "b3lyp differ in their local correlation (VWN5 versus VWN3) and "
            "gave total energies 0.24 hartree apart under identical strings; "
            "compare differences across programs, never totals, and say "
            "which variant each program means."
        ),
    ),
    GuideV1(
        guide_id="pyscf",
        title="PySCF: a library backend with one structure per result",
        tier="T1",
        activation_terms=("pyscf", "gpu4pyscf", "libxc"),
        programs=("pyscf",),
        workspace_kinds=("pyscf_hdf5",),
        body=(
            "PySCF is a library, not a binary: the host writes the driver, "
            "runs it in the registered PySCF interpreter, and the HDF5 "
            "result (pyscf_hdf5) is the program's typed account; its log is "
            "never read. Four executable stages, one node each: sp, opt, "
            "hess, td. Every quantity belongs to one structure, the final "
            "one: the SCF is re-converged there before anything is read. "
            "supplied_positions is what the run was handed; "
            "reached_positions (opt only) is the last geometry evaluated, "
            "converged or not; for sp, hess and td supplied and final "
            "coincide and the host checks it. A hess is judged by the 20 "
            "cm-1 rule like every program (one imaginary mode types the node "
            "failed_wrong_stationary_point); its D3/D4 and SMD cavity "
            "blocks are finite differences. No imaginary mode means none was "
            "found, not stationarity -- the outcome reports the gradient at "
            "the Hessian geometry. functional is the name the project asked "
            "for; b3lyp and b3lypg are one libxc functional (VWN3) and "
            "b3lyp5 the VWN5 form, so compare differences across programs, "
            "never totals. td (TDA or TDDFT on a Kohn-Sham reference) gives "
            "roots ascending within one manifold at this geometry: singlet "
            "or triplet on a closed shell, the one unrestricted manifold on "
            "an open shell, which prints no per-root <S^2>. A root is an "
            "index, never a state identity; select it with ref and indices. "
            "An opt carrying excited_state_root follows root k by index, "
            "re-evaluates the spectrum at the reached geometry and reports "
            "the gap to the ground state and to the neighbouring root -- a "
            "small gap is a sensor fact, not a verdict -- and no Hessian "
            "exists for that minimum here; gas phase only. PCM/SMD give td "
            "energies at a fixed geometry under PySCF's non-equilibrium "
            "response, whose dielectric (1.78 for every solvent) the "
            "artifact records. An unconverged root or CC amplitude set ends "
            "the node typed; td_max_cycle and cc_max_cycle are the repair "
            "controls. ab_initio mp2, ccsd and ccsd(t) run on an HF "
            "reference: sp for all, opt for mp2 and ccsd. reference_energy "
            "and correlation_energy are the program's own components and "
            "correlation_energy includes the triples; PySCF correlates "
            "every electron unless frozen_core (a count or auto) says "
            "otherwise, where ORCA and Gaussian freeze core by default. On "
            "every result the dipole, populations, orbital energies and "
            "<S^2> belong to the SCF reference and energy to the surface the "
            "job computed on; inspect_run says which beside each selector. "
            "Not available here: transition states, IRC, scans, a ccsd(t) "
            "optimisation, EOM or CASSCF, GPU execution, double hybrids, "
            "mixed basis/ECP, a Hessian on any excited or correlated "
            "surface or for ROHF/open-shell NLC; only the geometric "
            "optimiser is installed. A non-converged optimisation still "
            "writes its last structure. Free energies come from the host's "
            "RRHO engine (derive_thermochemistry; state T and p), never "
            "PySCF's thermo; the receipt states the isotope-averaged masses "
            "and the symmetry number behind them."
        ),
    ),
    GuideV1(
        guide_id="recovery",
        title="answering a run that failed or landed on the wrong stationary point",
        tier="T4",
        terminal_states=tuple(sorted(_RECOVERY_STATES)),
        tools=("bind_reached_geometry",),
        body=(
            "A failed run is evidence, and the wake context's repair_menu "
            "names, for each way a node ended, the ordinary route that "
            "answers it; the host names the route and the next run's physics "
            "grades it. Read the run's typed outcome (inspect_run) and the "
            "native findings before choosing. A wrong stationary point calls "
            "for reading which atoms carry the offending mode "
            "(vibrational_mode_atom_participation, checking "
            "vibrational_mode_degeneracy_group first), then displacing along "
            "it or editing the coordinate it moves; an SCF failure is a "
            "state question before it is a solver question; a convergence "
            "failure or timeout restarts from the reached geometry inside "
            "the remaining budget. A revision may change the structure or a "
            "setting the project exposes; it may not change identity, "
            "electronic state, or conditions -- those return to the human. "
            "A re-run of a failed node takes a fresh node id: its earlier "
            "directory is evidence and the plan refuses an id that already "
            "holds outputs. "
            "Recovering the structure does not recover numbers computed "
            "from the rejected one: re-derive and re-claim them. Standing by "
            "a result with a cited validation receipt is also an answer; "
            "leaving the failure unanswered is the one thing that is not. "
            "Beside every repair route stands a disposition: the ending may "
            "itself be the finding. A saddle where a minimum was promised is "
            "a stationary point of that surface with an energy -- an "
            "inversion, a symmetry breaking, a hidden reaction coordinate -- "
            "and the run outcome's anomalies name its mode and which heavy "
            "atoms carry it; an SCF that will not settle may be an "
            "instability; a geometry that walked away may have found another "
            "basin. Read the anomaly, say what the structure is in the "
            "decision citing its receipt, and then repair, stand by, or "
            "both. What the host records here is an observation, never a "
            "verdict; naming it is yours."
        ),
    ),
    GuideV1(
        guide_id="saddle",
        title="transition states, imaginary modes, and intrinsic reaction coordinates",
        tier="T5",
        tools=("characterise_stationary_point",),
        operations=("transition_state_crossover_temperature",),
        activation_terms=(
            "transition state",
            "saddle",
            "barrier",
            "irc",
            "reaction coordinate",
            "activation",
            "mechanism",
            "competing",
            "selectivity",
            "endo",
            "exo",
        ),
        jobtypes=("ts", "irc"),
        body=(
            "A transition-state search promises exactly one imaginary mode "
            "under the 20 cm-1 convention; the host judges every executed "
            "result on that promise and a mismatch is a typed failure, not a "
            "result to report. Seed a search from a validated "
            "frequency-bearing producer's Hessian where one exists; a "
            "hand-built guess is a starting structure. Which channel a "
            "saddle belongs to is decided from the atoms that carry its "
            "imaginary mode, never from how the guess was built or named. "
            "An IRC consumes the transition state's own geometry and "
            "analytic Hessian as role-distinct producer edges, one per "
            "direction; the path goes to an XYZ sidecar, so whether the "
            "saddle connects two particular minima is an observation you "
            "make from that trajectory, not a host-rendered claim. A "
            "barrier is stated relative to a reference you name and defend; "
            "a difference between two barriers at a small basis without "
            "dispersion licenses a direction, rarely a magnitude, and the "
            "delivery says so."
        ),
    ),
)

GUIDES_BY_ID: Mapping[str, GuideV1] = {
    guide.guide_id: guide for guide in GUIDES
}

#: Tools and operations that belong to some leaf; everything else is stem.
LEAF_TOOLS: Mapping[str, str] = {
    tool: guide.guide_id for guide in GUIDES for tool in guide.tools
}
LEAF_OPERATIONS: Mapping[str, str] = {
    operation: guide.guide_id
    for guide in GUIDES
    for operation in guide.operations
}


def guide_for_tool(tool_name: str) -> str | None:
    return LEAF_TOOLS.get(tool_name)


def guides_from_text(text: str) -> tuple[str, ...]:
    """Guides the task text asks for, by substring in the pack convention."""

    lowered = " ".join(str(text or "").lower().split())
    # A term matches as whole words: "base" inside "database" opened the
    # constants guide on a task with no constant in it (live, 2026-09-02).
    found = [
        guide.guide_id
        for guide in GUIDES
        if any(_term_in_text(term, lowered) for term in guide.activation_terms)
    ]
    return tuple(sorted(set(found)))


def _term_in_text(term: str, lowered: str) -> bool:
    pattern = r"(?<![a-z0-9])" + re.escape(term.lower()) + r"(?![a-z0-9])"
    return re.search(pattern, lowered) is not None


def guides_from_workspace(kinds: Iterable[str]) -> tuple[str, ...]:
    present = set(kinds)
    return tuple(
        sorted(
            guide.guide_id
            for guide in GUIDES
            if present.intersection(guide.workspace_kinds)
        )
    )


def guides_from_plan(
    *,
    jobtypes: Iterable[str] = (),
    operations: Iterable[str] = (),
    constants: Iterable[str] = (),
    tools: Iterable[str] = (),
    programs: Iterable[str] = (),
) -> tuple[str, ...]:
    """Guides a planned DAG needs: by jobtype, operation, tool, program."""

    found: set[str] = set()
    # Two programs in one DAG is the crossprogram guide's own signal; the
    # task text that named two programs never matched its phrase-shaped
    # terms (audit, 2026-09-03).
    if len({str(item).lower() for item in programs if str(item)}) >= 2:
        found.add("crossprogram")
    jobs = {str(item).lower() for item in jobtypes}
    ops = {str(item) for item in operations}
    named_programs = {str(item).lower() for item in programs if str(item)}
    for guide in GUIDES:
        if jobs.intersection(guide.jobtypes):
            found.add(guide.guide_id)
        if ops.intersection(guide.operations):
            found.add(guide.guide_id)
        if named_programs.intersection(guide.programs):
            found.add(guide.guide_id)
    if any(constants):
        found.add("constants")
    for tool in tools:
        owner = LEAF_TOOLS.get(str(tool))
        if owner:
            found.add(owner)
    return tuple(sorted(found))


def guides_from_states(states: Iterable[str]) -> tuple[str, ...]:
    present = {str(item) for item in states}
    found = {
        guide.guide_id
        for guide in GUIDES
        if present.intersection(guide.terminal_states)
    }
    if "failed_wrong_stationary_point" in present:
        found.update({"structure", "saddle"})
    return tuple(sorted(found))


def guide_index_sentence(active: Iterable[str] = ()) -> str:
    """The stem's one sentence about guides, naming what can be opened."""

    opened = set(active)
    listing = "; ".join(
        f"{guide.guide_id} ({guide.title})"
        + (" -- open" if guide.guide_id in opened else "")
        for guide in GUIDES
    )
    return (
        " Guides are family units of tools, operations and guidance the "
        "host opens when the task, the plan, or a previous run calls for "
        "them, and that you may open yourself with open_guide(guide_id); an "
        "opened guide's tools and operations join the surface on the next "
        "turn. Available: " + listing + "."
    )


__all__ = [
    "GUIDES",
    "GUIDES_BY_ID",
    "LEAF_OPERATIONS",
    "LEAF_TOOLS",
    "GuideV1",
    "guide_for_tool",
    "guide_index_sentence",
    "guides_from_plan",
    "guides_from_states",
    "guides_from_text",
    "guides_from_workspace",
]
