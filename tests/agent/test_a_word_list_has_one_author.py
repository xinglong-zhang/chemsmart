"""A list of jobtype or program words has one author.

PySCF ``td`` was registered in the reader's declaration and refused by a
workspace scan that hand-listed ``{"sp", "opt", "hess"}`` beside it; the
witness bank went red, the suite never did.  The general invariant is
that a set of jobtype or program words is the declaration everything
else derives from, a program's own native spelling, or an import of one
-- never a second copy that drifts.

This walks the syntax trees of the host layers and finds every literal
set, tuple or list whose members are all jobtype or program words.  Each
one must be recorded below with the reason it may stand: it is the
declaration, it is a program's own spelling, or it is a copy with the
commit that derives it named.  A literal nobody recorded fails, and a
record entry whose site is gone fails too, so the record shrinks as
copies are derived and can never quietly grow.
"""

from __future__ import annotations

import ast
import pathlib

import pytest

pytestmark = pytest.mark.capability("gate:resolver.one_answer_per_question")

REPO = pathlib.Path(__file__).resolve().parents[2]
ROOTS = (
    "chemsmart/agent",
    "chemsmart/analysis",
    "chemsmart/jobs/pyscf",
    "chemsmart/jobs/xtb",
    "chemsmart/settings",
    "chemsmart/io/pyscf",
)
JOBTYPE_WORDS = frozenset(
    {
        "freq",
        "hess",
        "irc",
        "ircf",
        "ircr",
        "link",
        "modred",
        "neb",
        "opt",
        "opt_freq",
        "qmmm",
        "scan",
        "sp",
        "td",
        "ts",
    }
)
PROGRAM_WORDS = frozenset({"gaussian", "orca", "pyscf", "xtb"})
WORDS = JOBTYPE_WORDS | PROGRAM_WORDS

_DECLARATION = "declaration: this is the authority others derive from"
_GAUSSIAN_NATIVE = "Gaussian's own ircf/ircr spellings inside one reader"
_GAS_SOLV = "the gas/solv project shape of the route-building programs"
_READER_PROGRAMS = "to derive from the reader registry (deferred, recorded)"
_SURFACING_SIGNAL = "declaration: which reference a planned job type surfaces"
_A4 = "derives in A.4: keyed on the Hessian stage the artifact records"

#: (path, enclosing symbol, assignment target, words) -> why it may stand.
RECORDED: dict[tuple[str, str, str | None, tuple[str, ...]], str] = {
    (
        "chemsmart/agent/capabilities.py",
        "<module>",
        "_VALIDATED_PROGRAMS",
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _READER_PROGRAMS,
    (
        "chemsmart/agent/knowledge.py",
        "<module>",
        "_PYSCF_SUBSTITUTION_JOB_TYPES",
        (
            "freq",
            "hess",
            "irc",
            "link",
            "modred",
            "neb",
            "opt",
            "opt_freq",
            "qmmm",
            "scan",
            "sp",
            "td",
            "ts",
        ),
    ): _DECLARATION
    + ": which PySCF job types carry which Gaussian job "
    "family, a judgement about chemistry; whether they execute is read "
    "from the capability registry",
    (
        "chemsmart/agent/terminal_states.py",
        "<module>",
        "FIXED_GEOMETRY_CURVATURE_JOBTYPES",
        ("freq", "hess"),
    ): _DECLARATION
    + ": the job types that measure a structure's curvature without "
    "moving it, and so inherit the promise of whatever produced that "
    "structure rather than making one of their own",
    (
        "chemsmart/agent/terminal_states.py",
        "<module>",
        "GEOMETRY_SEARCH_JOBTYPES",
        ("opt", "ts"),
    ): _DECLARATION
    + ": the job types that search for the structure",
    (
        "chemsmart/agent/terminal_states.py",
        "<module>",
        "STATIONARY_POINT_PROMISES",
        ("freq", "hess", "opt", "ts"),
    ): _DECLARATION
    + ": what each job type promises about imaginary "
    "modes; the coverage test below holds the searching job types "
    "inside it",
    (
        "chemsmart/agent/terminal_states.py",
        "<module>",
        "PATH_ACCOUNT_OBSERVATIONS",
        ("pyscf",),
    ): _DECLARATION
    + ": where each program's result validation files its own "
    "account of a path it walked. The key beside the program is that "
    "validator's field name, which only that validator can declare, so "
    "a second walking program is one line here rather than a second "
    "reader of the same shape",
    (
        "chemsmart/agent/terminal_states.py",
        "<module>",
        "SURFACE_SAMPLING_JOBTYPES",
        ("scan",),
    ): _DECLARATION
    + ": job types that sample a surface, whose points "
    "are not claimed stationary",
    (
        "chemsmart/agent/commands.py",
        "<module>",
        "_COORDINATE_DRIVEN_JOBTYPES",
        ("modred", "scan"),
    ): _DECLARATION
    + ": the jobtypes carrying a driven or held coordinate",
    (
        "chemsmart/agent/execution.py",
        "handoff_optimized_native_geometry",
        None,
        ("gaussian", "orca"),
    ): "the log-parsing programs; "
    + _READER_PROGRAMS,
    (
        "chemsmart/agent/execution.py",
        "<module>",
        "HESSIAN_CONSUMER_ROLES",
        ("irc",),
    ): _DECLARATION
    + ": which stage consumes the final Hessian of a "
    "converged transition state",
    (
        "chemsmart/agent/execution.py",
        "<module>",
        "HESSIAN_CONSUMER_ROLES",
        ("ts",),
    ): _DECLARATION
    + ": which stage produces that Hessian, and which "
    "stage a starting Hessian is fed to",
    (
        "chemsmart/agent/exposure.py",
        "<module>",
        "JOBTYPE_REFERENCES",
        ("irc", "scan", "ts"),
    ): _SURFACING_SIGNAL,
    (
        "chemsmart/agent/exposure.py",
        "<module>",
        "PROGRAM_REFERENCES",
        ("pyscf",),
    ): _SURFACING_SIGNAL,
    (
        "chemsmart/agent/knowledge.py",
        "assess_typed_program_substitution",
        None,
        ("gaussian", "pyscf"),
    ): _DECLARATION
    + ": the one admitted substitution pair",
    (
        "chemsmart/agent/live_session.py",
        "_scan_orca_result_artifacts",
        None,
        ("freq", "irc", "modred", "opt", "scan", "sp", "td", "ts"),
    ): "to derive from the ORCA reader's declared jobtypes (deferred)",
    (
        "chemsmart/agent/live_session.py",
        "_scan_gaussian_result_artifacts",
        None,
        ("freq", "irc", "link", "opt", "sp", "td", "ts"),
    ): "to derive from the Gaussian reader's declared jobtypes (deferred)",
    (
        "chemsmart/agent/live_session.py",
        "<module>",
        "_ROUTE_PROGRAM_STAGE_SECTIONS",
        ("link", "neb", "td"),
    ): _DECLARATION
    + ": route stages with their own project section",
    (
        "chemsmart/agent/live_session.py",
        "_validate_bounded_envelope_against_registry",
        "result_validators",
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _READER_PROGRAMS,
    (
        "chemsmart/agent/program_verifiers.py",
        "validate_preview_workspace",
        None,
        ("gaussian", "orca"),
    ): _GAS_SOLV,
    (
        "chemsmart/agent/projects.py",
        "project_section_application_observation",
        None,
        ("gaussian", "orca"),
    ): _GAS_SOLV,
    (
        "chemsmart/agent/projects.py",
        "_loader_project_section_sources",
        None,
        ("gaussian", "orca"),
    ): _GAS_SOLV,
    (
        "chemsmart/agent/terminal_states.py",
        "_native_failure",
        None,
        ("gaussian", "orca", "pyscf", "xtb"),
    ): "to derive from the native-failure tables (deferred)",
    (
        "chemsmart/agent/tool_runtime.py",
        "_neutral_sensor_facts",
        None,
        ("freq", "hess"),
    ): _A4,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/_evaluate_execution_outputs",
        None,
        ("gaussian", "orca", "xtb"),
    ): "the log-parsing readers; "
    + _READER_PROGRAMS,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/_evaluate_execution_outputs",
        None,
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _READER_PROGRAMS,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/_evaluate_execution_outputs",
        "expected_directions",
        ("ircf",),
    ): _GAUSSIAN_NATIVE,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/_evaluate_execution_outputs",
        "expected_directions",
        ("ircr",),
    ): _GAUSSIAN_NATIVE,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/_evaluate_execution_outputs",
        None,
        ("ircf", "ircr"),
    ): _GAUSSIAN_NATIVE,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/_evaluate_execution_outputs",
        "expected_directions",
        ("ircf", "ircr"),
    ): _GAUSSIAN_NATIVE,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/_evaluate_execution_outputs",
        None,
        ("irc", "td"),
    ): "Gaussian prints a jobtype word its request never used",
    (
        "chemsmart/agent/workflows.py",
        "CommandNodeIntentV1/_validate_internal_coordinates",
        None,
        ("scan", "ts"),
    ): _DECLARATION
    + ": the jobtypes a driven coordinate may ride",
    (
        "chemsmart/analysis/result_readers.py",
        "_irc_structures",
        None,
        ("irc", "ircf", "ircr"),
    ): _GAUSSIAN_NATIVE,
    (
        "chemsmart/jobs/pyscf/settings.py",
        "<module>",
        "PYSCF_JOBTYPES",
        ("hess", "irc", "opt", "sp", "td", "ts"),
    ): _DECLARATION,
    (
        "chemsmart/agent/execution.py",
        "<module>",
        "PATH_ENDPOINT_PRODUCER_STAGES",
        ("irc",),
    ): _DECLARATION
    + ": the path-walking stages whose endpoint may hand on, admitted per "
    "program by the reader's own declaration of a reached structure",
    (
        "chemsmart/agent/execution.py",
        "<module>",
        "CONSTRAINED_GEOMETRY_PRODUCER_STAGES",
        ("modred",),
    ): _DECLARATION
    + ": the stages that relax around a held coordinate and end on the "
    "structure at it, admitted per program by the reader's own "
    "declaration of a reached structure",
    (
        "chemsmart/agent/terminal_states.py",
        "<module>",
        "START_POINT_PROMISES",
        ("irc",),
    ): _DECLARATION
    + ": what a job type promises about the geometry it was handed",
    (
        "chemsmart/jobs/pyscf/settings.py",
        "<module>",
        "PYSCF_MOVING_STAGES",
        ("irc", "opt", "ts"),
    ): _DECLARATION
    + ": the stages that move the geometry they were handed",
    (
        "chemsmart/jobs/pyscf/settings.py",
        "<module>",
        "PYSCF_ANALYTIC_HESSIAN_STAGES",
        ("hess", "irc", "ts"),
    ): _DECLARATION
    + ": the stages whose first act can be PySCF's analytic Hessian",
    (
        "chemsmart/jobs/pyscf/validation.py",
        "<module>",
        "FIXED_GEOMETRY_JOBTYPES",
        ("hess", "sp", "td"),
    ): _DECLARATION
    + ": the jobtypes held to the geometry handed in",
    (
        "chemsmart/jobs/xtb/settings.py",
        "XTBJobSettings",
        "JOBTYPES",
        ("hess", "opt", "sp"),
    ): _DECLARATION,
    (
        "chemsmart/settings/capabilities.py",
        "loader_project_section_names",
        None,
        ("gaussian", "orca"),
    ): _GAS_SOLV,
    (
        "chemsmart/agent/bootstrap.py",
        "<module>",
        "_CONFORMANCE_COORDINATES",
        ("modred", "scan"),
    ): "table keyed by the coordinate-driven job types, whose keys the "
    "conformance probe reads from that declaration",
    (
        "chemsmart/agent/live_session.py",
        "<module>",
        "_READER_SUFFIXES",
        ("gaussian", "orca", "pyscf", "xtb"),
    ): "per-program file suffixes, a program-native fact; "
    + _READER_PROGRAMS,
    (
        "chemsmart/agent/live_session.py",
        "<module>",
        "_CONFORMANCE_PROJECT_SHAPES",
        ("gaussian", "orca", "pyscf"),
    ): _DECLARATION
    + ": each loader's own project shape",
    (
        "chemsmart/agent/live_session.py",
        "_pyscf_conformance_sections",
        None,
        ("hess", "irc", "opt", "sp", "td", "ts"),
    ): "table keyed by the PySCF vocabulary; the coverage test below "
    "holds it to PYSCF_JOBTYPES, so a new jobtype cannot lose its "
    "conformance section",
    (
        "chemsmart/agent/program_verifiers.py",
        "<module>",
        "PROGRAM_PREVIEW_VERIFIERS",
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _DECLARATION
    + ": one preview verifier per program",
    (
        "chemsmart/agent/projects.py",
        "<module>",
        "_VOCABULARY_PROVENANCE",
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _DECLARATION
    + ": where each program's vocabulary comes from",
    (
        "chemsmart/agent/workspace_record.py",
        "<module>",
        "_RESULT_SUFFIXES",
        ("gaussian", "orca", "pyscf", "xtb"),
    ): "byte-identical twin of live_session._READER_SUFFIXES; "
    + _READER_PROGRAMS,
    (
        "chemsmart/analysis/result_readers.py",
        "_irc_direction",
        "directions",
        ("irc", "ircf", "ircr"),
    ): _GAUSSIAN_NATIVE,
    (
        "chemsmart/settings/gaussian.py",
        "gaussian_jobtype_settings_classes",
        None,
        ("irc", "link", "qmmm", "td"),
    ): _DECLARATION
    + ": Gaussian's per-jobtype settings classes, now read by the project "
    "reader as well as by the builder, which is why it is a function",
    (
        "chemsmart/settings/orca.py",
        "YamlORCAProjectSettingsBuilder/_project_settings_for_job",
        "settings_mapping",
        ("irc", "neb", "qmmm", "ts"),
    ): _DECLARATION
    + ": ORCA's per-jobtype settings classes",
    (
        "chemsmart/settings/pyscf.py",
        "<module>",
        "PYSCF_STAGE_SOURCES",
        ("hess", "irc", "opt", "sp", "td", "ts"),
    ): "table keyed by the PySCF vocabulary; the coverage test below "
    "holds it to PYSCF_JOBTYPES",
    (
        "chemsmart/io/pyscf/output.py",
        "PySCFOutput/get_molecule",
        "molecule",
        ("opt",),
    ): "derives in C.1 from the stationary-point stage set",
    (
        "chemsmart/io/pyscf/output.py",
        "PySCFOutput/excited_state_record",
        None,
        ("opt", "td"),
    ): "the two stages that can own a followed root: the optimisation "
    "that walks it and the response stage a Hessian differentiates",
    (
        "chemsmart/jobs/pyscf/settings.py",
        "<module>",
        "PYSCF_EXCITED_SURFACE_JOBTYPES",
        ("hess", "opt"),
    ): _DECLARATION
    + ": the job types whose own surface can be an "
    "excited root",
    (
        "chemsmart/analysis/result_readers.py",
        "<module>",
        "RESULT_READERS",
        ("ircf", "ircr"),
    ): _DECLARATION
    + ": the result words one planned Gaussian irc stage writes under, "
    "which is what lets every consumer derive them instead of "
    "spelling the pair out again",
    # R10 Q1 claims: append entries below this line
    # R10 Q1 claims: end
    # R10 Q2 one name, one physics: append entries below this line
    (
        "chemsmart/analysis/result_readers.py",
        "<module>",
        "PRINTED_THERMOCHEMISTRY_CONVENTIONS",
        ("gaussian", "xtb"),
    ): _DECLARATION
    + ": what each program's own printed free energy is, measured by "
    "oracle O1 (R10 Q5); a program missing here is described as its own, "
    "never assumed to be the host's",
    # R10 Q2 one name, one physics: end
    # R10 Q3 knowledge: append entries below this line
    # R10 Q3 knowledge: end
    # R10 Q4 composition: append entries below this line
    (
        "chemsmart/agent/execution.py",
        "<module>",
        "STRUCTURE_HANDOFF_PROGRAMS",
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _DECLARATION
    + ": the programs a structure handoff exists for; the review, the "
    "bounded admission and the frontier each spelled the four out and "
    "now read this one",
    # R10 Q4 composition: end
    # R11 truth: append entries below this line
    # R11 truth: end
    # R11 evidence: append entries below this line
    # R11 evidence: end
    # R11 behaviour: append entries below this line
    # R11 behaviour: end
}


def _mapping_words(node: ast.AST) -> tuple[str, ...] | None:
    """Return the key words of a mapping keyed entirely by the vocabulary."""

    if not isinstance(node, ast.Dict) or not node.keys:
        return None
    if any(key is None for key in node.keys):  # ``**other`` merges in
        return None
    if not all(
        isinstance(key, ast.Constant) and isinstance(key.value, str)
        for key in node.keys
    ):
        return None
    values = [key.value for key in node.keys]
    if all(value in WORDS for value in values):
        return tuple(sorted(set(values)))
    return None


def _literal_words(node: ast.AST) -> tuple[str, ...] | None:
    if isinstance(node, (ast.Set, ast.Tuple, ast.List)):
        elements = node.elts
    elif (
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id in {"frozenset", "set", "tuple"}
        and len(node.args) == 1
        and isinstance(node.args[0], (ast.Set, ast.Tuple, ast.List))
    ):
        elements = node.args[0].elts
    else:
        return None
    if not elements or not all(
        isinstance(item, ast.Constant) and isinstance(item.value, str)
        for item in elements
    ):
        return None
    values = [item.value for item in elements]
    if all(value in WORDS for value in values):
        return tuple(sorted(set(values)))
    return None


def _word_lists(path: pathlib.Path) -> list[tuple[str, str | None, tuple]]:
    found: list[tuple[str, str | None, tuple]] = []

    def visit(node, enclosing, target):
        if isinstance(
            node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        ):
            enclosing = enclosing + [node.name]
        binds = isinstance(node, (ast.Assign, ast.AnnAssign))
        if binds:
            bound = (
                node.targets[0]
                if isinstance(node, ast.Assign)
                else node.target
            )
            if isinstance(bound, ast.Name):
                target = bound.id
            elif isinstance(bound, ast.Attribute):
                target = bound.attr
        words = _literal_words(node) or _mapping_words(node)
        if words is not None:
            found.append(("/".join(enclosing) or "<module>", target, words))
            return
        for child in ast.iter_child_nodes(node):
            # A name reaches the literal it binds even through a wrapper
            # call such as ``MappingProxyType({...})``, so the record can
            # say which table it means; a later statement that binds
            # nothing starts again with no name.
            inherits = binds or not isinstance(node, ast.stmt)
            visit(child, enclosing, target if inherits else None)

    visit(ast.parse(path.read_text()), [], None)
    return found


def _inventory() -> set[tuple[str, str, str | None, tuple[str, ...]]]:
    seen = set()
    for root in ROOTS:
        for path in sorted((REPO / root).rglob("*.py")):
            relative = path.relative_to(REPO).as_posix()
            for enclosing, target, words in _word_lists(path):
                seen.add((relative, enclosing, target, words))
    return seen


def _render(keys) -> str:
    ordered = sorted(
        keys, key=lambda key: (key[0], key[1], key[2] or "", key[3])
    )
    return "\n".join(
        f"  {path}  {enclosing}  {target or '-'}  {list(words)}"
        for path, enclosing, target, words in ordered
    )


def test_every_literal_word_list_is_a_declaration_or_recorded():
    found = _inventory()
    recorded = set(RECORDED)
    unrecorded = found - recorded
    assert not unrecorded, (
        "literal jobtype/program word lists with no author on record. "
        "Derive each from its declaration, or record it above with the "
        "reason it may stand:\n" + _render(unrecorded)
    )


def test_the_record_holds_no_site_the_tree_has_lost():
    stale = set(RECORDED) - _inventory()
    assert not stale, (
        "recorded word lists that are no longer in the tree; delete their "
        "entries so the record keeps shrinking:\n" + _render(stale)
    )


def test_the_record_carries_a_reason_for_every_entry():
    for key, reason in RECORDED.items():
        assert reason.strip(), key


def test_a_table_keyed_by_the_pyscf_vocabulary_covers_all_of_it():
    """A jobtype the vocabulary declares reaches every table keyed by it.

    Asserted on the live objects, never on their source text: round 2's
    loss was a jobtype the reader declared and one hand-list did not
    carry, and a table missing a key is the same loss wearing a mapping.
    """

    from chemsmart.agent.live_session import _pyscf_conformance_sections
    from chemsmart.jobs.pyscf.settings import PYSCF_JOBTYPES
    from chemsmart.settings.pyscf import PYSCF_STAGE_SOURCES

    assert set(PYSCF_STAGE_SOURCES) == set(PYSCF_JOBTYPES)
    assert set(_pyscf_conformance_sections(multiplicity=1)) == set(
        PYSCF_JOBTYPES
    )
    assert set(_pyscf_conformance_sections(multiplicity=2)) == set(
        PYSCF_JOBTYPES
    )


def test_the_stationary_point_declarations_agree_with_each_other():
    """A job type that searches for a structure promises something about
    it, and a job type that samples a surface promises nothing."""

    from chemsmart.agent.terminal_states import (
        GEOMETRY_SEARCH_JOBTYPES,
        STATIONARY_POINT_PROMISES,
        SURFACE_SAMPLING_JOBTYPES,
        expected_imaginary_mode_count,
    )

    assert GEOMETRY_SEARCH_JOBTYPES <= set(STATIONARY_POINT_PROMISES)
    assert not SURFACE_SAMPLING_JOBTYPES & set(STATIONARY_POINT_PROMISES)
    for jobtype in SURFACE_SAMPLING_JOBTYPES:
        assert expected_imaginary_mode_count(jobtype) is None
