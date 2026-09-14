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
_GUIDE_SIGNAL = "declaration: one guide's activation signal"
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
        "chemsmart/agent/bootstrap.py",
        "bootstrap_program_conformance",
        None,
        ("modred", "scan"),
    ): "derives in 0.4 from commands._COORDINATE_DRIVEN_JOBTYPES",
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
        "chemsmart/agent/guides.py",
        "<module>",
        "GUIDES",
        ("ts",),
    ): _GUIDE_SIGNAL,
    (
        "chemsmart/agent/guides.py",
        "<module>",
        "GUIDES",
        ("scan",),
    ): _GUIDE_SIGNAL,
    (
        "chemsmart/agent/guides.py",
        "<module>",
        "GUIDES",
        ("td",),
    ): _GUIDE_SIGNAL,
    (
        "chemsmart/agent/guides.py",
        "<module>",
        "GUIDES",
        ("pyscf",),
    ): _GUIDE_SIGNAL,
    (
        "chemsmart/agent/guides.py",
        "<module>",
        "GUIDES",
        ("irc", "ts"),
    ): _GUIDE_SIGNAL,
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
        "chemsmart/agent/program_verifiers.py",
        "_validate_gaussian_irc_bundle",
        None,
        ("ircf", "ircr"),
    ): _GAUSSIAN_NATIVE,
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
        "CommandCompiledToolHostV1/_bounded_deferred_target_ids",
        None,
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _READER_PROGRAMS,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/_execute_approved_program_node",
        None,
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _READER_PROGRAMS,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/build_execution_review",
        None,
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _READER_PROGRAMS,
    (
        "chemsmart/agent/tool_runtime.py",
        "CommandCompiledToolHostV1/_admit_bounded_workflow",
        None,
        ("gaussian", "orca", "pyscf", "xtb"),
    ): _READER_PROGRAMS,
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
        ("hess", "opt", "sp", "td"),
    ): _DECLARATION,
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
    ): "table keyed by the coordinate-driven jobtypes; 0.4 keys it off "
    "commands._COORDINATE_DRIVEN_JOBTYPES",
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
        ("hess", "opt", "sp", "td"),
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
        "YamlGaussianProjectSettingsBuilder/_project_settings_for_job",
        "settings_mapping",
        ("irc", "link", "qmmm", "td"),
    ): _DECLARATION
    + ": Gaussian's per-jobtype settings classes",
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
        ("hess", "opt", "sp", "td"),
    ): "table keyed by the PySCF vocabulary; the coverage test below "
    "holds it to PYSCF_JOBTYPES",
    (
        "chemsmart/io/pyscf/output.py",
        "PySCFOutput/get_molecule",
        "molecule",
        ("opt",),
    ): "derives in C.1 from the stationary-point stage set",
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
