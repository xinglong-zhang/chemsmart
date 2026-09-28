"""Task-code audit: what in the tree is named for a task, and can a route reach it?

    PYTHONPATH=<tree> python task_code_audit.py [--scan-root DIR ...]

Reads the tree ``chemsmart`` resolves to (printed first) and never runs an
engine or a provider. Four sections:

1. The expression vocabulary: every operation with its family, whether it
   carries a convention, its admissible input counts and fixed unit, and
   whether its name is a task.
2. For each task-named operation, what its implementation holds (AST): the
   numeric literals in its evaluator branch and in the aggregation function
   it calls, the default values of that function's parameters, and the
   module-level names both read. A protocol constant hidden inside an
   operation shows up here; a physical constant is named.
3. Task job words in the program registry: for every job word that names a
   task rather than a stage, whether the registry declares it, whether it is
   a preview or execution engine-job pair, and what the host says when a
   support overlay tries to admit it (the only door to compile); and the
   capability answer for a run-level task command (``run pka``).
4. A lexicon scan of code under ``chemsmart/agent`` and
   ``chemsmart/analysis``: identifiers (function, class, argument, variable
   and attribute names) and non-docstring string literals that name a
   task, a molecule or a protocol, with file:line. Identifiers are the
   strong signal; string literals are prose or registry data and are
   reported separately so the two are never mixed.
"""

from __future__ import annotations

import argparse
import ast
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

TASK_NAMED_OPERATIONS = (
    "gibbs_to_pka",
    "gibbs_to_redox_potential",
    "exponential_cbs_limit",
    "scf_exponential_cbs_limit",
    "scf_inverse_power_cbs_limit",
    "correlation_inverse_power_cbs_limit",
    "transition_state_crossover_temperature",
)

TASK_JOB_WORDS = (
    "pka",
    "dias",
    "nci",
    "resp",
    "wbi",
    "crest",
    "qrc",
    "traj",
    "neb",
    "link",
    "userjob",
)

#: Words that name a task, a molecule or a protocol. A hit is a place to
#: read, never a verdict.
LEXICON = re.compile(
    r"(pka|pcet|redox|fukui|bde\b|dissociation|proton_affinity|"
    r"square_scheme|crossover|acidity|tautomer|hammett|malonaldehyde|"
    r"phenol|binol|acetamide|benzene|ethylene|formaldehyde|methanol|"
    r"ammonia|losartan|indomethacin|atorvastatin|ferrocene|hcn\b|hnc\b|"
    r"helgaker|karton|martin_cbs|g4mp2|cbs_qb3|w1\b)",
    re.IGNORECASE,
)


def _operation_blocks(tree: ast.AST) -> dict[str, list[ast.If]]:
    """``if operation == "x"`` / ``if operation in {...}`` branches."""

    blocks: dict[str, list[ast.If]] = defaultdict(list)
    for node in ast.walk(tree):
        if not isinstance(node, ast.If) or not isinstance(
            node.test, ast.Compare
        ):
            continue
        test = node.test
        if not (
            isinstance(test.left, ast.Name) and test.left.id == "operation"
        ):
            continue
        comparator = test.comparators[0]
        names: list[str] = []
        if isinstance(test.ops[0], ast.Eq) and isinstance(
            comparator, ast.Constant
        ):
            names = [str(comparator.value)]
        elif isinstance(test.ops[0], ast.In) and isinstance(
            comparator, (ast.Set, ast.Tuple, ast.List)
        ):
            names = [
                str(item.value)
                for item in comparator.elts
                if isinstance(item, ast.Constant)
            ]
        for name in names:
            blocks[name].append(node)
    return blocks


def _literals_and_names(nodes) -> tuple[list, list[str]]:
    literals, names = [], set()
    for root in nodes:
        for node in ast.walk(root):
            if (
                isinstance(node, ast.Constant)
                and isinstance(node.value, (int, float))
                and not isinstance(node.value, bool)
            ):
                literals.append(node.value)
            elif isinstance(node, ast.Name) and node.id.isupper():
                names.add(node.id)
            elif isinstance(node, ast.Attribute) and node.attr.startswith("_"):
                names.add(node.attr)
    return literals, sorted(names)


def _called_aggregation_functions(nodes, available: set[str]) -> list[str]:
    called = set()
    for root in nodes:
        for node in ast.walk(root):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id in available:
                    called.add(node.func.id)
    return sorted(called)


def section_operations(chemsmart_root: Path) -> None:
    from chemsmart.analysis import quantity_expressions as qe

    print("\n== 1. expression vocabulary ==")
    ops = sorted(qe._OPERATIONS)
    fixed = getattr(qe, "_FIXED_DIMENSION_OPERATION_UNITS", {})
    for op in ops:
        print(
            f"  {op:40s} family={qe.OPERATION_FAMILIES.get(op, '?'):12s} "
            f"convention={'yes' if op in qe.CONVENTION_OPERATIONS else 'no ':3s} "
            f"inputs={qe.OPERATION_INPUT_COUNTS[op]!s:14s} "
            f"unit={fixed.get(op, '-'):8s} "
            f"{'TASK-NAMED' if op in TASK_NAMED_OPERATIONS else ''}"
        )
    print(f"  total {len(ops)} operations")

    print("\n== 2. what each task-named operation holds ==")
    source = Path(qe.__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    blocks = _operation_blocks(tree)
    from chemsmart.analysis import aggregation

    agg_source = Path(aggregation.__file__).read_text(encoding="utf-8")
    agg_tree = ast.parse(agg_source)
    agg_functions = {
        node.name: node
        for node in agg_tree.body
        if isinstance(node, ast.FunctionDef)
    }
    agg_constants = {}
    for node in agg_tree.body:
        if isinstance(node, ast.Assign) and isinstance(
            node.targets[0], ast.Name
        ):
            agg_constants[node.targets[0].id] = ast.get_source_segment(
                agg_source, node.value
            )
    for op in TASK_NAMED_OPERATIONS:
        branch = blocks.get(op, [])
        literals, names = _literals_and_names(branch)
        called = _called_aggregation_functions(branch, set(agg_functions))
        lines = [f"{b.lineno}-{b.end_lineno}" for b in branch]
        print(f"  {op}: evaluator lines {lines}")
        print(f"      literals in branch: {literals}")
        print(f"      module constants read: {names}")
        for name in called:
            fn = agg_functions[name]
            fn_literals, fn_names = _literals_and_names([fn])
            defaults = {}
            args = fn.args
            positional = args.args[len(args.args) - len(args.defaults) :]
            for arg, value in zip(positional, args.defaults):
                defaults[arg.arg] = ast.get_source_segment(agg_source, value)
            for arg, value in zip(args.kwonlyargs, args.kw_defaults):
                if value is not None:
                    defaults[arg.arg] = ast.get_source_segment(
                        agg_source, value
                    )
            print(
                f"      calls aggregation.{name} (lines {fn.lineno}-"
                f"{fn.end_lineno}): literals {fn_literals}; defaults "
                f"{defaults}; constants {fn_names}"
            )
            for constant in fn_names:
                if constant in agg_constants:
                    print(f"          {constant} = {agg_constants[constant]}")
        for constant in names:
            if constant in agg_constants:
                print(f"      {constant} = {agg_constants[constant]}")


def section_registry() -> None:
    from chemsmart.agent.capabilities import (
        CapabilityQueryV1,
        ProgramSupportRuleV1,
        SupportLevel,
        build_support_overlay,
        load_program_capabilities,
        query_capability,
    )

    print("\n== 3. task job words in the program registry ==")
    from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

    registry = load_program_capabilities()
    for program, hub in sorted(PROGRAM_CAPABILITIES.items()):
        hub_words = [w for w in TASK_JOB_WORDS if w in hub.jobtypes]
        agent = registry.get(program)
        agent_words = (
            [w for w in TASK_JOB_WORDS if w in agent.jobtypes] if agent else []
        )
        print(
            f"  {program}: hub CLI task words {hub_words}; agent registry "
            f"{'absent' if agent is None else 'task words ' + str(agent_words)}"
        )
    for capability in registry.programs:
        hub = PROGRAM_CAPABILITIES.get(capability.program)
        declared = [
            w
            for w in TASK_JOB_WORDS
            if w in capability.jobtypes or (hub and w in hub.jobtypes)
        ]
        if not declared:
            continue
        preview = {j for _e, j in capability.preview_engine_job_pairs}
        execute = {j for _e, j in capability.execution_engine_job_pairs}
        # A rule carrying well-formed (placeholder) conformance digests, so
        # the probe reaches the overlay's broadening check -- the check that
        # decides whether any conformance evidence could admit the pair.
        placeholder = "0" * 64
        for word in declared:
            try:
                rule = ProgramSupportRuleV1(
                    program=capability.program,
                    support_level=SupportLevel.PREVIEW_ONLY,
                    rule_ids=("c5.audit.probe",),
                    allowed_jobtypes=(word,),
                    allowed_engines=("cpu",),
                    allowed_engine_job_pairs=(("cpu", word),),
                    compiler_evidence_sha256=placeholder,
                    preview_evidence_sha256=placeholder,
                    preflight_evidence_sha256=placeholder,
                )
                build_support_overlay(
                    overlay_id="c5-audit-probe",
                    registry=registry,
                    rules=(rule,),
                )
                probe = "ADMITTED by an overlay"
            except Exception as exc:  # the refusal is the finding
                probe = f"refused: {type(exc).__name__}: {exc}"
            print(
                f"  {capability.program}:{word:8s} declared jobtype; "
                f"preview pair={'yes' if word in preview else 'no'}; "
                f"execution pair={'yes' if word in execute else 'no'}; "
                f"overlay probe {probe}"
            )
    receipt = query_capability(
        CapabilityQueryV1(program="pka", engine="cpu", jobtype="analyze")
    )
    print(f"  run-level 'pka' command as a program: status={receipt.status}")


def _docstring_nodes(tree: ast.AST) -> set[int]:
    ids = set()
    for node in ast.walk(tree):
        if isinstance(
            node,
            (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
        ):
            body = getattr(node, "body", [])
            if (
                body
                and isinstance(body[0], ast.Expr)
                and isinstance(body[0].value, ast.Constant)
                and isinstance(body[0].value.value, str)
            ):
                ids.add(id(body[0].value))
    return ids


def section_lexicon(roots: list[Path], repo: Path) -> None:
    print("\n== 4. lexicon scan of code (identifiers vs string literals) ==")
    ident_hits: dict[str, list[str]] = defaultdict(list)
    string_hits: Counter = Counter()
    string_examples: dict[str, list[str]] = defaultdict(list)
    for root in roots:
        for path in sorted(root.rglob("*.py")):
            try:
                source = path.read_text(encoding="utf-8")
                tree = ast.parse(source)
            except (OSError, SyntaxError, UnicodeDecodeError):
                continue
            docstrings = _docstring_nodes(tree)
            rel = path.relative_to(repo)
            for node in ast.walk(tree):
                name = None
                if isinstance(
                    node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
                ):
                    name = node.name
                elif isinstance(node, ast.arg):
                    name = node.arg
                elif isinstance(node, ast.Name):
                    name = node.id
                elif isinstance(node, ast.Attribute):
                    name = node.attr
                if name and LEXICON.search(name):
                    ident_hits[name].append(f"{rel}:{node.lineno}")
                if (
                    isinstance(node, ast.Constant)
                    and isinstance(node.value, str)
                    and id(node) not in docstrings
                ):
                    match = LEXICON.search(node.value)
                    if match:
                        key = f"{rel}"
                        string_hits[key] += 1
                        if len(string_examples[key]) < 3:
                            string_examples[key].append(
                                f"{node.lineno}: ...{node.value[max(0, match.start() - 40):match.end() + 40]!r}..."
                            )
    print("  identifiers (name: occurrences, first sites):")
    for name in sorted(ident_hits, key=lambda n: -len(ident_hits[n])):
        sites = ident_hits[name]
        print(f"    {name}: {len(sites)}  {sites[:4]}")
    print("  non-docstring string literals (file: count, examples):")
    for key in sorted(string_hits, key=lambda k: -string_hits[k]):
        print(f"    {key}: {string_hits[key]}")
        for example in string_examples[key]:
            print(f"        {example[:200]}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--scan-root", action="append", default=[])
    args = parser.parse_args(argv)
    import chemsmart

    package = Path(chemsmart.__file__).resolve().parent
    repo = package.parent
    print("chemsmart from", chemsmart.__file__)
    section_operations(package)
    section_registry()
    roots = [Path(item).resolve() for item in args.scan_root] or [
        package / "agent",
        package / "analysis",
    ]
    section_lexicon(roots, repo)
    return 0


if __name__ == "__main__":
    sys.exit(main())
