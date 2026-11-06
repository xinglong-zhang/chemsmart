#!/usr/bin/env python3
"""
Find duplicated or redefined constants, regex patterns, functions and
methods introduced by a change, for the xinglong-zhang-code-reviewer skill.

Usage (from the repository root)::

    python .agents/skills/xinglong-zhang-code-reviewer/scripts/find_duplicates.py --base main
    python .agents/skills/xinglong-zhang-code-reviewer/scripts/find_duplicates.py --base main --head pr-branch
    python .agents/skills/xinglong-zhang-code-reviewer/scripts/find_duplicates.py --all

``--base REF`` reports definitions on lines changed relative to ``REF``
in the working tree (``git diff REF``, so uncommitted work counts).
Add ``--head REV`` to review a commit or branch that is not checked out:
the code is then read at ``REV`` and compared with its merge base with
``REF`` (``git diff REF...REV``), as a pull request shows it. ``--all``
reports every duplicate in the codebase.

Findings (each with file:line and the definition it duplicates):

- ``redefinition``: a name assigned or defined twice in the same module or
  class body. Python keeps the last definition silently. Ruff's F811
  catches unused function redefinitions; this also covers constants and
  definitions used between the two.
- ``same-constant``: a numeric constant with the same value as one
  defined elsewhere, under any name; ``near-constant`` when the values
  agree to 1e-6 relative but differ, which usually means one quantity
  from two sources (e.g. two CODATA releases). When ``ase`` is installed,
  values in ``ase.units`` (and their reciprocals) are compared too, since
  CHEMSMART derives its unit conversions from them.
- ``regex-duplicate``: a regex already defined elsewhere;
  ``regex-outside-repattern``: a new regex literal not in
  ``chemsmart/utils/repattern.py``.
- ``same-body``: a function or method that does the same as another one
  under any name: identical after removing names, docstrings, and the
  names of parameters and local variables. Contributors often write the
  same routine independently, so names are not compared.
- ``similar-body``: a function or method whose normalised code mostly
  overlaps another one's (Jaccard similarity of 3-token windows >= 0.75):
  a copy with small edits, or a shared routine waiting to be factored out.
- ``same-name``: a new public class whose name already exists in another
  module (informational; same-named functions in different modules are
  normal).

The script only uses the standard library and parses files without
importing them. It reports candidates; check each before raising it.
"""

import argparse
import ast
import hashlib
import math
import re
import subprocess
import sys
import warnings
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

REPATTERN_MODULE = "chemsmart/utils/repattern.py"
SOURCE_ROOTS = ("chemsmart",)
REGEX_FUNCTIONS = {
    "compile",
    "search",
    "match",
    "fullmatch",
    "findall",
    "finditer",
    "sub",
    "subn",
    "split",
}
# Numbers too common to be meaningful duplicates
MIN_SIGNIFICANT_DIGITS = 4
MIN_BODY_STATEMENTS = 3
MIN_BODY_TOKENS = 25  # ignore tiny functions for similarity
SIMILARITY_THRESHOLD = 0.75


@dataclass(frozen=True)
class Definition:
    kind: str  # constant, regex, function, class
    name: str
    path: str
    line: int
    end_line: int
    scope: str  # "module" or the enclosing class name
    key: object = None  # value, pattern or body hash
    accessor: bool = False  # @x.setter/.deleter/.getter or @overload
    shingles: frozenset = frozenset()  # 3-token windows of normalised code
    direct: bool = True  # a direct statement of the module/class body

    @property
    def where(self):
        """``path:line``, or just the source for reference values."""
        return f"{self.path}:{self.line}" if self.line else self.path


def significant_digits(value):
    text = repr(abs(value)).lower().split("e")[0].replace(".", "")
    return len(text.strip("0"))


class _LocalNames(ast.NodeTransformer):
    """Rename parameters and local variables to v0, v1, ... in order of
    first appearance, so renamed copies compare equal. Attribute names
    (``self.natoms``), globals, imports and builtins keep their names."""

    def __init__(self, local_names):
        self.local_names = local_names
        self.mapping = {}

    def _placeholder(self, name):
        if name not in self.mapping:
            self.mapping[name] = f"v{len(self.mapping)}"
        return self.mapping[name]

    def visit_Name(self, node):
        if node.id in self.local_names:
            node.id = self._placeholder(node.id)
        return node

    def visit_arg(self, node):
        node.arg = self._placeholder(node.arg)
        node.annotation = None
        return node


def _local_names(node):
    names = {a.arg for a in ast.walk(node.args) if isinstance(a, ast.arg)}
    for child in ast.walk(node):
        if isinstance(child, ast.Name) and isinstance(child.ctx, ast.Store):
            names.add(child.id)
    return names


def _normalised_body(node):
    """The function's parameters and body without docstring, annotations,
    or local names; ``None`` if too short to be meaningful."""
    body = list(node.body)
    if (
        body
        and isinstance(body[0], ast.Expr)
        and isinstance(body[0].value, ast.Constant)
        and isinstance(body[0].value.value, str)
    ):
        body = body[1:]
    if len(body) < MIN_BODY_STATEMENTS:
        return None
    clone = ast.parse(ast.unparse(node)).body[0]  # independent copy
    clone.body = clone.body[len(node.body) - len(body) :]
    clone.returns = None
    renamer = _LocalNames(_local_names(clone))
    renamer.visit(clone.args)
    for statement in clone.body:
        renamer.visit(statement)
    return clone


_PLACEHOLDER = re.compile(r"v\d+")


def _tokens(node):
    """Node-type tokens for similarity. All local-name placeholders become
    one token, so an extra local variable does not shift every later one."""
    tokens = []
    for child in ast.walk(node):
        token = type(child).__name__
        for field in ("id", "attr", "arg"):
            value = getattr(child, field, None)
            if isinstance(value, str):
                if field != "attr" and _PLACEHOLDER.fullmatch(value):
                    value = "v"
                token += f":{value}"
        if isinstance(child, ast.Constant):
            token += f":{child.value!r}"[:40]
        tokens.append(token)
    return tokens


def body_signature(node):
    """(hash, 3-token shingles) of a normalised function, or (None, ())."""
    clone = _normalised_body(node)
    if clone is None:
        return None, frozenset()
    dumped = ast.dump(clone.args, annotate_fields=False) + "".join(
        ast.dump(statement, annotate_fields=False) for statement in clone.body
    )
    tokens = [t for statement in clone.body for t in _tokens(statement)]
    shingles = (
        frozenset(zip(tokens, tokens[1:], tokens[2:]))
        if len(tokens) >= MIN_BODY_TOKENS
        else frozenset()
    )
    return hashlib.sha1(dumped.encode()).hexdigest(), shingles


class Collector(ast.NodeVisitor):
    """Collect definitions from one module."""

    def __init__(self, path):
        self.path = path
        self.scope = ["module"]
        self.definitions = []
        self._direct = set()  # ids of direct module/class body statements

    def visit_Module(self, node):
        self._direct.update(id(statement) for statement in node.body)
        self.generic_visit(node)

    def _add(
        self, kind, name, node, key=None, accessor=False, shingles=frozenset()
    ):
        self.definitions.append(
            Definition(
                kind,
                name,
                self.path,
                node.lineno,
                getattr(node, "end_lineno", node.lineno),
                self.scope[-1],
                key,
                accessor,
                shingles,
                id(node) in self._direct,
            )
        )

    def visit_ClassDef(self, node):
        if len(self.scope) == 1:
            self._add("class", node.name, node)
        self._direct.update(id(statement) for statement in node.body)
        self.scope.append(node.name)
        self.generic_visit(node)
        self.scope.pop()

    def visit_FunctionDef(self, node):
        if len(self.scope) <= 2:  # module functions and methods
            digest, shingles = body_signature(node)
            self._add(
                "function",
                node.name,
                node,
                digest,
                accessor=_is_accessor(node),
                shingles=shingles,
            )
        self.scope.append(f"{self.scope[-1]}.{node.name}()")
        self.generic_visit(node)
        self.scope.pop()

    visit_AsyncFunctionDef = visit_FunctionDef

    def _constant(self, target_names, value_node, node):
        if len(self.scope) > 2 or self.scope[-1].endswith("()"):
            return  # local variables are not shared definitions
        value = _literal_number(value_node)
        for name in target_names:
            if value is not None:
                self._add("constant", name, node, value)
            elif _is_pattern_name(name) and _literal_string(value_node):
                self._add("regex", name, node, _literal_string(value_node))

    def visit_Assign(self, node):
        names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        self._constant(names, node.value, node)
        self.generic_visit(node)

    def visit_AnnAssign(self, node):
        if isinstance(node.target, ast.Name) and node.value is not None:
            self._constant([node.target.id], node.value, node)
        self.generic_visit(node)

    def visit_Call(self, node):
        func = node.func
        if (
            isinstance(func, ast.Attribute)
            and func.attr in REGEX_FUNCTIONS
            and isinstance(func.value, ast.Name)
            and func.value.id in {"re", "regex"}
            and node.args
        ):
            pattern = _literal_string(node.args[0])
            if pattern is not None:
                self._add("regex", f"re.{func.attr}(...)", node, pattern)
        self.generic_visit(node)


def _is_accessor(node):
    """Decorators that legitimately reuse a name in the same scope."""
    for decorator in node.decorator_list:
        if isinstance(decorator, ast.Attribute) and decorator.attr in {
            "setter",
            "deleter",
            "getter",
        }:
            return True
        name = getattr(decorator, "id", getattr(decorator, "attr", None))
        if name == "overload":
            return True
    return False


def _literal_number(node):
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        inner = _literal_number(node.operand)
        return None if inner is None else -inner
    if (
        isinstance(node, ast.Constant)
        and isinstance(node.value, (int, float))
        and not isinstance(node.value, bool)
    ):
        value = float(node.value)
        if math.isfinite(value):
            return value
    return None


def _distinctive(value):
    """Numbers specific enough that equal values suggest one quantity."""
    return significant_digits(value) >= MIN_SIGNIFICANT_DIGITS


def _literal_string(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _is_pattern_name(name):
    lowered = name.lower()
    return lowered.endswith(("_pattern", "_regex", "_re")) or lowered in {
        "pattern",
        "regex",
    }


def _git(root, *args):
    return subprocess.run(
        ["git", *args], cwd=root, capture_output=True, text=True, check=True
    ).stdout


def _sources(root, rev=None):
    """Yield (relative path, source text) for Python files under
    SOURCE_ROOTS, from the working tree or from commit ``rev``."""
    if rev is None:
        for source_root in SOURCE_ROOTS:
            for path in sorted((root / source_root).rglob("*.py")):
                yield path.relative_to(root).as_posix(), path.read_text(
                    encoding="utf-8"
                )
        return
    listed = _git(
        root, "ls-tree", "-r", "--name-only", rev, "--", *SOURCE_ROOTS
    )
    for relative in sorted(listed.split()):
        if relative.endswith(".py"):
            yield relative, _git(root, "show", f"{rev}:{relative}")


def collect(root, rev=None):
    definitions = []
    for relative, text in _sources(root, rev):
        try:
            with warnings.catch_warnings():
                # e.g. invalid escape sequences in existing strings
                warnings.simplefilter("ignore", SyntaxWarning)
                tree = ast.parse(text)
        except (SyntaxError, UnicodeDecodeError) as exc:
            print(f"skipped {relative}: {exc}", file=sys.stderr)
            continue
        collector = Collector(relative)
        collector.visit(tree)
        definitions += collector.definitions
    return definitions


def reference_constants():
    """``ase.units`` values and reciprocals, if ``ase`` is importable."""
    try:
        from ase import units
    except ImportError:
        return []
    references = []
    for name in dir(units):
        value = getattr(units, name)
        if (
            name.startswith("_")
            or isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value == 0
        ):
            continue
        for label, number in ((name, float(value)), (f"1/{name}", 1 / value)):
            if _distinctive(number):
                references.append(
                    Definition(
                        "constant",
                        f"ase.units.{label}",
                        "ase.units",
                        0,
                        0,
                        "module",
                        number,
                    )
                )
    return references


def changed_lines(root, base, head=None):
    """Map of path -> set of added/changed line numbers: the working tree
    versus ``base``, or ``head`` versus its merge base with ``base``."""
    target = base if head is None else f"{base}...{head}"
    diff = _git(root, "diff", "-U0", target, "--", *SOURCE_ROOTS)
    lines = defaultdict(set)
    path = None
    for line in diff.splitlines():
        if line.startswith("+++ "):
            path = line[6:] if line.startswith("+++ b/") else None
        elif line.startswith("@@") and path:
            match = re.search(r"\+(\d+)(?:,(\d+))?", line)
            start, count = int(match[1]), int(match[2] or 1)
            lines[path].update(range(start, start + count))
    return lines


def touches(definition, changed):
    rows = changed.get(definition.path)
    if not rows:
        return False
    if definition.kind in {"function", "class"}:
        # a new or rewritten def/class line, not any edit in its body
        return definition.line in rows
    return any(
        row in rows for row in range(definition.line, definition.end_line + 1)
    )


def find(definitions, selected, threshold=SIMILARITY_THRESHOLD):
    findings = []
    by_scope = defaultdict(list)
    for d in definitions:
        if d.kind in {"constant", "function", "class"} or (
            d.kind == "regex" and not d.name.startswith("re.")
        ):
            by_scope[(d.path, d.scope, d.name)].append(d)

    for d in selected:
        if d.kind == "regex" and d.name.startswith("re."):
            pass
        else:
            same = by_scope[(d.path, d.scope, d.name)]
            earlier = [o for o in same if o.line < d.line and o.direct]
            if d.direct and earlier and not _overload(d, earlier[-1]):
                findings.append(
                    (
                        "redefinition",
                        d,
                        earlier[-1],
                        "same scope; Python keeps the last definition",
                    )
                )

        others = [
            o
            for o in definitions
            if o is not d and o.kind == d.kind and o.key is not None
        ]
        if d.kind == "constant" and d.key is not None and _distinctive(d.key):
            for o in others:
                if o.path == d.path and o.name == d.name:
                    continue
                if o.key == d.key:
                    findings.append(
                        ("same-constant", d, o, f"value {d.key!r}")
                    )
                elif math.isclose(o.key, d.key, rel_tol=1e-6):
                    findings.append(
                        (
                            "near-constant",
                            d,
                            o,
                            f"{d.key!r} vs {o.key!r}: one quantity, "
                            "two sources?",
                        )
                    )
        elif d.kind == "regex":
            duplicates = [o for o in others if o.key == d.key]
            for o in duplicates:
                findings.append(("regex-duplicate", d, o, repr(d.key)))
            if d.path != REPATTERN_MODULE and not duplicates:
                findings.append(
                    (
                        "regex-outside-repattern",
                        d,
                        None,
                        f"{d.key!r}; move to {REPATTERN_MODULE}",
                    )
                )
        elif d.kind == "function" and d.key is not None:
            for o in others:
                if (o.path, o.line) == (d.path, d.line):
                    continue
                if o.key == d.key:
                    note = (
                        "identical body"
                        if o.name == d.name
                        else "identical body under another name"
                    )
                    findings.append(("same-body", d, o, note))
                elif d.shingles and o.shingles:
                    score = _similarity(d.shingles, o.shingles, threshold)
                    if score >= threshold:
                        findings.append(
                            ("similar-body", d, o, f"similarity {score:.2f}")
                        )

        if (
            d.kind == "class"
            and d.scope == "module"
            and not d.name.startswith("_")
        ):
            for o in definitions:
                if (
                    o.kind == d.kind
                    and o.scope == "module"
                    and o.name == d.name
                    and o.path != d.path
                ):
                    findings.append(("same-name", d, o, "same public name"))
    return _unique(findings)


def _similarity(first, second, threshold=SIMILARITY_THRESHOLD):
    """Jaccard similarity; 0 for pairs too different in size to reach
    ``threshold`` (Jaccard <= smaller/larger set size)."""
    small, large = sorted((len(first), len(second)))
    if small < threshold * large:
        return 0.0
    return len(first & second) / len(first | second)


def _report(findings):
    """Report lines; functions with identical bodies print as one group."""
    groups = defaultdict(dict)
    lines = []
    for kind, d, o, note in findings:
        if kind == "same-body":
            for member in (d, o):
                groups[d.key][(member.path, member.line)] = member
            if len(groups[d.key]) == 2 and d.key not in {
                key for key, _ in lines if key is not None
            }:
                lines.append((d.key, None))
            continue
        other = f" duplicates {o.name} at {o.where}" if o else ""
        lines.append((None, f"{kind}: {d.name} at {d.where}{other} ({note})"))
    for key, text in lines:
        if key is None:
            yield text
        else:
            members = sorted(groups[key].values(), key=lambda m: m.where)
            places = ", ".join(f"{m.name} at {m.where}" for m in members)
            yield (
                f"same-body: {len(members)} functions with identical bodies "
                f"(names and local variables ignored): {places}"
            )


def _overload(definition, earlier):
    """A property setter/deleter or ``@overload`` reusing a name."""
    return definition.kind == "function" and (
        definition.accessor or earlier.accessor
    )


def _unique(findings):
    seen = set()
    result = []
    for kind, d, o, note in findings:
        pair = (kind,) + tuple(
            sorted([(d.path, d.line), (o.path, o.line) if o else ("", 0)])
        )
        if pair in seen:
            continue
        seen.add(pair)
        result.append((kind, d, o, note))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--base", help="git ref to diff against")
    group.add_argument("--all", action="store_true", help="whole codebase")
    parser.add_argument(
        "--head", help="commit or branch to review (default: working tree)"
    )
    parser.add_argument("--root", default=".", help="repository root")
    parser.add_argument(
        "--similarity",
        type=float,
        default=SIMILARITY_THRESHOLD,
        help="minimum similarity for similar-body (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    if args.head and not args.base:
        parser.error("--head needs --base")

    root = Path(args.root).resolve()
    definitions = collect(root, args.head)
    if args.all:
        selected = definitions
    else:
        changed = changed_lines(root, args.base, args.head)
        selected = [d for d in definitions if touches(d, changed)]

    findings = find(
        definitions + reference_constants(), selected, args.similarity
    )
    order = [
        "redefinition",
        "same-constant",
        "near-constant",
        "regex-duplicate",
        "regex-outside-repattern",
        "same-body",
        "similar-body",
        "same-name",
    ]
    findings.sort(key=lambda f: (order.index(f[0]), f[1].path, f[1].line))
    for line in _report(findings):
        print(line)
    print(
        f"{len(findings)} candidate(s) in {len(selected)} definition(s) "
        "checked",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
