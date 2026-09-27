"""Census E: where the Agent's scientific work leaves the typed evidence layer.

Reads, never writes, under the named roots, and imports nothing from
chemsmart, so it runs on any Python 3.8+ (a cluster slot job, or a laptop
reading a mirror). It is descriptive: it counts and locates, with
denominators, and grades nothing. Whether an exit changed a delivered
conclusion is decided by reading the transcript behind it, never here.

Sources, per root (``os.walk`` following directory links, each real
directory once; excluded globs and private stores pruned before descent):

- ``public-transcript-*.json``: every assistant tool call and its reply,
  every user (task or wake) message, every assistant message. One session
  is one ``transcript_sha256``; a copy of a transcript counts once.
- ``events.jsonl`` beside a transcript: the observed model
  (``provider_turn_observed``), the provider-turn count and the ending
  (``runtime_terminated``) -- a session with no provider turn or one that
  ended on ``turn_deadline_exceeded`` is infrastructure and is kept out of
  every behavioural denominator.
- ``.chemsmart-agent/goals/*/ledger.jsonl``: the settlement of every goal
  of the workspace a session belongs to, and the goal's executed route
  from its run streams (``goals/*/runs/*/events.jsonl``).
- ``.chemsmart-agent/runs/*/projects/*.yaml``: promoted project files
  (R10 Q28's hatch census, adopted unchanged).

Detectors (rows in ``rows.jsonl``; one row per detected event):

- ``hatch``, ``unknown_key``, ``project_file``: R10 Q28's census of
  native escape hatches in Agent-authored project sections, adopted
  verbatim (``HATCH_KEYS``, ``FREE_WORD_KEYS``, ``UNKNOWN_KEY``).
- ``string_exit``: a native-input, path or shell string in any string
  leaf of any tool call, with the leaf's key path (prose fields and
  operative fields are reported apart) and the reply's status.
- ``refusal``: every rejected tool reply, classified by message into the
  vocabulary classes (selector, operation, unit, thermochemistry kind,
  shape) or the reuse classes (unknown receipt, workflow, run, artifact),
  with the requested selector/operation names parsed out. With
  ``--declarations`` (a JSON table of a tree's reader declarations) each
  requested selector is placed: declared for this program and job type,
  for another job type of this program, for another program, or nowhere.
- ``declared_gap``: an analysis node the Agent planned
  ``blocked_unsupported`` (with its reason) and every
  ``plan_unsupported_external`` call: the Agent's own statement that the
  vocabulary lacks something.
- ``citation``: a receipt digest or artifact id a tool call cites that no
  earlier reply of the same session carried -- shown only by the task or
  wake message (``context``) or by nothing (``unknown``) -- with the reply.
- ``number``: a unit-bearing number (at least two significant digits) in
  the Agent's prose, classified by the first rule that matches: ``bound``
  (a value some typed reply of the session carried, at the displayed
  precision), ``context`` (a value only the task or wake message carried),
  ``own_argument`` (a number the Agent itself put in an earlier tool
  argument), ``prose_computed`` (one unit factor applied to a seen value,
  or a sum or difference of two seen values, optionally converted), or
  ``unmatched``. Surfaces: ``finding`` and ``unreachable`` statements and
  the session's last assistant message are *delivered* prose; decision
  text fields and mid-session assistant text are *exploration*. The same
  matcher run against another session's seen values (a fixed derangement)
  is the control: an accidental match rate for every class.
- ``route``: per session the accepted plans' (program, jobtype) nodes,
  analysis kinds and input edges; per goal the executed (program, jobtype)
  and the geometry handoffs and data edges its run streams recorded.

Outputs: ``OUT_DIR/rows.jsonl``, ``OUT_DIR/sessions.jsonl`` and
``OUT_DIR/summary.json``. Rows carry paths and excerpts of private
transcripts: they stay beside the records, never in Git.

    python evidence_census.py OUT_DIR ROOT [ROOT ...]
        [--exclude ROOT:GLOB ...] [--declarations FILE.json]
"""

from __future__ import annotations

import bisect
import collections
import fnmatch
import hashlib
import json
import os
import re
import sys

# --------------------------------------------------------------------------
# R10 Q28's census (instruments/q28/census/census.py, sha256 423e919b...),
# adopted verbatim: the hatch keys, the free-word keys, the loader's
# unknown-key refusal and the project-authoring tools of both generations.
# --------------------------------------------------------------------------

HATCH_KEYS = (
    "input_string",
    "route_to_be_written",
    "additional_route_parameters",
    "additional_opt_options_in_route",
    "additional_opt_options",
    "append_additional_info",
    "additional_solvent_options",
    "solvent_options",
    "custom_solvent",
    "dieze_tag",
    "link_route",
)
#: Typed names whose values are free words written into native input.
FREE_WORD_KEYS = ("scf_algorithm", "scf_tol", "guess", "stable")
ALL_KEYS = HATCH_KEYS + FREE_WORD_KEYS

#: The loader's refusal of a key no settings class owns (utils.update_dict_
#: with_existing_keys): the model asked for something with no typed field.
UNKNOWN_KEY = re.compile(r"Keyword `([^`]+)` is not in list of keywords")

PROJECT_TOOLS = (
    "project_yaml",
    "render_project_yaml",
    "establish_project",
)

try:
    import yaml  # noqa: F401

    def load_yaml(text):
        return yaml.safe_load(text)

except ImportError:  # pragma: no cover - the cluster env has PyYAML
    load_yaml = None


def nonempty(value):
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (list, tuple, dict)):
        return bool(value)
    return value is not False


# --------------------------------------------------------------------------
# Walking
# --------------------------------------------------------------------------

#: Private stores are pruned before descent whatever the roots say.
PRUNE_PREFIX = ("sealed", "grade-")
PRUNE_EXACT = frozenset({"seals", "grading", "claude", ".codex", "scratch"})


def walk(root, excludes):
    """(dirpath, filenames) under root: links followed, each real dir once."""

    seen = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:
            dirnames[:] = []
            continue
        seen.add(real)
        rel = os.path.relpath(dirpath, root)
        keep = []
        for name in dirnames:
            if name in PRUNE_EXACT or name.startswith(PRUNE_PREFIX):
                continue
            child = os.path.normpath(os.path.join(rel, name))
            if any(fnmatch.fnmatch(child, pattern) for pattern in excludes):
                continue
            keep.append(name)
        dirnames[:] = keep
        yield dirpath, filenames


def sha(text):
    return hashlib.sha256(text.encode("utf-8", "replace")).hexdigest()


def load_jsonl(path):
    rows = []
    try:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    try:
                        rows.append(json.loads(line))
                    except ValueError:
                        continue
    except OSError:
        pass
    return rows


def parse_json(text):
    if not isinstance(text, str):
        return text
    try:
        return json.loads(text)
    except ValueError:
        return None


def leaves(value, path=()):
    """(key path, leaf) for every scalar leaf of a JSON value."""

    if isinstance(value, dict):
        for key, item in value.items():
            yield from leaves(item, path + (str(key),))
    elif isinstance(value, list):
        for item in value:
            yield from leaves(item, path + ("[]",))
    else:
        yield path, value


# --------------------------------------------------------------------------
# X: native input, paths and shell in any tool argument
# --------------------------------------------------------------------------

NATIVE_PATTERNS = (
    ("orca_simple_line", re.compile(r"(?m)^\s*!\s*[A-Za-z]")),
    (
        "orca_block",
        re.compile(
            r"(?i)%(pal|maxcore|scf|geom|tddft|cpcm|method|basis|freq|output|"
            r"elprop|mdci|casscf|irc|neb|plots|eprnmr|rel|coords|mp2|loc)\b"
            r"[^%]*?\bend\b",
            re.S,
        ),
    ),
    ("gaussian_route", re.compile(r"(?m)^\s*#[pPnNtT]?\s+[A-Za-z]")),
    (
        "gaussian_link0",
        re.compile(r"(?i)%(chk|oldchk|mem|nprocshared|nproc|rwf)\s*="),
    ),
    (
        "dollar_block",
        re.compile(
            r"(?m)^\s*\$(constrain|fix|scan|opt|chrg|spin|wall|end|md)\b"
        ),
    ),
    (
        "python_code",
        re.compile(
            r"from pyscf|import pyscf|gto\.M\(|scf\.(RHF|UHF|ROHF)\(|"
            r"dft\.(RKS|UKS)\(|\.kernel\(\)|import numpy|subprocess\."
        ),
    ),
)
PATH_PATTERNS = (
    (
        "absolute_path",
        re.compile(
            r"(?:^|[\s\"'=:(,])(/home/|/project/|/Users/|/tmp/|/scratch/|"
            r"/lustre/|/opt/|/usr/|/private/|/var/|~/)\S*"
        ),
    ),
    ("relative_path", re.compile(r"(?:^|[\s\"'=:(,])\.\.?/\S+")),
    (
        "native_file",
        re.compile(
            r"\b[\w.+-]+\.(log|out|inp|com|gjf|chk|fchk|gbw|hess|h5|engrad|"
            r"trj|molden|cube|densities)\b"
        ),
    ),
    ("xyz_file", re.compile(r"\b[\w.+-]+\.xyz\b")),
)
SHELL_PATTERNS = (
    (
        "command_line",
        re.compile(
            r"(?m)^\s*(\$\s+)?(bash|sh|zsh|cat|grep|sed|awk|python3?|chemsmart|"
            r"sbatch|srun|module|conda|pip|ls|cd|cp|mv|rm|orca|g16|xtb)\s+[-./\w]"
        ),
    ),
    ("pipe", re.compile(r"\|\s*(grep|head|tail|awk|sed|sort|wc|cut|tr)\b")),
    ("and_chain", re.compile(r"&&\s*[a-z]")),
    ("subshell", re.compile(r"\$\([^)]+\)|`[^`\n]+`")),
    ("cli_word", re.compile(r"\bchemsmart\s+(run|sub)\b")),
)

#: Leaf keys whose value is the Agent's own prose (a mention there is
#: reported apart from an operative argument).
PROSE_KEYS = frozenset(
    {
        "method_rationale",
        "statement",
        "reason",
        "assumptions",
        "alternatives",
        "uncertainties",
        "diagnostics",
        "rationale",
        "meaning",
        "expectation_basis",
        "blocked_reason",
        "note",
        "notes",
        "description",
        "question",
        "hypothesis",
        "approach",
        "outcome",
        "text",
        "summary",
        "justification",
        "basis",
        "mechanism",
        "failure_update_rule",
        "success_criterion",
        "observation",
        "interpretation",
        "caveat",
        "limitations",
        "purpose",
        "why",
        # A search_capabilities query is words the host indexes, never an
        # argument anything executes; the first CUHK run tagged twelve of
        # them ("orca project yaml ...", "from pyscf hdf5") as operative.
        "query",
    }
)


def string_exits(arguments):
    """Every native/path/shell match in the string leaves of arguments."""

    found = []
    for path, value in leaves(arguments):
        if not isinstance(value, str) or len(value) < 2:
            continue
        key = next((p for p in reversed(path) if p != "[]"), "")
        prose = key in PROSE_KEYS
        for family, patterns in (
            ("native", NATIVE_PATTERNS),
            ("path", PATH_PATTERNS),
            ("shell", SHELL_PATTERNS),
        ):
            for name, pattern in patterns:
                match = pattern.search(value)
                if match:
                    start = max(0, match.start() - 60)
                    found.append(
                        {
                            "family": family,
                            "pattern": name,
                            "leaf": ".".join(path),
                            "prose_field": prose,
                            "excerpt": value[start : match.end() + 80],
                        }
                    )
    return found


# --------------------------------------------------------------------------
# Refusals: vocabulary and reuse classes
# --------------------------------------------------------------------------


def _names(text):
    return tuple(re.findall(r"'([^']+)'", text or ""))


REFUSAL_RULES = (
    # (class, family, regex); first match wins; groups parsed below.
    (
        "selector_undeclared_for_jobtype",
        "vocabulary",
        re.compile(
            r"selector\(s\) \[(?P<sel>[^\]]*)\] (?:are|is) not declared for "
            r"'?(?P<program>\w+)'? jobtype '(?P<jobtype>[\w-]+)'"
        ),
    ),
    (
        "selector_undeclared_for_jobtype",
        "vocabulary",
        re.compile(
            r"requests selector\(s\) \[(?P<sel>[^\]]*)\] that are not declared "
            r"for '(?P<program>\w+)' jobtype '(?P<jobtype>[\w-]+)'"
        ),
    ),
    (
        "selector_not_in_enum",
        "vocabulary",
        re.compile(r"\.selector is '(?P<sel1>[^']+)', which is not one of"),
    ),
    (
        "reader_does_not_provide",
        "vocabulary",
        re.compile(
            r"(?P<program>\w+) result reader does not provide '(?P<sel1>[^']+)'"
        ),
    ),
    (
        "value_absent_in_result",
        "vocabulary",
        re.compile(
            r"(?P<program>\w+) result contains no '(?P<sel1>[^']+)' value"
        ),
    ),
    (
        "value_absent_in_result",
        "vocabulary",
        re.compile(r"prints no <S\^2> expectation value"),
    ),
    (
        "operation_not_in_enum",
        "vocabulary",
        re.compile(r"\.operation is '(?P<op>[^']+)', which is not one of"),
    ),
    (
        "operation_unknown",
        "vocabulary",
        re.compile(
            r"(?i)unknown (?:expression )?operation '?(?P<op>[\w-]+)'?"
        ),
    ),
    (
        "field_not_accepted",
        "schema",
        re.compile(
            r"supplied \[(?P<fields>[^\]]*)\], which this object does not accept"
        ),
    ),
    (
        "unit_unsupported",
        "vocabulary",
        re.compile(r"unsupported unit: '(?P<unit>[^']*)'"),
    ),
    (
        "unit_unsupported",
        "vocabulary",
        re.compile(r"declares unsupported unit '(?P<unit>[^']*)'"),
    ),
    (
        "unit_unreachable",
        "vocabulary",
        re.compile(r"convert to '(?P<unit>[^']*)' is unreachable"),
    ),
    (
        "thermochemistry_kind_not_derived",
        "vocabulary",
        re.compile(
            r"declares quantity_kind '(?P<kind>[^']+)', which thermochemistry "
            r"does not derive"
        ),
    ),
    (
        "shape",
        "vocabulary",
        re.compile(
            r"requires one (matrix|vector)|is not numerical|requires the values|"
            r"accepts only input_ids|got \d+ inputs|requires exactly|"
            r"repeated source quantities"
        ),
    ),
    (
        "receipt_unknown",
        "reuse",
        re.compile(
            r"references an unknown receipt|cites an unknown receipt|"
            r"cites an unknown postprocessing receipt|"
            r"receipt_is_one_the_host_minted|is not in this \w+ receipt"
        ),
    ),
    (
        "workflow_unknown",
        "reuse",
        re.compile(
            r"unknown scientific workflow ID|workflow\.id_names_a_plan|"
            r"plan\.amend_needs_a_recorded_plan"
        ),
    ),
    (
        "run_unknown",
        "reuse",
        re.compile(r"records no run '|run_reference_names_one_recorded_run"),
    ),
    (
        "artifact_unknown",
        "reuse",
        re.compile(
            r"unknown trusted artifact ID|artifact\.id_is_registered|"
            r"artifact_id is '[^']*', which does not match"
        ),
    ),
    (
        "unknown_key",
        "native",
        UNKNOWN_KEY,
    ),
    (
        "native_word_refused",
        "native",
        re.compile(r"native_words_have_typed_settings|native word"),
    ),
    (
        "cli_option_absent",
        "native",
        re.compile(
            r"live Click scope (?P<scope>\S+) has no (?P<option>\w+) option"
        ),
    ),
)


def classify_refusal(message):
    for name, family, rule in REFUSAL_RULES:
        match = rule.search(message or "")
        if not match:
            continue
        groups = {k: v for k, v in match.groupdict().items() if v}
        requested = []
        if "sel" in groups:
            requested.extend(_names(groups["sel"]))
        if "sel1" in groups:
            requested.append(groups["sel1"])
        if name == "value_absent_in_result" and "sel1" not in groups:
            requested.append("spin_square")
        declared_here = ()
        tail = (message or "").split("Declared here:", 1)
        if len(tail) == 2:
            declared_here = _names(tail[1].split("]", 1)[0])
        return {
            "class": name,
            "family": family,
            "requested": tuple(requested),
            "program": groups.get("program"),
            "jobtype": groups.get("jobtype"),
            "operation": groups.get("op"),
            "unit": groups.get("unit"),
            "kind": groups.get("kind"),
            "declared_here": declared_here,
        }
    return {"class": "other", "family": "other", "requested": ()}


def place_selector(selector, program, jobtype, declarations):
    """Where a tree's readers declare a selector, relative to the request."""

    if not declarations:
        return None
    programs = declarations.get("programs") or {}
    entry = programs.get(str(program or "")) or {}
    jobtypes = entry.get("jobtypes") or {}
    if jobtype and selector in (jobtypes.get(jobtype) or ()):
        return "declared_here_today"
    if selector in (entry.get("selectors") or ()):
        return "declared_other_jobtype"
    if any(
        selector in (other.get("selectors") or ())
        for name, other in programs.items()
        if name != program
    ):
        return "declared_other_program"
    return "declared_nowhere"


# --------------------------------------------------------------------------
# Numbers
# --------------------------------------------------------------------------

#: A number as prose writes it: an ASCII or Unicode minus, thousands
#: grouped by commas or thin spaces (-147,707.15 is one number, not the
#: fragment 707.15), and an exponent that may carry a Unicode minus.
NUMBER = re.compile(
    r"(?<![\w.,  ])[-−+]?(?:\d{1,3}(?:[,  ]\d{3})+|\d+)"
    r"(?:\.\d+)?(?:[eE][-−+]?\d+)?(?!\w)(?!\.\d)(?![,  ]\d)"
)
_SEPARATORS = str.maketrans({",": None, " ": None, " ": None, "−": "-"})


def clean_number(text):
    return text.translate(_SEPARATORS)


UNIT_AFTER = re.compile(
    r"\s*(?:±\s*[\d.]+\s*)?(?P<unit>kcal\s*/?\s*mol|kcal\s*mol|kJ\s*/?\s*mol|"
    r"kJ\s*mol|meV|eV|mEh|mHa|cm\s*(?:-1|⁻¹|\^-1|\^\{-1\})|"
    r"Å|[Aa]ngstr(?:o|ö)m|E_?h\b|[Hh]artrees?|a\.u\.|nm\b|ppm|"
    r"[Dd]ebye|°|deg(?:rees)?\b|GHz|MHz)"
)
#: Conversion factors a chemist applies by hand (target = factor * value).
FACTORS = (
    627.509474,
    2625.49964,
    27.211386,
    219474.631,
    4.184,
    23.060548,
    96.485332,
    8065.544,
    0.01196266,
    0.00285914,
    0.529177,
    1000.0,
)
FACTORS = FACTORS + tuple(1.0 / f for f in FACTORS)
#: Reciprocal conversions: target = constant / value.
RECIPROCAL = (1239.84198, 1.0e7)
TYPED_TOOLS_EXCLUDED_FROM_SEEN = frozenset(
    {"consult_domain_skill", "open_guide"}
)


def significant_digits(text):
    digits = re.sub(r"[^0-9]", "", text.split("e")[0].split("E")[0])
    return len(digits.lstrip("0"))


def decimals(text):
    mantissa = text.split("e")[0].split("E")[0]
    return len(mantissa.split(".", 1)[1]) if "." in mantissa else 0


def numbers_in_text(text):
    """Every number in free text, as floats (any unit, any precision)."""

    out = []
    for match in NUMBER.finditer(text or ""):
        try:
            out.append(float(clean_number(match.group(0))))
        except ValueError:
            continue
    return out


def target_numbers(text):
    """Unit-bearing numbers with at least two significant digits."""

    out = []
    for match in NUMBER.finditer(text or ""):
        raw = clean_number(match.group(0))
        after = UNIT_AFTER.match(text, match.end())
        if not after:
            continue
        if significant_digits(raw) < 2:
            continue
        try:
            value = float(raw)
        except ValueError:
            continue
        start = max(0, match.start() - 70)
        out.append(
            {
                "value": value,
                "text": raw,
                "decimals": decimals(raw),
                "unit": after.group("unit"),
                "context": text[start : after.end() + 50].replace("\n", " "),
            }
        )
    return out


def numeric_values(value, strings=True):
    """Every number a JSON value carries: numeric leaves and, with
    strings, numbers written inside its string leaves."""

    out = []
    for _path, leaf in leaves(value):
        if isinstance(leaf, bool):
            continue
        if isinstance(leaf, (int, float)):
            out.append(float(leaf))
        elif strings and isinstance(leaf, str) and len(leaf) < 20000:
            out.extend(numbers_in_text(leaf))
    return out


def string_numbers(value):
    """Numbers written inside the string leaves of a JSON value only: the
    host's words (reports, meanings, diagnoses) rather than its values."""

    out = []
    for _path, leaf in leaves(value):
        if isinstance(leaf, str) and len(leaf) < 20000:
            out.extend(numbers_in_text(leaf))
    return out


#: A typed value is a record that names what it is and carries its value.
VALUE_KEYS = ("value", "source_value", "display_value")
ID_KEYS = frozenset(
    {
        "quantity_id",
        "claim_id",
        "output_id",
        "observable_id",
        "selector",
        "node_id",
    }
)


def quantity_values(value, scalars_only=False):
    """The values of typed quantity records inside a reply: every dict that
    carries an identifying key and a value key contributes the numbers under
    its value keys (scalars, and vector or matrix elements unless
    scalars_only), plus ``energy_hartree`` wherever it appears."""

    out = []

    def numbers(item):
        if isinstance(item, bool):
            return
        if isinstance(item, (int, float)):
            out.append(float(item))
        elif isinstance(item, list) and not scalars_only:
            for element in item:
                numbers(element)

    def visit(node):
        if isinstance(node, dict):
            if ID_KEYS.intersection(node) and any(
                k in node for k in VALUE_KEYS
            ):
                for key in VALUE_KEYS:
                    if key in node:
                        numbers(node[key])
            if "energy_hartree" in node:
                numbers(node["energy_hartree"])
            for item in node.values():
                if isinstance(item, (dict, list)):
                    visit(item)
        elif isinstance(node, list):
            for item in node:
                visit(item)

    visit(value)
    return out


def scalar_values(value):
    return quantity_values(value, scalars_only=True)


class SeenSet:
    """Sorted absolute magnitudes, for tolerance lookups."""

    def __init__(self, values):
        self.values = sorted({round(abs(v), 12) for v in values if v == v})

    def near(self, target, tolerance):
        """The seen magnitude within tolerance of target, or None."""

        index = bisect.bisect_left(self.values, target - tolerance)
        if (
            index < len(self.values)
            and self.values[index] <= target + tolerance
        ):
            return self.values[index]
        return None

    def __len__(self):
        return len(self.values)


def tolerance_of(number):
    return 0.5 * 10 ** (-number["decimals"]) * (1 + 1e-9) + 1e-12


def matches_direct(number, seen):
    return seen.near(abs(number["value"]), tolerance_of(number))


def matches_computed(number, singles, pairs):
    """One factor on a seen value; a sum or difference of two seen scalar
    values, optionally times one factor; or a reciprocal conversion.
    Returns (how, witness) or None; the witness names the operands."""

    x = abs(number["value"])
    tol = tolerance_of(number)
    for factor in FACTORS:
        hit = singles.near(x / factor, tol / factor)
        if hit is not None:
            return "factor", {"operand": hit, "factor": factor}
    for constant in RECIPROCAL:
        if x > 0:
            hit = singles.near(constant / x, tol * constant / (x * x) + 1e-12)
            if hit is not None:
                return "reciprocal", {"operand": hit, "constant": constant}
    values = pairs.values
    if not values:
        return None
    for factor in (1.0,) + FACTORS:
        y = x / factor
        t = tol / factor
        for v in values:
            # v - w = +-y  ->  w = v - y  or  w = v + y ; v + w = y -> w = y - v
            for w in (v - y, v + y, y - v):
                if w < 0 or abs(w - v) <= t:
                    continue
                hit = pairs.near(w, t)
                if hit is not None:
                    how = (
                        "difference"
                        if factor == 1.0
                        else "difference_converted"
                    )
                    return how, {"operands": (v, hit), "factor": factor}
    return None


def classify_number(number, direct, singles, pairs):
    """The first direct class whose seen set holds the number, else a
    one-step computation over typed scalars, else unmatched; with the
    witness that decided it."""

    for name, seen in direct:
        hit = matches_direct(number, seen)
        if hit is not None:
            return name, {"matched": hit}
    found = matches_computed(number, singles, pairs)
    if found:
        return "prose_computed:" + found[0], found[1]
    return "unmatched", {}


def precision_band(number):
    digits = significant_digits(number["text"])
    return "sig2-3" if digits <= 3 else "sig4+"


# --------------------------------------------------------------------------
# Sessions
# --------------------------------------------------------------------------

HEX64 = re.compile(r"\b[0-9a-f]{64}\b")
ARTIFACT_ID = re.compile(
    r"\b(?:orca|gaussian|pyscf|xtb|xyz|geometry|project|result|artifact)"
    r"-[a-z0-9-]*[0-9a-f]{12,16}\b"
)
CITATION_KEYS = re.compile(r"receipt|evidence_ref|sha256|digest")
ARTIFACT_KEYS = re.compile(
    r"artifact_id|artifact_ids|source_artifact|result_id"
)
DECISION_TEXT_KEYS = (
    "method_rationale",
    "assumptions",
    "alternatives",
    "uncertainties",
    "diagnostics",
)


def reply_view(content):
    data = parse_json(content)
    if not isinstance(data, dict):
        return {"status": None, "message": str(content)[:600]}
    status = data.get("status")
    message = ""
    if status == "rejected":
        error = data.get("error")
        if isinstance(error, dict):
            message = str(
                error.get("message") or error.get("diagnosis") or error
            )
        elif error:
            message = str(error)
        else:
            message = str(data.get("message") or "")
    return {
        "status": status,
        "error_class": data.get("error_class"),
        "message": message[:3000],
        "data": data,
    }


def stream_facts(directory):
    rows = load_jsonl(os.path.join(directory, "events.jsonl"))
    models = collections.Counter()
    provider_turns = 0
    terminal = {}
    for row in rows:
        kind = row.get("kind")
        payload = row.get("payload") or {}
        if kind == "provider_turn_observed":
            provider_turns += 1
            models[
                f"{payload.get('requested_model')}->{payload.get('observed_model')}"
            ] += 1
        elif kind == "runtime_terminated":
            terminal = {
                "terminal_state": payload.get("terminal_state"),
                "reason": str(payload.get("reason") or "")[:300],
            }
    return {
        "stream_events": len(rows),
        "provider_turns": provider_turns,
        "models": dict(models),
        **terminal,
    }


def workspace_of(path):
    marker = os.sep + ".chemsmart-agent" + os.sep
    return path.split(marker, 1)[0] if marker in path else None


class Goals:
    """Settlements and executed routes of every goal, by workspace."""

    def __init__(self):
        self.by_workspace = {}

    def of(self, workspace):
        if workspace is None:
            return ()
        if workspace in self.by_workspace:
            return self.by_workspace[workspace]
        goals = []
        root = os.path.join(workspace, ".chemsmart-agent", "goals")
        try:
            names = sorted(os.listdir(root))
        except OSError:
            names = []
        for name in names:
            ledger = load_jsonl(os.path.join(root, name, "ledger.jsonl"))
            if not ledger:
                continue
            settled = [r for r in ledger if r.get("kind") == "goal_settled"]
            payload = settled[-1].get("payload") if settled else {}
            goals.append(
                {
                    "goal_dir": name,
                    "settlement": (payload or {}).get("state"),
                    "reasons": [
                        str(r)[:300]
                        for r in (payload or {}).get("reasons") or ()
                    ],
                    "cycles": sum(
                        1 for r in ledger if r.get("kind") == "run_started"
                    ),
                    "route": executed_route(os.path.join(root, name, "runs")),
                }
            )
        self.by_workspace[workspace] = tuple(goals)
        return self.by_workspace[workspace]


def executed_route(runs_dir):
    nodes = collections.Counter()
    edges = collections.Counter()
    try:
        cycles = sorted(os.listdir(runs_dir))
    except OSError:
        cycles = []
    for cycle in cycles:
        for row in load_jsonl(os.path.join(runs_dir, cycle, "events.jsonl")):
            kind = row.get("kind")
            payload = row.get("payload") or {}
            if kind == "program_result_verified":
                record = payload.get("record") or {}
                observations = record.get("observations") or {}
                nodes[
                    f"{record.get('program')}:"
                    f"{observations.get('jobtype') or record.get('jobtype')}"
                ] += 1
            elif kind in (
                "optimized_geometry_handed_off",
                "workflow_data_edge_bound",
            ):
                edges[kind] += 1
            elif kind == "workflow_analysis_node_settled":
                edges["analysis:" + str(payload.get("analysis_kind"))] += 1
    return {"executed": dict(nodes), "edges": dict(edges)}


def planned_route(arguments):
    nodes = collections.Counter()
    edges = collections.Counter()
    analyses = collections.Counter()
    for node in arguments.get("calculation_nodes") or ():
        if not isinstance(node, dict):
            continue
        nodes[f"{node.get('program')}:{node.get('jobtype')}"] += 1
        for item in node.get("inputs") or ():
            if isinstance(item, dict):
                edges[
                    str(item.get("source_kind") or item.get("kind") or "?")
                ] += 1
    for node in arguments.get("analysis_nodes") or ():
        if not isinstance(node, dict):
            continue
        analyses[str(node.get("analysis_kind"))] += 1
    return {
        "nodes": dict(nodes),
        "input_edges": dict(edges),
        "analysis": dict(analyses),
    }


SEEN_KINDS = (
    "typed",
    "typed_other",
    "host_text",
    "typed_scalars",
    "context",
    "context_scalars",
    "own",
)


def _number_row(base, seen, **fields):
    """A number row carrying how much of each seen set preceded it."""

    return {
        **base,
        "detector": "number",
        **fields,
        "_seen": {kind: len(seen[kind]) for kind in SEEN_KINDS},
    }


def session_rows(path, data, declarations, goals):
    """Rows, a session record and the seen values of one distinct transcript."""

    messages = data.get("transcript") or data.get("messages") or []
    digest = data.get("transcript_sha256") or sha(
        json.dumps(messages, sort_keys=True)
    )
    directory = os.path.dirname(path)
    facts = stream_facts(directory)
    workspace = workspace_of(path)
    goal_list = goals.of(workspace)
    base = {"path": path, "session": digest}
    rows = []

    replies = {}
    for message in messages:
        if isinstance(message, dict) and message.get("role") == "tool":
            replies[message.get("tool_call_id")] = message.get("content")

    seen = {kind: [] for kind in SEEN_KINDS}
    minted_ids, context_ids = set(), set()
    counts = collections.Counter()
    key_use = {}
    last_assistant_text = ""
    last_assistant_index = None
    routes = []

    for index, message in enumerate(messages):
        if not isinstance(message, dict):
            continue
        role = message.get("role")
        content = message.get("content")
        if role == "user":
            text = content if isinstance(content, str) else json.dumps(content)
            parsed = parse_json(text)
            seen["context"].extend(
                numeric_values(parsed)
                if parsed is not None
                else numbers_in_text(text)
            )
            if parsed is not None:
                seen["context_scalars"].extend(scalar_values(parsed))
            context_ids.update(HEX64.findall(text))
            context_ids.update(ARTIFACT_ID.findall(text))
            counts["user_messages"] += 1
            continue
        if role == "tool":
            parsed = parse_json(content)
            name = message.get("name") or ""
            if name not in TYPED_TOOLS_EXCLUDED_FROM_SEEN:
                if parsed is not None:
                    seen["typed"].extend(quantity_values(parsed))
                    seen["typed_other"].extend(
                        numeric_values(parsed, strings=False)
                    )
                    seen["typed_scalars"].extend(scalar_values(parsed))
                    seen["host_text"].extend(string_numbers(parsed))
                else:
                    seen["host_text"].extend(numbers_in_text(str(content)))
            text = content if isinstance(content, str) else json.dumps(content)
            minted_ids.update(HEX64.findall(text))
            minted_ids.update(ARTIFACT_ID.findall(text))
            continue
        if role != "assistant":
            continue
        counts["assistant_messages"] += 1
        text = content if isinstance(content, str) else ""
        calls = message.get("tool_calls") or []
        if text.strip():
            last_assistant_text = text
            last_assistant_index = index
            for number in target_numbers(text):
                rows.append(
                    _number_row(
                        base,
                        seen,
                        surface="assistant_midsession",
                        delivered=False,
                        message_index=index,
                        number=number,
                    )
                )
        for call_index, call in enumerate(calls):
            counts["tool_calls"] += 1
            function = call.get("function") or {}
            name = function.get("name") or ""
            raw = function.get("arguments") or ""
            arguments = parse_json(raw) if isinstance(raw, str) else raw
            if not isinstance(arguments, dict):
                arguments = {}
            reply = reply_view(replies.get(call.get("id")))
            counts[f"reply:{reply['status']}"] += 1
            counts[f"tool:{name}"] += 1
            where = {
                **base,
                "message_index": index,
                "call_index": call_index,
                "tool": name,
                "reply_status": reply["status"],
            }

            # R10 Q28's hatch census, verbatim in substance.
            sections = arguments.get("sections")
            for match in UNKNOWN_KEY.finditer(
                str(replies.get(call.get("id")) or "")
            ):
                rows.append(
                    {
                        **where,
                        "detector": "unknown_key",
                        "program": arguments.get("program"),
                        "key": match.group(1),
                    }
                )
            if isinstance(sections, dict):
                counts["authoring_calls"] += 1
                program = str(arguments.get("program") or "?")
                for section, settings in sections.items():
                    if not isinstance(settings, dict):
                        continue
                    for key in settings:
                        key_use.setdefault(program, {}).setdefault(key, 0)
                        key_use[program][key] += 1
                    for key in ALL_KEYS:
                        if key in settings and nonempty(settings[key]):
                            rows.append(
                                {
                                    **where,
                                    "detector": "hatch",
                                    "program": arguments.get("program"),
                                    "section": section,
                                    "key": key,
                                    "value": settings[key],
                                    "reply_message": reply["message"][:600],
                                }
                            )
            elif isinstance(raw, str) and any(
                key in raw for key in HATCH_KEYS
            ):
                rows.append(
                    {
                        **where,
                        "detector": "hatch_other_tool",
                        "arguments_excerpt": raw[:600],
                    }
                )

            # X: native, path and shell strings anywhere in the arguments.
            for hit in string_exits(arguments):
                rows.append(
                    {
                        **where,
                        "detector": "string_exit",
                        **hit,
                        "reply_message": reply["message"][:400],
                    }
                )

            # V and R: every refusal, classified.
            if reply["status"] == "rejected":
                counts["refusals"] += 1
                verdict = classify_refusal(reply["message"])
                placed = [
                    (
                        selector,
                        place_selector(
                            selector,
                            verdict.get("program") or arguments.get("program"),
                            verdict.get("jobtype"),
                            declarations,
                        ),
                    )
                    for selector in verdict.get("requested") or ()
                ]
                rows.append(
                    {
                        **where,
                        "detector": "refusal",
                        **verdict,
                        "placed": placed,
                        "error_class": reply.get("error_class"),
                        "message": reply["message"][:1200],
                        "arguments_excerpt": (
                            raw if isinstance(raw, str) else json.dumps(raw)
                        )[:1200],
                    }
                )

            # V: the Agent's own statements of a vocabulary gap.
            if name == "plan_unsupported_external":
                rows.append(
                    {
                        **where,
                        "detector": "declared_gap",
                        "via": name,
                        "arguments_excerpt": str(raw)[:1500],
                    }
                )
            for node in arguments.get("analysis_nodes") or ():
                if isinstance(node, dict) and (
                    str(node.get("support_state") or "")
                    == "blocked_unsupported"
                    or str(node.get("blocked_reason") or "").strip()
                ):
                    rows.append(
                        {
                            **where,
                            "detector": "declared_gap",
                            "via": "blocked_unsupported",
                            "node_id": node.get("node_id"),
                            "analysis_kind": node.get("analysis_kind"),
                            "support_state": node.get("support_state"),
                            "blocked_reason": str(
                                node.get("blocked_reason") or ""
                            )[:1500],
                        }
                    )

            # R: citations of host-minted digests or ids no earlier reply
            # of this session carried.
            for leaf_path, value in leaves(arguments):
                if not isinstance(value, str):
                    continue
                key = ".".join(leaf_path)
                cited = []
                if CITATION_KEYS.search(key):
                    cited.extend(HEX64.findall(value))
                if ARTIFACT_KEYS.search(key):
                    cited.extend(ARTIFACT_ID.findall(value))
                for item in cited:
                    counts["citations"] += 1
                    if item in minted_ids:
                        counts["citations_in_session"] += 1
                        continue
                    origin = "context" if item in context_ids else "unknown"
                    counts[f"citations_{origin}"] += 1
                    rows.append(
                        {
                            **where,
                            "detector": "citation",
                            "origin": origin,
                            "leaf": key,
                            "cited": item,
                            "reply_message": reply["message"][:600],
                        }
                    )

            # N: decision prose (delivered findings, exploration fields).
            if name == "record_scientific_decision":
                delivered = reply["status"] == "ok"
                for finding in arguments.get("findings") or ():
                    if not isinstance(finding, dict):
                        continue
                    for number in target_numbers(
                        str(finding.get("statement") or "")
                    ):
                        rows.append(
                            _number_row(
                                where,
                                seen,
                                surface="finding",
                                delivered=delivered,
                                finding_id=finding.get("finding_id"),
                                number=number,
                            )
                        )
                for item in arguments.get("unreachable_observable_ids") or ():
                    if not isinstance(item, dict):
                        continue
                    for number in target_numbers(
                        str(item.get("statement") or "")
                    ):
                        rows.append(
                            _number_row(
                                where,
                                seen,
                                surface="unreachable",
                                delivered=delivered,
                                number=number,
                            )
                        )
                for field in DECISION_TEXT_KEYS:
                    value = arguments.get(field)
                    for text_item in (
                        value if isinstance(value, list) else [value]
                    ):
                        if isinstance(text_item, dict):
                            text_item = json.dumps(text_item)
                        for number in target_numbers(str(text_item or "")):
                            rows.append(
                                _number_row(
                                    where,
                                    seen,
                                    surface="decision_text",
                                    delivered=False,
                                    field=field,
                                    number=number,
                                )
                            )

            # S: the accepted plan's shape.
            if name == "plan_scientific_workflow" and reply["status"] == "ok":
                routes.append(planned_route(arguments))

            # The model's own numbers join the "own" set after its call.
            seen["own"].extend(numeric_values(arguments))

    # The last assistant message is delivered prose: its mid-session row
    # becomes the delivered row (same numbers, judged on everything seen).
    final_rows = []
    for row in rows:
        if (
            row.get("detector") == "number"
            and row.get("surface") == "assistant_midsession"
            and row.get("message_index") == last_assistant_index
        ):
            continue
        final_rows.append(row)
    if last_assistant_index is not None:
        for number in target_numbers(last_assistant_text):
            final_rows.append(
                _number_row(
                    base,
                    seen,
                    surface="final_message",
                    delivered=True,
                    message_index=last_assistant_index,
                    number=number,
                )
            )
    policy = False
    for message in messages:
        if isinstance(message, dict) and message.get("role") == "user":
            parsed = parse_json(message.get("content"))
            if isinstance(parsed, dict) and parsed.get(
                "analysis_completion_policy"
            ):
                policy = True
                break
    record = {
        **base,
        **facts,
        "run_dir": os.path.basename(directory),
        "workspace": workspace,
        "goals": [
            {k: g[k] for k in ("goal_dir", "settlement", "cycles")}
            for g in goal_list
        ],
        "analysis_completion_policy": policy,
        "counts": dict(counts),
        "key_use": key_use,
        "planned_routes": routes,
    }
    return final_rows, record, seen


def classify_numbers(rows, seen, control_seen=None):
    """Classify every number row against the seen values that preceded it
    (a number is judged on what the session had seen when it was written);
    with control_seen, the same matcher against another session's."""

    key = "control_class" if control_seen is not None else "class"
    for row in rows:
        if row.get("detector") != "number":
            continue
        cut = row.get("_seen") or {}
        if control_seen is None:
            source = {
                kind: seen[kind][: cut.get(kind, 0)] for kind in SEEN_KINDS
            }
        else:
            source = control_seen
        direct = (
            ("bound", SeenSet(source["typed"])),
            ("bound_other", SeenSet(source["typed_other"])),
            ("host_text", SeenSet(source["host_text"])),
            ("context", SeenSet(source["context"])),
            ("own_argument", SeenSet(source["own"])),
        )
        scalars = source["typed_scalars"] + source["context_scalars"]
        row[key], row[key + "_witness"] = classify_number(
            row["number"], direct, SeenSet(scalars), SeenSet(scalars)
        )
        row["precision"] = precision_band(row["number"])


# --------------------------------------------------------------------------
# Project files (Q28, verbatim in substance)
# --------------------------------------------------------------------------


def project_file_rows(path):
    try:
        text = open(path, encoding="utf-8").read()
    except OSError:
        return []
    rows = []
    data = None
    if load_yaml is not None:
        try:
            data = load_yaml(text)
        except Exception:  # noqa: BLE001
            data = None
    if isinstance(data, dict):
        for section, settings in data.items():
            if not isinstance(settings, dict):
                continue
            for key in ALL_KEYS:
                if key in settings and nonempty(settings[key]):
                    rows.append(
                        {
                            "detector": "project_file",
                            "path": path,
                            "file_sha256": sha(text),
                            "section": section,
                            "key": key,
                            "value": settings[key],
                        }
                    )
    return rows


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------


def dump_declarations(path):
    """Write the importing tree's reader declarations and operation names.

    The one place this file imports chemsmart, and only when asked: the
    census itself stays importable on a machine with no chemsmart, and the
    table it places selectors against is computed from the tree it names,
    never written by hand.
    """

    import chemsmart
    from chemsmart.agent.tool_specs import OPERATION_DESCRIPTIONS
    from chemsmart.analysis.result_readers import RESULT_READERS

    table = {
        "chemsmart_file": chemsmart.__file__,
        "programs": {
            program: {
                "selectors": sorted(reader.selectors),
                "jobtypes": {
                    jobtype: sorted(selectors)
                    for jobtype, selectors in reader.jobtype_selectors
                },
            }
            for program, reader in sorted(RESULT_READERS.items())
        },
        "operations": sorted(OPERATION_DESCRIPTIONS),
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(table, handle, indent=1, sort_keys=True)
    print(f"declarations from {chemsmart.__file__} -> {path}")


def main(argv):
    if len(argv) == 3 and argv[1] == "--dump-declarations":
        dump_declarations(argv[2])
        return
    out_dir = argv[1]
    roots, excludes, declarations = [], [], None
    items = iter(argv[2:])
    for item in items:
        if item == "--exclude":
            excludes.append(next(items))
        elif item == "--declarations":
            declarations = json.load(open(next(items), encoding="utf-8"))
        else:
            roots.append(item)
    os.makedirs(out_dir, exist_ok=True)
    goals = Goals()
    seen_digests = set()
    seen_project_files = set()
    sessions, all_rows, seen_sets = [], [], []
    copies = 0
    for root in roots:
        mine = [
            p.split(":", 1)[1] for p in excludes if p.split(":", 1)[0] == root
        ]
        for dirpath, files in walk(root, mine):
            for name in sorted(files):
                path = os.path.join(dirpath, name)
                if name.startswith("public-transcript") and name.endswith(
                    ".json"
                ):
                    try:
                        data = json.load(open(path, encoding="utf-8"))
                    except Exception:  # noqa: BLE001
                        continue
                    if not isinstance(data, dict):
                        data = {"transcript": data}
                    digest = data.get("transcript_sha256")
                    if digest and digest in seen_digests:
                        copies += 1
                        continue
                    rows, record, seen = session_rows(
                        path, data, declarations, goals
                    )
                    seen_digests.add(record["session"])
                    record["root"] = root
                    sessions.append(record)
                    seen_sets.append((rows, seen))
                elif (
                    name.endswith((".yaml", ".yml"))
                    and os.sep + "projects" + os.sep in path
                    and ".chemsmart-agent" in path
                ):
                    for row in project_file_rows(path):
                        if row["file_sha256"] in seen_project_files:
                            continue
                        seen_project_files.add(row["file_sha256"])
                        row["root"] = root
                        all_rows.append(row)
    # Numbers: real classification, then the control against a fixed
    # derangement (session i judged on session i+1's seen values).
    order = sorted(range(len(sessions)), key=lambda i: sessions[i]["session"])
    for position, i in enumerate(order):
        rows, seen = seen_sets[i]
        partner = (
            seen_sets[order[(position + 1) % len(order)]][1]
            if len(order) > 1
            else None
        )
        control_rows = [dict(r) for r in rows if r.get("detector") == "number"]
        classify_numbers(control_rows, seen, control_seen=partner or seen)
        classify_numbers(rows, seen)
        controls = iter(control_rows)
        for row in rows:
            if row.get("detector") == "number":
                row["control_class"] = next(controls).get("control_class")
        for row in rows:
            row.pop("_seen", None)
            row["root"] = sessions[i]["root"]
            row["models"] = sessions[i].get("models")
        all_rows.extend(rows)
    goal_records = {}
    for workspace, goal_list in goals.by_workspace.items():
        for goal in goal_list:
            goal_records[f"{workspace}::{goal['goal_dir']}"] = goal
    with open(
        os.path.join(out_dir, "rows.jsonl"), "w", encoding="utf-8"
    ) as handle:
        for row in all_rows:
            handle.write(json.dumps(row, sort_keys=True, default=str) + "\n")
    with open(
        os.path.join(out_dir, "sessions.jsonl"), "w", encoding="utf-8"
    ) as handle:
        for record in sessions:
            handle.write(
                json.dumps(record, sort_keys=True, default=str) + "\n"
            )
    with open(
        os.path.join(out_dir, "goals.json"), "w", encoding="utf-8"
    ) as handle:
        json.dump(goal_records, handle, indent=1, sort_keys=True, default=str)
    summary = summarise(
        sessions, all_rows, goal_records, roots, excludes, copies
    )
    with open(
        os.path.join(out_dir, "summary.json"), "w", encoding="utf-8"
    ) as handle:
        json.dump(summary, handle, indent=1, sort_keys=True, default=str)
    print(json.dumps(summary["headline"], indent=1, sort_keys=True))


def infrastructure(record):
    return record.get(
        "provider_turns", 0
    ) == 0 or "turn_deadline_exceeded" in str(record.get("reason") or "")


def summarise(sessions, rows, goal_records, roots, excludes, copies):
    behavioural = [s for s in sessions if not infrastructure(s)]
    live = {s["session"] for s in behavioural}
    by = collections.defaultdict(collections.Counter)
    for row in rows:
        if row.get("session") and row["session"] not in live:
            by["excluded_infrastructure_rows"][row.get("detector")] += 1
            continue
        detector = row.get("detector")
        if detector == "refusal":
            by["refusal_family"][row.get("family")] += 1
            by["refusal_class"][row.get("class")] += 1
            for selector, placement in row.get("placed") or ():
                by["selector_placement"][str(placement)] += 1
        elif detector == "number":
            surface = f"{row.get('surface')}:{row.get('precision')}"
            by[f"number:{surface}"][row.get("class")] += 1
            by[f"number_control:{surface}"][row.get("control_class")] += 1
        elif detector == "string_exit":
            by["string_exit"][
                f"{row.get('family')}:{row.get('pattern')}:"
                f"{'prose' if row.get('prose_field') else 'operative'}"
            ] += 1
        elif detector == "citation":
            by["citation_origin"][
                f"{row.get('origin')}:{row.get('reply_status')}"
            ] += 1
        elif detector == "declared_gap":
            by["declared_gap"][row.get("via")] += 1
        elif detector == "hatch_other_tool":
            by[detector][row.get("tool")] += 1
        elif detector in ("hatch", "unknown_key", "project_file"):
            by[detector][row.get("key")] += 1
    counts = collections.Counter()
    models = collections.Counter()
    for record in behavioural:
        counts.update(record.get("counts") or {})
        for model in record.get("models") or {}:
            models[model] += 1
    settlements = collections.Counter(
        str(goal.get("settlement")) for goal in goal_records.values()
    )
    hatch_sessions = {
        r["session"]
        for r in rows
        if r.get("detector") == "hatch" and r.get("session") in live
    }
    headline = {
        "roots": roots,
        "excludes": excludes,
        "transcript_copies_skipped": copies,
        "sessions_distinct": len(sessions),
        "sessions_infrastructure": len(sessions) - len(behavioural),
        "sessions_behavioural": len(behavioural),
        "models_by_session": dict(models),
        "tool_calls": counts.get("tool_calls", 0),
        "authoring_calls": counts.get("authoring_calls", 0),
        "refusals": counts.get("refusals", 0),
        "citations": counts.get("citations", 0),
        "citations_in_session": counts.get("citations_in_session", 0),
        "sessions_with_hatch": len(hatch_sessions),
        "goals": len(goal_records),
        "goal_settlements": dict(settlements),
    }
    return {
        "headline": headline,
        "tallies": {k: dict(v) for k, v in sorted(by.items())},
    }


if __name__ == "__main__":
    main(sys.argv)
