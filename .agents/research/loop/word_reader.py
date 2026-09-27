"""Is each word the host signed true of the records it cites?

    python .agents/research/loop/word_reader.py OUT_DIR label=ROOT [...]
    python .agents/research/loop/word_reader.py OUT_DIR label=ROOT [...] \
        --words RESULTS_JSONL      # read the settlement words resign.py
                                   # signed instead of the archived ones

An independent reader: it imports nothing from ``chemsmart``, so it cannot
restate the signer. For every goal ledger under the roots (the same
discovery and deduplication as ``resign.py``) it reads the ledger and every
stream the goal wrote -- the planning sessions the ledger names (or, for a
ledger older than ``session_stream_recorded``, the live streams inside the
goal's lifetime, marked inferred), the analysis-evidence sessions, and the
cycle run streams -- and checks each class of signed word
(``signed_words.py``) against those records:

- W1 settlement: an achieved word over a delivery whose newest completion
  did not pass, over a declared id whose latest typed word is a verified
  refusal, over a declared id nobody delivered; a plain ``achieved`` over a
  standing host anomaly, a falsified expectation or an answered failed
  criterion (the first word must not hide what the run found); an
  ``unreachable_from_evidence`` whose refused ids have no verified refusal,
  were delivered after it, or leave another declared id undelivered; an
  ``exhausted`` with every budget left.
- W3 ``qualified`` rows under a word that is not achieved, or under a
  flagged achieved word.
- W4 the executor's analysis word against its run stream.
- W5 a ``passed`` completion that lists an unanswered failed criterion.
- W6 an expectation verdict against the arithmetic of its own row.
- W7 a verified refusal whose blocked-node basis is in no plan the goal
  recorded, or whose refused selector a later extraction read.
- W8 a finding's standing against the declarations it had, and each
  relation's truth against the values it names.
- W9 a category's answer against the extraction receipt it names.
- W10 a characterisation's order against its observed mode count.
- W11 a free energy derived from a result its verification did not pass.

Every check reports ok / flag / insufficient; a flag is a candidate to be
read against its records, never a verdict. It also lists the dissent
markers each goal holds (falsified expectations and diagnostics, answered
criteria, verified and unverified refusals, unrequested findings,
superseded findings and observables, doubts, rejected routes) and whether
the settlement's reasons name them. Written for R11 episode `truth`;
Q24's ``signed_word_violations`` (tests/agent) is the ancestor of the W1,
W3 and W4 checks, re-derived here so the census reads timestamps rather
than cycle order and can see ledgers that name no session stream.
"""

from __future__ import annotations

import collections
import json
import os
import re
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

ACHIEVED = {"achieved", "achieved_with_observations"}
EXCLUDED_PARTS = ("sealed", "/private", "/claude/", "/grade-", "/grading/")
HESS_GRADIENT = 4.5e-4


def rows(path: Path) -> list[dict]:
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    out = []
    for line in text.splitlines():
        if line.strip():
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return out


def when(text: str) -> str:
    """ISO timestamps compare as strings once normalised to UTC seconds."""

    try:
        stamp = datetime.fromisoformat(str(text).replace("Z", "+00:00"))
    except ValueError:
        return str(text)
    if stamp.tzinfo is not None:
        stamp = stamp.astimezone(timezone.utc).replace(tzinfo=None)
    return stamp.strftime("%Y-%m-%dT%H:%M:%S.%f")


def discover(specs: list[str]) -> list[dict]:
    found: list[dict] = []
    for spec in specs:
        label, root = spec.split("=", 1)
        for folder, dirs, files in os.walk(root):
            dirs[:] = [
                d
                for d in dirs
                if not any(p in d for p in ("sealed", "private"))
                and d not in {"scratch", "claude"}
            ]
            if "ledger.jsonl" not in files:
                continue
            ledger = Path(folder) / "ledger.jsonl"
            text = str(ledger)
            if "/.chemsmart-agent/goals/" not in text:
                continue
            if any(part in text for part in EXCLUDED_PARTS):
                continue
            found.append(
                {
                    "label": label,
                    "agent": str(ledger.parents[2]),
                    "goal_id": ledger.parent.name,
                }
            )
    unique: dict[tuple[str, str], dict] = {}
    for item in found:
        goal_file = Path(item["agent"]) / "goals" / item["goal_id"] / "goal.json"
        try:
            digest = json.loads(goal_file.read_text()).get("goal_sha256", "")
        except (OSError, json.JSONDecodeError):
            digest = ""
        item["goal_sha256"] = digest
        key = (item["goal_id"], digest or item["agent"])
        if key not in unique:
            unique[key] = item
    return sorted(unique.values(), key=lambda i: (i["label"], i["agent"]))


def walk(obj, depth=0):
    """Every mapping inside a payload, depth-first."""

    if depth > 12:
        return
    if isinstance(obj, dict):
        yield obj
        for value in obj.values():
            yield from walk(value, depth + 1)
    elif isinstance(obj, list):
        for value in obj:
            yield from walk(value, depth + 1)


class Goal:
    """One goal's records, read once."""

    def __init__(self, agent: Path, goal_id: str, word: dict | None = None):
        self.agent = agent
        self.goal_id = goal_id
        self.ledger = rows(agent / "goals" / goal_id / "ledger.jsonl")
        try:
            self.record = json.loads(
                (agent / "goals" / goal_id / "goal.json").read_text()
            )
        except (OSError, json.JSONDecodeError):
            self.record = {}
        settled = [r for r in self.ledger if r.get("kind") == "goal_settled"]
        self.settled = settled[-1] if settled else None
        payload = (self.settled or {}).get("payload") or {}
        self.word = str(payload.get("state") or "")
        self.reasons = [str(item) for item in payload.get("reasons") or ()]
        self.qualified = [
            r.get("payload") or {}
            for r in self.ledger
            if r.get("kind") == "qualified"
        ]
        if word is not None:
            # The word another tree signed over these records.
            self.word = str(word.get("state") or "")
            self.reasons = [str(item) for item in word.get("reasons") or ()]
            self.qualified = [{"id": q} for q in word.get("qualified") or ()]
        self.streams = self._streams()
        self.events: list[tuple[str, str, str, dict]] = []
        for role, path, inferred in self.streams:
            for event in rows(path):
                self.events.append(
                    (
                        when(event.get("timestamp") or ""),
                        role,
                        str(path),
                        event,
                    )
                )
        self.events.sort(key=lambda item: item[0])
        self.declared = self._declared()

    # -- streams --------------------------------------------------------

    def _streams(self) -> list[tuple[str, Path, bool]]:
        out: list[tuple[str, Path, bool]] = []
        seen: set[str] = set()

        def add(role: str, path: Path, inferred: bool = False) -> None:
            key = str(path.resolve()) if path.exists() else str(path)
            if key in seen or not path.exists():
                return
            seen.add(key)
            out.append((role, path, inferred))

        named = False
        for entry in self.ledger:
            payload = entry.get("payload") or {}
            kind = entry.get("kind")
            if kind == "session_stream_recorded" and payload.get("run_id"):
                named = True
                add(
                    "session",
                    self.agent / "runs" / str(payload["run_id"]) / "events.jsonl",
                )
            elif kind == "analysis_evidence_recorded" and payload.get("evidence"):
                add(
                    "session",
                    self.agent
                    / Path(*str(payload["evidence"]).split("/"))
                    / "events.jsonl",
                )
            elif kind == "run_recorded" and payload.get("run"):
                add(
                    "run",
                    self.agent
                    / Path(*str(payload["run"]).split("/"))
                    / "events.jsonl",
                )
        for stream in sorted(
            (self.agent / "goals" / self.goal_id / "runs").glob(
                "cycle-*/events.jsonl"
            )
        ):
            add("run", stream)
        if not named:
            created = [r for r in self.ledger if r.get("kind") == "goal_created"]
            start = end = ""
            if created:
                try:
                    opened = datetime.fromisoformat(str(created[0].get("at")))
                    start = (opened - timedelta(hours=2)).strftime(
                        "%Y%m%dT%H%M%S"
                    )
                except ValueError:
                    start = ""
            if self.settled:
                end = (
                    str(self.settled.get("at") or "")
                    .replace("-", "")
                    .replace(":", "")
                    .split(".")[0]
                    .split("+")[0][:15]
                )
            for stream in sorted((self.agent / "runs").glob("live-*/events.jsonl")):
                name = stream.parent.name
                began = name.split("-")[1][:15] if "-" in name else ""
                if start and began and began < start:
                    continue
                if end and began and began > end:
                    continue
                add("session", stream, inferred=True)
        return out

    @property
    def inferred(self) -> bool:
        return any(inferred for _r, _p, inferred in self.streams)

    def of(self, *kinds: str):
        for stamp, role, path, event in self.events:
            if event.get("kind") in kinds:
                yield stamp, role, path, event.get("payload") or {}

    # -- declarations and deliveries ------------------------------------

    def _declared(self) -> dict[str, dict]:
        declared: dict[str, dict] = {}
        for entry in self.ledger:
            if entry.get("kind") != "observables_declared":
                continue
            for item in (entry.get("payload") or {}).get("observables") or ():
                oid = str(item.get("observable_id") or "")
                if oid and oid not in declared:
                    declared[oid] = dict(item)
                    declared[oid]["_at"] = when(entry.get("at") or "")
        for stamp, _role, _path, payload in self.of("requested_observable_declared"):
            record = payload.get("record") or payload
            oid = str(record.get("observable_id") or payload.get("observable_id") or "")
            if oid and oid not in declared:
                declared[oid] = dict(record)
                declared[oid]["_at"] = stamp
        return declared

    def retired(self) -> set[str]:
        out = set()
        for item in self.declared.values():
            if item.get("supersedes_observable_id"):
                out.add(str(item["supersedes_observable_id"]))
        for entry in self.ledger:
            if entry.get("kind") == "observables_declared":
                for item in (entry.get("payload") or {}).get("observables") or ():
                    if item.get("supersedes_observable_id"):
                        out.add(str(item["supersedes_observable_id"]))
        return out

    def required(self) -> set[str]:
        """Declared, requested (not diagnostic), not retired."""

        return {
            oid
            for oid, item in self.declared.items()
            if str(item.get("role") or "requested") != "diagnostic"
        } - self.retired()

    def claims(self):
        """(stamp, claim) for every claim any stream rendered."""

        for stamp, _role, _path, payload in self.of("analysis_claims_recorded"):
            for claim in (payload.get("record") or {}).get("claims") or ():
                yield stamp, claim

    def claimed_ids(self) -> dict[str, str]:
        latest: dict[str, str] = {}
        for stamp, claim in self.claims():
            for key in ("claim_id", "quantity_id"):
                name = str(claim.get(key) or "")
                if name:
                    latest[name] = max(latest.get(name, ""), stamp)
        # A finding that answers a declared category delivers it.
        for stamp, _r, _p, payload in self.of("scientific_decision_recorded"):
            for finding in payload.get("findings") or ():
                name = str(finding.get("answers_observable_id") or "")
                if name and finding.get("answer"):
                    latest[name] = max(latest.get(name, ""), stamp)
        record = self.agent / "workspace-record.jsonl"
        for row in rows(record):
            if str(row.get("goal_id") or "") not in {"", self.goal_id}:
                continue
            if str(row.get("kind") or "") != "claim":
                continue
            for key in ("claim_id", "quantity_id", "observable_id"):
                name = str(row.get(key) or "")
                if name:
                    latest.setdefault(name, "")
        return latest

    def refusals(self):
        """(stamp, path, item) for every typed refusal recorded."""

        for stamp, _role, path, payload in self.of("scientific_decision_recorded"):
            for item in payload.get("unreachable_observables") or ():
                yield stamp, path, item

    def completions(self):
        return list(self.of("analysis_completion_evaluated"))


# -- the checks -------------------------------------------------------------


class Report:
    def __init__(self) -> None:
        self.counts: dict[tuple[str, str], collections.Counter] = (
            collections.defaultdict(collections.Counter)
        )
        self.flags: list[dict] = []

    def add(self, goal: Goal, cls: str, check: str, verdict: str, detail=""):
        self.counts[(cls, check)][verdict] += 1
        if verdict == "flag":
            self.flags.append(
                {
                    "goal_id": goal.goal_id,
                    "agent": str(goal.agent),
                    "word": goal.word,
                    "class": cls,
                    "check": check,
                    "detail": detail,
                }
            )


def check_settlement(goal: Goal, report: Report) -> list[str]:
    """W1; returns the flags raised on the word itself."""

    raised: list[str] = []
    if not goal.word:
        return raised
    required = goal.required()
    claimed = goal.claimed_ids()
    completions = goal.completions()
    latest = str(completions[-1][3].get("status") or "") if completions else ""

    def flag(check: str, detail: str) -> None:
        raised.append(check)
        report.add(goal, "W1", check, "flag", detail)

    if goal.word in ACHIEVED:
        if latest != "passed" and (latest or required):
            flag(
                "achieved_over_uncertified",
                f"newest completion {latest or 'absent'}"
                + (" (session streams inferred)" if goal.inferred else ""),
            )
        else:
            report.add(goal, "W1", "achieved_over_uncertified", "ok")
        # Latest typed word per id: a claim, or a verified refusal of
        # presence (a refused precision stands on its number).
        refused_last = []
        for stamp, _path, item in goal.refusals():
            oid = str(item.get("observable_id") or "")
            presence = item.get("blocked_node_id") or item.get("selector")
            tolerance = (goal.declared.get(oid) or {}).get("required_tolerance")
            if (
                item.get("verified")
                and presence
                and tolerance is None
                and oid in required
                and stamp > claimed.get(oid, "")
            ):
                refused_last.append(oid)
        if refused_last:
            flag(
                "achieved_over_later_refusal",
                ", ".join(sorted(set(refused_last))),
            )
        else:
            report.add(goal, "W1", "achieved_over_later_refusal", "ok")
        verified = {
            str(item.get("observable_id") or "")
            for _s, _p, item in goal.refusals()
            if item.get("verified")
        }
        undelivered = sorted(
            oid for oid in required if oid not in claimed and oid not in verified
        )
        if undelivered:
            flag("achieved_over_undelivered", ", ".join(undelivered))
        else:
            report.add(goal, "W1", "achieved_over_undelivered", "ok")
    if goal.word == "achieved":
        hidden = []
        standing: dict[tuple[str, str], str] = {}
        for entry in goal.ledger:
            if entry.get("kind") != "anomalies_observed":
                continue
            for item in (entry.get("payload") or {}).get("anomalies") or ():
                key = (str(item.get("signal_id")), str(item.get("node_id")))
                standing[key] = str(item.get("status") or "")
        hidden += [
            f"anomaly {signal} on {node} ({status})"
            for (signal, node), status in sorted(standing.items())
            if status != "refuted"
        ]
        last_row: dict[str, dict] = {}
        for _stamp, _r, _p, payload in completions:
            for row in payload.get("declared_observable_predictions") or ():
                oid = str(row.get("observable_id") or "")
                if oid and row.get("delivered_claim_id"):
                    last_row[oid] = row
        hidden += [
            f"falsified expectation {oid}"
            for oid, row in sorted(last_row.items())
            if row.get("agreement") == "diverged" and oid not in goal.retired()
        ]
        if completions:
            hidden += [
                f"answered criterion {item}"
                for item in completions[-1][3].get("anomaly_output_ids") or ()
                if str(item).startswith("failed_criterion:")
                and ":answered:" in str(item)
            ]
        if hidden:
            flag("achieved_hides_what_the_run_found", "; ".join(hidden))
        else:
            report.add(goal, "W1", "achieved_hides_what_the_run_found", "ok")
    if goal.word in ACHIEVED | {"unreachable_from_evidence"}:
        # Every pre-registered expectation the physics left, for an id the
        # goal still delivers, is named by the word's reasons: the latest
        # row per id across every completion of the goal, any cycle.
        text = " ".join(goal.reasons)
        last: dict[str, dict] = {}
        for _stamp, _r, _p, payload in completions:
            for row in payload.get("declared_observable_predictions") or ():
                oid = str(row.get("observable_id") or "")
                if oid and row.get("delivered_claim_id"):
                    last[oid] = row
        # An id whose latest typed word is a verified refusal delivers no
        # number, so an earlier claim's verdict no longer describes it.
        # A refused precision leaves its number standing (the host's own
        # rule, driver._claims_a_later_refusal_supersedes), so only an id
        # declared without a tolerance is superseded by its refusal.
        refused_after: set[str] = set()
        for stamp, _path, item in goal.refusals():
            oid = str(item.get("observable_id") or "")
            tolerance = (goal.declared.get(oid) or {}).get("required_tolerance")
            if (
                item.get("verified")
                and tolerance is None
                and stamp > claimed.get(oid, "")
            ):
                refused_after.add(oid)
        unnamed = sorted(
            oid
            for oid, row in last.items()
            if row.get("agreement") == "diverged"
            and oid not in goal.retired()
            and oid not in refused_after
            and f"falsified_expectation:{oid}" not in text
        )
        if unnamed:
            flag("reasons_name_each_falsified_expectation", ", ".join(unnamed))
        else:
            report.add(goal, "W1", "reasons_name_each_falsified_expectation", "ok")
    if goal.word == "unreachable_from_evidence":
        text = " ".join(goal.reasons)
        named = sorted(
            oid for oid in goal.declared if f"{oid} -- " in text
        )
        if not named:
            report.add(goal, "W1", "unreachable_names_its_refusals", "insufficient")
        else:
            verified_at: dict[str, str] = {}
            for stamp, _path, item in goal.refusals():
                if item.get("verified"):
                    oid = str(item.get("observable_id") or "")
                    verified_at[oid] = max(verified_at.get(oid, ""), stamp)
            missing = [oid for oid in named if oid not in verified_at]
            later = [
                oid
                for oid in named
                if oid in verified_at
                and (goal.declared.get(oid) or {}).get("required_tolerance") is None
                and claimed.get(oid, "") > verified_at[oid]
            ]
            others = sorted(
                oid
                for oid in required
                if oid not in named and oid not in claimed and oid not in verified_at
            )
            if missing:
                flag("unreachable_without_verified_refusal", ", ".join(missing))
            if later:
                flag("unreachable_over_later_delivery", ", ".join(later))
            if others:
                flag("unreachable_leaves_other_ids_undelivered", ", ".join(others))
            if not (missing or later or others):
                report.add(goal, "W1", "unreachable_names_its_refusals", "ok")
    if goal.word == "exhausted":
        envelope = goal.record.get("envelope") or {}
        calls = int(envelope.get("max_engine_calls") or 0)
        wall = float(envelope.get("episode_wall_time_seconds") or 0.0)
        revisions = int(goal.record.get("max_revisions") or 0)
        for entry in goal.ledger:
            payload = entry.get("payload") or {}
            if entry.get("kind") == "run_recorded":
                calls -= int(payload.get("engine_calls_consumed") or 0)
                wall -= float(payload.get("engine_wall_seconds") or 0.0)
            elif entry.get("kind") in {"revision_admitted", "rewake_opened"}:
                revisions -= 1
        if not goal.record:
            report.add(goal, "W1", "exhausted_has_a_spent_budget", "insufficient")
        elif calls > 0 and wall > 0 and revisions > 0:
            flag(
                "exhausted_with_budget_left",
                f"calls {calls}, wall {wall:.0f} s, revisions {revisions}",
            )
        else:
            report.add(goal, "W1", "exhausted_has_a_spent_budget", "ok")
    return raised


def check_qualified(goal: Goal, report: Report, raised: list[str]) -> None:
    if not goal.qualified:
        return
    if goal.word not in ACHIEVED:
        report.add(goal, "W3", "qualified_under_its_word", "flag",
                   f"{len(goal.qualified)} rows under {goal.word or 'unsettled'}")
    elif raised:
        report.add(goal, "W3", "qualified_under_its_word", "flag",
                   f"{len(goal.qualified)} rows under a word flagged {raised}")
    else:
        report.add(goal, "W3", "qualified_under_its_word", "ok")


def check_executor(goal: Goal, report: Report) -> None:
    for entry in goal.ledger:
        if entry.get("kind") != "run_recorded":
            continue
        payload = entry.get("payload") or {}
        run = str(payload.get("run") or "")
        if not run:
            continue
        stream = goal.agent / Path(*run.split("/")) / "events.jsonl"
        word = payload.get("analysis_status")
        result_file = stream.parent / "execution-result.json"
        if word is None and result_file.is_file():
            try:
                word = json.loads(result_file.read_text()).get("analysis_status")
            except (OSError, json.JSONDecodeError):
                word = None
        if word is None:
            report.add(goal, "W4", "executor_word_is_its_walk", "insufficient")
            continue
        events = rows(stream)
        receipts = [e for e in events if e.get("kind") == "analysis_completion_evaluated"]
        refused = any(e.get("kind") == "workflow_analysis_completion_refused" for e in events)
        ran = [
            (e.get("payload") or {}).get("node_id")
            for e in events
            if e.get("kind") == "workflow_analysis_node_settled"
            and (e.get("payload") or {}).get("state") != "blocked_unsupported"
        ]
        problem = ""
        if word == "completed" and not receipts:
            problem = "completed over a stream with no completion receipt"
        elif word == "partial" and not (receipts or refused):
            problem = "partial over neither a completion receipt nor a refusal"
        elif word == "" and (receipts or ran):
            problem = f"'' over {len(receipts)} receipt(s) and nodes {ran[:4]}"
        report.add(
            goal, "W4", "executor_word_is_its_walk",
            "flag" if problem else "ok",
            f"cycle {payload.get('cycle')}: {problem}" if problem else "",
        )


def check_completions(goal: Goal, report: Report) -> None:
    for stamp, _role, path, payload in goal.completions():
        status = str(payload.get("status") or "")
        unanswered = [
            item
            for item in payload.get("anomaly_output_ids") or ()
            if str(item).startswith("failed_criterion:") and ":unanswered:" in str(item)
        ]
        if status == "passed" and unanswered:
            report.add(goal, "W5", "passed_lists_no_unanswered_criterion", "flag",
                       f"{Path(path).parent.name}: {unanswered}")
        else:
            report.add(goal, "W5", "passed_lists_no_unanswered_criterion", "ok")
        for row in payload.get("declared_observable_predictions") or ():
            check_prediction(goal, report, row)


def _number(value):
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def check_prediction(goal: Goal, report: Report, row: dict) -> None:
    """W6: the verdict from the row's own numbers."""

    recorded = str(row.get("agreement") or "")
    value = _number(row.get("delivered_value"))
    if value is None or not row.get("delivered_claim_id"):
        report.add(goal, "W6", "verdict_is_its_arithmetic",
                   "ok" if recorded == "not_comparable" else "insufficient")
        return
    declared_unit = str(
        (goal.declared.get(str(row.get("observable_id") or "")) or {}).get("unit")
        or ""
    )
    if "delivered_value_in_declared_unit" in row:
        comparable = _number(row.get("delivered_value_in_declared_unit"))
    elif row.get("band_untestable"):
        comparable = None
    elif declared_unit and str(row.get("delivered_unit") or "") != declared_unit:
        # A row written before the host converted units: the band was in
        # another unit and is untestable here; a sign is unit-free.
        comparable = None
    else:
        comparable = value
    verdicts = []
    sign = str(row.get("expected_sign") or "")
    if sign and not row.get("sign_implied_by_band"):
        verdicts.append(value != 0.0 and (("positive" if value > 0 else "negative") == sign))
    low, high = row.get("expected_low"), row.get("expected_high")
    if _number(low) is not None and _number(high) is not None:
        verdicts.append(None if comparable is None else (float(low) <= comparable <= float(high)))
    resolution = _number(row.get("method_resolution"))
    if resolution is not None and comparable is not None and abs(comparable) <= resolution:
        expected = {"indeterminate"}
    elif any(v is False for v in verdicts):
        expected = {"diverged"}
    elif verdicts and None not in verdicts:
        expected = {"agreed", "agreed_as_approximation"}
    else:
        expected = {"not_comparable"}
    if recorded in expected:
        report.add(goal, "W6", "verdict_is_its_arithmetic", "ok")
    else:
        report.add(goal, "W6", "verdict_is_its_arithmetic", "flag",
                   f"{row.get('observable_id')}: recorded {recorded}, "
                   f"arithmetic {sorted(expected)} (value {value}, "
                   f"band {low}..{high}, sign {sign or '-'})")


def _plans_hold_blocked(goal: Goal, node_id: str, output_id: str, path: str | None):
    """Where a blocked_unsupported node with that output was recorded.

    Returns (places, described): the streams whose records hold the node
    blocked with that output, and whether any record describes the node's
    outputs at all. Tool arguments are recorded only as a digest, so a node
    built through the plan-draft constructors can be named by id in the
    records with its outputs nowhere: that is an insufficient record, not a
    false basis.
    """

    places = []
    described = False
    for _stamp, _role, stream, event in goal.events:
        if path is not None and stream != path:
            continue
        payload = event.get("payload") or {}
        for mapping in walk(payload):
            if mapping.get("node_id") != node_id:
                continue
            if "outputs" in mapping or "support_state" in mapping:
                described = True
            if mapping.get("support_state") != "blocked_unsupported" and not (
                event.get("kind") == "workflow_analysis_node_settled"
                and mapping.get("state") == "blocked_unsupported"
            ):
                continue
            outputs = mapping.get("outputs")
            if outputs is None or any(
                isinstance(o, dict) and o.get("output_id") == output_id
                for o in outputs
            ):
                places.append(stream)
                break
    return places, described


def check_refusals(goal: Goal, report: Report) -> None:
    """W7: a verified refusal's basis is in the records."""

    for stamp, path, item in goal.refusals():
        if not item.get("verified"):
            continue
        oid = str(item.get("observable_id") or "")
        node = str(item.get("blocked_node_id") or "")
        selector = str(item.get("selector") or "")
        basis = str(item.get("basis") or "")
        if node and "blocked_unsupported in this session's plan" in basis:
            here, _ = _plans_hold_blocked(goal, node, oid, path)
            anywhere, described = _plans_hold_blocked(goal, node, oid, None)
            if here:
                report.add(goal, "W7", "blocked_basis_in_its_plan", "ok")
            elif anywhere:
                report.add(goal, "W7", "blocked_basis_in_its_plan", "flag",
                           f"{oid}: node {node} recorded blocked only in another stream")
            elif not described:
                report.add(goal, "W7", "blocked_basis_in_its_plan", "insufficient",
                           f"{oid}: no record describes {node}'s outputs")
            else:
                report.add(goal, "W7", "blocked_basis_in_its_plan", "flag",
                           f"{oid}: no recorded plan holds {node} blocked with that output")
        if selector and "declares selector" in basis and "no program" in basis:
            read_later = []
            for later, _role, _p, payload in goal.of("result_quantities_extracted"):
                if later <= stamp:
                    continue
                bindings = payload.get("selector_bindings") or {}
                values = (
                    bindings.values() if isinstance(bindings, dict)
                    else [b.get("selector") if isinstance(b, dict) else b for b in bindings]
                )
                if selector in {str(v) for v in values}:
                    read_later.append(later[:19])
            report.add(goal, "W7", "refused_selector_not_read_later",
                       "flag" if read_later else "ok",
                       f"{oid}: selector {selector} read at {read_later[:3]}" if read_later else "")


def check_findings(goal: Goal, report: Report) -> None:
    """W8 and W9."""

    extractions: dict[str, dict] = {}
    for _stamp, _role, _path, payload in goal.of("result_quantities_extracted"):
        receipt = str(payload.get("receipt_sha256") or "")
        values = {}
        for quantity in (payload.get("record") or {}).get("quantities") or ():
            values[str(quantity.get("quantity_id"))] = quantity.get("value")
        extractions[receipt] = values

    def word_matches(word: dict) -> str:
        source = str(word.get("source_receipt_sha256") or "")
        if source not in extractions:
            return "insufficient"
        value = extractions[source].get(str(word.get("quantity_id") or ""))
        if value is None:
            return "insufficient"
        read = str(int(value)) if isinstance(value, (int, float)) and not isinstance(value, bool) and float(value).is_integer() else str(value)
        return "ok" if read.strip().casefold() == str(word.get("word")).strip().casefold() else "flag"

    for stamp, _role, _path, payload in goal.of("scientific_decision_recorded"):
        declared_then = {
            oid for oid, item in goal.declared.items() if item.get("_at", "") <= stamp
        }
        for finding in payload.get("findings") or ():
            fid = finding.get("finding_id")
            answers = str(finding.get("answers_observable_id") or "")
            operands = [
                side
                for relation in finding.get("relations") or ()
                for side in (relation.get("left") or {}, relation.get("right") or {})
                if side.get("claim_id")
            ]
            undeclared = [
                side for side in operands
                if str(side.get("claim_id")) not in declared_then
                and str(side.get("quantity_id") or "") not in declared_then
            ]
            expected = "answers" if answers else ("unrequested" if undeclared else "on_the_request")
            recorded = str(finding.get("standing") or "")
            if recorded == expected:
                report.add(goal, "W8", "standing_is_its_declarations", "ok")
            else:
                report.add(goal, "W8", "standing_is_its_declarations", "flag",
                           f"{fid}: recorded {recorded}, declarations say {expected} "
                           f"(undeclared operands {[s.get('claim_id') for s in undeclared]})")
            for relation in finding.get("relations") or ():
                left = relation.get("left") or {}
                right = relation.get("right") or {}
                op = relation.get("relation")
                lv = left.get("value")
                if isinstance(lv, str):
                    other = right.get("value", right.get("literal"))
                    if other is None:
                        report.add(goal, "W8", "relation_holds", "insufficient")
                        continue
                    equal = lv.strip().casefold() == str(other).strip().casefold()
                    truth = equal if op == "==" else (not equal if op == "!=" else None)
                else:
                    a = _number(lv)
                    b = _number(right.get("value_in_left_unit"))
                    if b is None and "literal" in right and not right.get("unit"):
                        b = _number(right.get("literal"))
                    if a is None or b is None:
                        report.add(goal, "W8", "relation_holds", "insufficient")
                        continue
                    truth = {"<": a < b, "<=": a <= b, ">": a > b, ">=": a >= b,
                             "==": a == b, "!=": a != b}.get(op)
                if truth is True:
                    report.add(goal, "W8", "relation_holds", "ok")
                else:
                    report.add(goal, "W8", "relation_holds", "flag",
                               f"{fid}: {left.get('claim_id')} {op} {right} read {truth}")
            for word in finding.get("answer") or ():
                verdict = word_matches(word)
                report.add(goal, "W9", "answer_is_the_word_read", verdict,
                           f"{fid}: {word.get('word')!r} vs receipt {str(word.get('source_receipt_sha256'))[:8]}"
                           if verdict == "flag" else "")
    for _stamp, _role, _path, payload in goal.completions():
        for oid, words in (payload.get("declared_categorical_answers") or {}).items():
            for word in words or ():
                verdict = word_matches(word)
                report.add(goal, "W9", "answer_is_the_word_read", verdict,
                           f"{oid}: {word.get('word')!r}" if verdict == "flag" else "")


def check_stationarity(goal: Goal, report: Report) -> None:
    """W10 and W11 (reader side)."""

    for _stamp, _role, _path, payload in goal.of("stationary_point_characterised"):
        record = payload.get("record") or {}
        if record.get("order_claimed") != record.get("observed_imaginary_modes"):
            report.add(goal, "W10", "order_is_its_mode_count", "flag", str(record)[:200])
        else:
            report.add(goal, "W10", "order_is_its_mode_count", "ok")
        report.add(goal, "W10", "receipt_names_its_stationarity",
                   "ok" if "stationarity" in record else "insufficient")
    verified: dict[str, tuple[str, str]] = {}
    for _stamp, _role, _path, payload in goal.of("program_result_verified"):
        record = payload.get("record") or {}
        status = str(payload.get("status") or record.get("state") or "")
        for artifact in record.get("output_artifacts") or ():
            verified[str(artifact.get("sha256") or "")] = (
                status, str(payload.get("node_id") or record.get("node_id") or ""))
    for _stamp, _role, _path, payload in goal.of("thermochemistry_derived"):
        source = str(payload.get("artifact_sha256") or (payload.get("record") or {}).get("artifact_sha256") or "")
        status, node = verified.get(source, ("", ""))
        if not status:
            report.add(goal, "W11", "free_energy_on_a_passing_result", "insufficient")
        elif status != "valid":
            report.add(goal, "W11", "free_energy_on_a_passing_result", "flag",
                       f"node {node} verified {status}; receipt {str(payload.get('receipt_sha256'))[:8]}")
        else:
            report.add(goal, "W11", "free_energy_on_a_passing_result", "ok")


def dissent(goal: Goal) -> list[dict]:
    """The dissent markers a goal holds, and whether its word names each."""

    text = " ".join(goal.reasons)
    markers: list[dict] = []

    def add(kind: str, name: str, extra: str = "") -> None:
        # A falsification is named by its own token, not by its id, which a
        # refusal or an earlier-cycle delivery line also prints.
        needle = (
            f"falsified_expectation:{name}"
            if kind in {"falsified_expectation", "falsified_diagnostic"}
            else name
        )
        markers.append({"kind": kind, "id": name, "named": bool(name) and needle in text,
                        "extra": extra})

    last_row: dict[str, dict] = {}
    for _stamp, _r, _p, payload in goal.completions():
        for row in payload.get("declared_observable_predictions") or ():
            oid = str(row.get("observable_id") or "")
            if oid and row.get("delivered_claim_id"):
                last_row[oid] = row
    for oid, row in sorted(last_row.items()):
        if row.get("agreement") == "diverged":
            add("falsified_diagnostic" if row.get("role") == "diagnostic" else "falsified_expectation",
                oid, "post-hoc" if row.get("declared_after_evidence") else "")
    seen_criteria = set()
    for _stamp, _r, _p, payload in goal.completions():
        for item in payload.get("anomaly_output_ids") or ():
            if str(item).startswith("failed_criterion:") and ":answered:" in str(item):
                label = str(item).split(":")[1]
                if label not in seen_criteria:
                    seen_criteria.add(label)
                    add("answered_criterion", label)
    for _stamp, _path, item in goal.refusals():
        add("verified_refusal" if item.get("verified") else "unverified_refusal",
            str(item.get("observable_id") or ""))
    for _stamp, _r, _p, payload in goal.of("scientific_decision_recorded"):
        for finding in payload.get("findings") or ():
            if finding.get("standing") == "unrequested":
                add("unrequested_finding", str(finding.get("finding_id") or ""))
            if finding.get("supersedes_finding_id"):
                add("superseded_finding", str(finding.get("finding_id") or ""),
                    str(finding.get("supersedes_finding_id")))
        for ref in (payload.get("record") or {}).get("evidence_refs") or ():
            if str(ref).startswith("doubt:"):
                add("doubt", str(ref)[6:14])
        for item in payload.get("menu_route_dispositions") or ():
            if str(item.get("disposition") or "") == "rejected":
                add("rejected_route", str(item.get("route_id") or item.get("route") or "")[:60])
    for oid in sorted(goal.retired()):
        add("superseded_observable", oid)
    return markers


def main() -> None:
    args = sys.argv[1:]
    out_dir = Path(args.pop(0))
    words: dict[tuple[str, str], dict] = {}
    if "--words" in args:
        index = args.index("--words")
        for line in Path(args[index + 1]).read_text().splitlines():
            result = json.loads(line)
            if result.get("replayed"):
                replayed = dict(result["replayed"])
                if replayed.get("kind") != "goal_settled" and result.get("path") == "run":
                    # The tree held the goal open: a recovery or reading row.
                    replayed["state"] = replayed.get("kind") or replayed.get("state")
                replayed["qualified"] = result.get("replayed_qualified") or []
                words[(result["agent"], result["goal_id"])] = replayed
        del args[index : index + 2]
    out_dir.mkdir(parents=True, exist_ok=True)
    report = Report()
    per_goal = []
    goals = discover(args)
    for item in goals:
        override = words.get((item["agent"], item["goal_id"]))
        if words and override is None:
            continue
        goal = Goal(Path(item["agent"]), item["goal_id"], override)
        raised = check_settlement(goal, report)
        check_qualified(goal, report, raised)
        check_executor(goal, report)
        check_completions(goal, report)
        check_refusals(goal, report)
        check_findings(goal, report)
        check_stationarity(goal, report)
        per_goal.append({
            **item,
            "word": goal.word,
            "streams": len(goal.streams),
            "inferred_streams": goal.inferred,
            "declared": len(goal.declared),
            "dissent": dissent(goal),
        })
    (out_dir / "goals.json").write_text(json.dumps(per_goal, indent=1))
    (out_dir / "flags.json").write_text(json.dumps(report.flags, indent=1))
    summary = []
    for (cls, check), counter in sorted(report.counts.items()):
        summary.append(f"{cls:4} {check:42} ok {counter['ok']:4} flag {counter['flag']:4} insufficient {counter['insufficient']:4}")
    words_hist = collections.Counter(g["word"] or "(unsettled)" for g in per_goal)
    summary.append("words: " + ", ".join(f"{k} {v}" for k, v in sorted(words_hist.items())))
    summary.append(f"goals read: {len(per_goal)}; with inferred session streams: "
                   f"{sum(1 for g in per_goal if g['inferred_streams'])}")
    (out_dir / "summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
