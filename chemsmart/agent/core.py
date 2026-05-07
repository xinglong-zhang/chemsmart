from __future__ import annotations

import importlib
import json
import os
import re
import uuid
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

from chemsmart.agent.providers import get_provider
from chemsmart.agent.registry import ToolRegistry
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.settings import ORCAJobSettings

_RISKY_TOOLS = {"run_local", "submit_hpc"}
_GAUSSIAN_ROUTE_RE = re.compile(r"^\s*#\s*\S+", re.MULTILINE)
_ORCA_ROUTE_RE = re.compile(r"^\s*!\s*\S+", re.MULTILINE)
_REFERENCE_RE = re.compile(
    r"^\$step(?P<index>\d+)(?P<path>(?:\.[A-Za-z_][A-Za-z0-9_]*)*)$"
)


class Step(BaseModel):
    tool: str
    args: dict[str, Any] = Field(default_factory=dict)
    rationale: str = ""


class Plan(BaseModel):
    steps: list[Step]
    rationale: str = ""
    estimated_cost: str | None = None


class CriticVerdict(BaseModel):
    verdict: Literal["ok", "warn", "reject"]
    issues: list[str] = Field(default_factory=list)
    rationale: str = ""


class SessionState(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    session_id: str
    cwd: str
    current_step_index: int = 0
    plan: Plan | None = None
    request: str | None = None
    env_snapshot: dict[str, str | None] = Field(default_factory=dict)

    def save(self, path: Path) -> None:
        path.write_text(self.model_dump_json(indent=2))

    @classmethod
    def load(cls, path: Path) -> "SessionState":
        return cls.model_validate_json(path.read_text())


class DecisionLog:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def write(
        self,
        kind: str,
        payload: dict[str, Any],
        rationale: str = "",
    ) -> None:
        entry = {
            "ts": datetime.now(UTC).isoformat(),
            "kind": kind,
            "payload": _json_safe(payload),
            "rationale": rationale,
        }
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(entry, sort_keys=True) + "\n")

    def read_all(self) -> list[dict[str, Any]]:
        if not self.path.exists():
            return []
        with self.path.open(encoding="utf-8") as handle:
            return [json.loads(line) for line in handle if line.strip()]


class AgentSession:
    def __init__(
        self,
        provider: Any | None = None,
        registry: ToolRegistry | None = None,
        session_root: str | os.PathLike[str] | None = None,
        transport: Any | None = None,
    ) -> None:
        self._provider = provider
        self.registry = registry or ToolRegistry.default()
        self.transport = transport
        self.session_root = Path(session_root or _default_session_root())
        self.session_root.mkdir(parents=True, exist_ok=True)
        self.state: SessionState | None = None
        self.session_dir: Path | None = None
        self.decision_log: DecisionLog | None = None

    @classmethod
    def resume(
        cls,
        session_id: str,
        **kwargs: Any,
    ) -> dict[str, Any]:
        session = cls(**_session_kwargs(kwargs))
        session._load_existing_session(session_id)
        assert session.state is not None
        assert session.decision_log is not None
        original_cwd = os.getcwd()
        os.chdir(session.state.cwd)
        try:
            completed_results = session._load_completed_results()
            return session._continue_run(
                completed_results=completed_results,
                dry_submit=kwargs.get("dry_submit", True),
                allow_remote_unknown=kwargs.get("allow_remote_unknown", False),
                allow_critic_override=kwargs.get(
                    "allow_critic_override", False
                ),
                rerender_plan=False,
            )
        finally:
            os.chdir(original_cwd)

    def run(
        self,
        request: str,
        dry_submit: bool = True,
        allow_remote_unknown: bool = False,
        allow_critic_override: bool = False,
    ) -> dict[str, Any]:
        session_id = _new_session_id()
        self.session_dir = self.session_root / session_id
        self.session_dir.mkdir(parents=True, exist_ok=True)
        self.decision_log = DecisionLog(
            self.session_dir / "decision_log.jsonl"
        )
        self.state = SessionState(
            session_id=session_id,
            cwd=os.path.abspath(os.getcwd()),
            current_step_index=0,
            request=request,
            env_snapshot=_env_snapshot(),
        )
        self._save_state()
        self.decision_log.write(
            "request", {"request": request}, rationale=request
        )

        plan = self._planner_call(request)
        self.state.plan = plan
        self._save_state()
        self.decision_log.write(
            "plan", plan.model_dump(), rationale=plan.rationale
        )

        return self._continue_run(
            completed_results=[],
            dry_submit=dry_submit,
            allow_remote_unknown=allow_remote_unknown,
            allow_critic_override=allow_critic_override,
            rerender_plan=True,
        )

    def _continue_run(
        self,
        completed_results: list[Any],
        dry_submit: bool,
        allow_remote_unknown: bool,
        allow_critic_override: bool,
        rerender_plan: bool,
    ) -> dict[str, Any]:
        assert self.state is not None
        assert self.state.plan is not None
        assert self.session_dir is not None
        assert self.decision_log is not None

        plan = self.state.plan
        results = list(completed_results)
        dry_run_result = self._find_prior_result(results, "dry_run_input")
        runtime_result = self._find_prior_result(results, "validate_runtime")

        risky_start = len(plan.steps)
        for step_index in range(
            self.state.current_step_index, len(plan.steps)
        ):
            step = plan.steps[step_index]
            if step.tool in _RISKY_TOOLS:
                risky_start = step_index
                break
            result = self._execute_step(step_index, step, results)
            results.append(result)
            self.state.current_step_index = step_index + 1
            self._save_state()
            if step.tool == "dry_run_input":
                dry_run_result = result
            elif step.tool == "validate_runtime":
                runtime_result = result

        preview_submit = None
        if (
            risky_start < len(plan.steps)
            and plan.steps[risky_start].tool == "submit_hpc"
        ):
            preview_submit = self._preview_submit_step(
                risky_start,
                plan.steps[risky_start],
                results,
            )

        verdict = self._get_logged_verdict() or self._critic_call(
            plan=plan,
            dry_run_result=dry_run_result,
        )
        verdict = self._apply_deterministic_gates(
            verdict=verdict,
            runtime_result=runtime_result,
            dry_run_result=dry_run_result,
            preview_submit=preview_submit,
        )
        self.decision_log.write(
            "critic_verdict",
            verdict.model_dump(),
            rationale=verdict.rationale,
        )

        should_block = self._should_block(
            verdict=verdict,
            allow_remote_unknown=allow_remote_unknown,
            allow_critic_override=allow_critic_override,
        )
        if should_block:
            return {
                "session_id": self.state.session_id,
                "session_dir": str(self.session_dir),
                "plan": plan,
                "plan_text": render_plan(plan) if rerender_plan else None,
                "critic_verdict": verdict,
                "completed_steps": self.state.current_step_index,
                "blocked": True,
                "dry_run_result": dry_run_result,
                "runtime_result": runtime_result,
                "preview_submit": preview_submit,
            }

        for step_index in range(risky_start, len(plan.steps)):
            if step_index < self.state.current_step_index:
                continue
            step = plan.steps[step_index]
            if (
                step.tool == "submit_hpc"
                and preview_submit is not None
                and dry_submit
            ):
                result = preview_submit
                self._record_final_preview_step(step_index, step, result)
            else:
                extra_kwargs = {}
                if step.tool == "submit_hpc":
                    extra_kwargs["execute"] = not dry_submit
                    if self.transport is not None:
                        extra_kwargs["transport"] = self.transport
                result = self._execute_step(
                    step_index,
                    step,
                    results,
                    extra_kwargs=extra_kwargs,
                )
            results.append(result)
            self.state.current_step_index = step_index + 1
            self._save_state()

        return {
            "session_id": self.state.session_id,
            "session_dir": str(self.session_dir),
            "plan": plan,
            "plan_text": render_plan(plan) if rerender_plan else None,
            "critic_verdict": verdict,
            "completed_steps": self.state.current_step_index,
            "blocked": False,
            "dry_run_result": dry_run_result,
            "runtime_result": runtime_result,
            "preview_submit": preview_submit,
            "results": results,
        }

    def _provider_instance(self) -> Any:
        if self._provider is None:
            self._provider = get_provider()
        return self._provider

    def _planner_call(self, request: str) -> Plan:
        provider = self._provider_instance()
        prompt = _load_prompt("planner.md")
        tool_defs = self.registry.openai_tool_defs()
        response = provider.chat(
            [
                {"role": "system", "content": prompt},
                {
                    "role": "user",
                    "content": json.dumps(
                        {
                            "request": request,
                            "tools": tool_defs,
                        }
                    ),
                },
            ],
            tools=None,
        )
        return Plan.model_validate(_parse_json_response(response))

    def _critic_call(
        self,
        plan: Plan,
        dry_run_result: dict[str, Any] | None,
    ) -> CriticVerdict:
        provider = self._provider_instance()
        prompt = _load_prompt("critic.md")
        response = provider.chat(
            [
                {"role": "system", "content": prompt},
                {
                    "role": "user",
                    "content": json.dumps(
                        {
                            "plan": plan.model_dump(),
                            "dry_run_input": dry_run_result,
                        }
                    ),
                },
            ],
            tools=None,
        )
        return CriticVerdict.model_validate(_parse_json_response(response))

    def _execute_step(
        self,
        step_index: int,
        step: Step,
        prior_results: list[Any],
        extra_kwargs: dict[str, Any] | None = None,
    ) -> Any:
        assert self.session_dir is not None
        assert self.decision_log is not None
        resolved_args = _resolve_refs(step.args, prior_results)
        if extra_kwargs:
            resolved_args.update(extra_kwargs)
        self.decision_log.write(
            "tool_call",
            {
                "step_index": step_index,
                "tool": step.tool,
                "args": _preview_value(resolved_args),
            },
            rationale=step.rationale,
        )
        result = self.registry.call(step.tool, resolved_args)
        if _is_tool_error(result):
            raise RuntimeError(result["error"]["message"])
        artifact_path = self._write_result_artifact(step_index, result)
        self.decision_log.write(
            "tool_result",
            {
                "step_index": step_index,
                "tool": step.tool,
                "artifact": artifact_path.name,
                "payload": _preview_value(result),
            },
            rationale=step.rationale,
        )
        return result

    def _preview_submit_step(
        self,
        step_index: int,
        step: Step,
        prior_results: list[Any],
    ) -> Any:
        assert self.decision_log is not None
        resolved_args = _resolve_refs(step.args, prior_results)
        resolved_args["execute"] = False
        self.decision_log.write(
            "tool_preview",
            {
                "step_index": step_index,
                "tool": step.tool,
                "args": _preview_value(resolved_args),
            },
            rationale=step.rationale,
        )
        result = self.registry.call(step.tool, resolved_args)
        if _is_tool_error(result):
            raise RuntimeError(result["error"]["message"])
        self.decision_log.write(
            "tool_preview_result",
            {
                "step_index": step_index,
                "tool": step.tool,
                "payload": _preview_value(result),
            },
            rationale=step.rationale,
        )
        return result

    def _record_final_preview_step(
        self,
        step_index: int,
        step: Step,
        result: Any,
    ) -> None:
        assert self.decision_log is not None
        artifact_path = self._write_result_artifact(step_index, result)
        self.decision_log.write(
            "tool_result",
            {
                "step_index": step_index,
                "tool": step.tool,
                "artifact": artifact_path.name,
                "payload": _preview_value(result),
                "from_preview": True,
            },
            rationale=step.rationale,
        )

    def _write_result_artifact(self, step_index: int, result: Any) -> Path:
        assert self.session_dir is not None
        artifact_path = self.session_dir / f"step_{step_index + 1:02d}.json"
        artifact_path.write_text(
            json.dumps(_json_safe(result), indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return artifact_path

    def _load_existing_session(self, session_id: str) -> None:
        self.session_dir = self.session_root / session_id
        state_path = self.session_dir / "session.json"
        self.state = SessionState.load(state_path)
        self.decision_log = DecisionLog(
            self.session_dir / "decision_log.jsonl"
        )

    def _load_completed_results(self) -> list[Any]:
        assert self.session_dir is not None
        assert self.state is not None
        results: list[Any] = []
        for step_index in range(self.state.current_step_index):
            artifact_path = self.session_dir / (
                f"step_{step_index + 1:02d}.json"
            )
            with artifact_path.open(encoding="utf-8") as handle:
                results.append(json.load(handle))
        return results

    def _save_state(self) -> None:
        assert self.session_dir is not None
        assert self.state is not None
        self.state.save(self.session_dir / "session.json")

    def _find_prior_result(self, results: list[Any], tool_name: str) -> Any:
        assert self.state is not None
        assert self.state.plan is not None
        return next(
            (
                result
                for step, result in zip(self.state.plan.steps, results)
                if step.tool == tool_name
            ),
            None,
        )

    def _get_logged_verdict(self) -> CriticVerdict | None:
        assert self.decision_log is not None
        verdict_entries = [
            entry
            for entry in self.decision_log.read_all()
            if entry.get("kind") == "critic_verdict"
        ]
        if not verdict_entries:
            return None
        return CriticVerdict.model_validate(verdict_entries[-1]["payload"])

    def _apply_deterministic_gates(
        self,
        verdict: CriticVerdict,
        runtime_result: dict[str, Any] | None,
        dry_run_result: dict[str, Any] | None,
        preview_submit: dict[str, Any] | None,
    ) -> CriticVerdict:
        issues = list(verdict.issues)
        rationale_parts = [verdict.rationale] if verdict.rationale else []
        final_verdict = verdict.verdict

        if runtime_result is not None:
            runtime_ok = runtime_result.get("ok")
            if runtime_ok == "fail":
                final_verdict = "reject"
                issues.extend(runtime_result.get("local_issues", []))
                rationale_parts.append("validate_runtime returned fail")
            elif runtime_ok == "partial":
                final_verdict = (
                    "warn" if final_verdict == "ok" else final_verdict
                )
                issues.extend(runtime_result.get("remote_unknown", []))
                rationale_parts.append("validate_runtime returned partial")

        malformed_issue = _malformed_input_issue(dry_run_result)
        if malformed_issue is not None and final_verdict != "reject":
            final_verdict = "warn"
            issues.append(malformed_issue)
            rationale_parts.append(
                "dry-run input route line failed basic validation"
            )

        if preview_submit is not None:
            duplicate_check = preview_submit.get("duplicate_check", {})
            if duplicate_check.get("duplicate"):
                final_verdict = "reject"
                message = (
                    duplicate_check.get("message")
                    or "duplicate submission detected"
                )
                issues.append(message)
                rationale_parts.append(
                    "submit_hpc duplicate check rejected the plan"
                )

        return CriticVerdict(
            verdict=final_verdict,
            issues=_dedupe_strings(issues),
            rationale="; ".join(part for part in rationale_parts if part),
        )

    @staticmethod
    def _should_block(
        verdict: CriticVerdict,
        allow_remote_unknown: bool,
        allow_critic_override: bool,
    ) -> bool:
        if verdict.verdict == "reject":
            return True
        if verdict.verdict == "ok":
            return False

        has_remote_unknown = any(
            issue.startswith("server.")
            or issue.endswith("on HPC")
            or issue == "ssh login reachable"
            for issue in verdict.issues
        )
        has_other_warn = any(
            not (
                issue.startswith("server.")
                or issue.endswith("on HPC")
                or issue == "ssh login reachable"
            )
            for issue in verdict.issues
        )
        if not verdict.issues:
            has_other_warn = True
        if has_remote_unknown and not allow_remote_unknown:
            return True
        if has_other_warn and not allow_critic_override:
            return True
        return False


def render_plan(plan: Plan) -> str:
    lines = ["Plan:"]
    if plan.rationale:
        lines.append(f"Rationale: {plan.rationale}")
    if plan.estimated_cost:
        lines.append(f"Estimated cost: {plan.estimated_cost}")
    for index, step in enumerate(plan.steps, start=1):
        lines.append(
            f"{index}. {step.tool} {json.dumps(_preview_value(step.args), sort_keys=True)}"
        )
        if step.rationale:
            lines.append(f"   - {step.rationale}")
    return "\n".join(lines)


def run_agent(request: str, **kwargs: Any) -> dict[str, Any]:
    return AgentSession(**_session_kwargs(kwargs)).run(
        request,
        dry_submit=kwargs.get("dry_submit", True),
        allow_remote_unknown=kwargs.get("allow_remote_unknown", False),
        allow_critic_override=kwargs.get("allow_critic_override", False),
    )


def _session_kwargs(kwargs: dict[str, Any]) -> dict[str, Any]:
    return {
        key: kwargs[key]
        for key in ("provider", "registry", "session_root", "transport")
        if key in kwargs
    }


def _default_session_root() -> str:
    return str(Path.home() / ".chemsmart" / "agent" / "sessions")


def _new_session_id() -> str:
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}-{uuid.uuid4().hex[:8]}"


def _env_snapshot() -> dict[str, str | None]:
    return {
        "AI_PROVIDER": os.environ.get("AI_PROVIDER"),
        "PWD": os.path.abspath(os.getcwd()),
    }


def _load_prompt(name: str) -> str:
    prompt_path = Path(__file__).with_name("prompts") / name
    return prompt_path.read_text(encoding="utf-8")


def _parse_json_response(response: Any) -> dict[str, Any]:
    if (
        isinstance(response, dict)
        and "parsed" in response
        and isinstance(response["parsed"], dict)
    ):
        return response["parsed"]
    if (
        isinstance(response, dict)
        and "json" in response
        and isinstance(response["json"], dict)
    ):
        return response["json"]
    if (
        isinstance(response, dict)
        and "content" in response
        and isinstance(response["content"], dict)
    ):
        return response["content"]
    text = _extract_text(response)
    text = text.strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"LLM returned invalid JSON: {exc}\nRaw: {text[:300]!r}"
        ) from exc


def _extract_text(response: Any) -> str:
    if isinstance(response, str):
        return response
    if isinstance(response, dict):
        content = response.get("content")
        if isinstance(content, str):
            return content
        if isinstance(content, list):
            parts: list[str] = []
            for item in content:
                if isinstance(item, dict) and item.get("type") == "text":
                    parts.append(item.get("text", ""))
            if parts:
                return "\n".join(parts)
        choices = response.get("choices") or []
        if choices:
            message = choices[0].get("message", {})
            content = message.get("content", "")
            if isinstance(content, str):
                return content
            if isinstance(content, list):
                return "\n".join(
                    item.get("text", "")
                    for item in content
                    if isinstance(item, dict)
                )
    raise ValueError("Could not extract text from provider response")


def _resolve_refs(value: Any, prior_results: list[Any]) -> Any:
    if isinstance(value, str):
        return _resolve_ref_string(value, prior_results)
    if isinstance(value, list):
        return [_resolve_refs(item, prior_results) for item in value]
    if isinstance(value, dict):
        return {
            key: _resolve_refs(item, prior_results)
            for key, item in value.items()
        }
    return value


def _resolve_ref_string(value: str, prior_results: list[Any]) -> Any:
    match = _REFERENCE_RE.match(value)
    if match is None:
        return value
    index = int(match.group("index")) - 1
    if index < 0 or index >= len(prior_results):
        raise IndexError(f"Step reference {value!r} is out of range")
    resolved = prior_results[index]
    for part in match.group("path").split("."):
        if not part:
            continue
        if isinstance(resolved, dict):
            resolved = resolved[part]
        else:
            resolved = getattr(resolved, part)
    return _restore_json_result(resolved)


def _json_safe(value: Any) -> Any:
    if isinstance(value, Molecule):
        positions = value.positions
        if hasattr(positions, "tolist"):
            positions = positions.tolist()
        return {
            "__chemsmart_type__": "molecule",
            "symbols": _json_safe(list(value.symbols)),
            "positions": _json_safe(positions),
            "charge": value.charge,
            "multiplicity": value.multiplicity,
            "frozen_atoms": _json_safe(value.frozen_atoms),
            "pbc_conditions": _json_safe(value.pbc_conditions),
            "translation_vectors": _json_safe(value.translation_vectors),
            "energy": value.energy,
            "forces": _json_safe(value.forces),
            "velocities": _json_safe(value.velocities),
            "info": _json_safe(value.info),
        }
    if isinstance(value, (GaussianJobSettings, ORCAJobSettings)):
        return {
            "__chemsmart_type__": "settings",
            "module": value.__class__.__module__,
            "class": value.__class__.__name__,
            **{
                key: _json_safe(item)
                for key, item in vars(value).items()
                if not key.startswith("_")
            },
        }
    if isinstance(value, Job):
        return {
            "__chemsmart_type__": "job",
            "module": value.__class__.__module__,
            "class": value.__class__.__name__,
            "molecule": _json_safe(value.molecule),
            "settings": _json_safe(value.settings),
            "label": value.label,
            "folder": value.folder,
        }
    if isinstance(value, BaseModel):
        return _json_safe(value.model_dump())
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, bytes):
        return {"type": "bytes", "length": len(value)}
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return {"type": value.__class__.__name__, "repr": repr(value)}


def _preview_value(value: Any) -> Any:
    return _json_safe(value)


def _restore_json_result(value: Any) -> Any:
    if isinstance(value, list):
        return [_restore_json_result(item) for item in value]
    if not isinstance(value, dict):
        return value

    marker = value.get("__chemsmart_type__")
    if marker == "molecule":
        return Molecule(
            symbols=value.get("symbols"),
            positions=value.get("positions"),
            charge=value.get("charge"),
            multiplicity=value.get("multiplicity"),
            frozen_atoms=value.get("frozen_atoms"),
            pbc_conditions=value.get("pbc_conditions"),
            translation_vectors=value.get("translation_vectors"),
            energy=value.get("energy"),
            forces=value.get("forces"),
            velocities=value.get("velocities"),
            info=value.get("info"),
        )

    if marker == "settings":
        settings_cls = _load_class(value["module"], value["class"])
        kwargs = {
            key: _restore_json_result(item)
            for key, item in value.items()
            if key not in {"__chemsmart_type__", "module", "class"}
        }
        return settings_cls(**kwargs)

    if marker == "job":
        job_cls = _load_class(value["module"], value["class"])
        job = job_cls(
            molecule=_restore_json_result(value["molecule"]),
            settings=_restore_json_result(value["settings"]),
            label=value["label"],
            jobrunner=None,
        )
        if value.get("folder"):
            job.set_folder(value["folder"])
        return job

    return {key: _restore_json_result(item) for key, item in value.items()}


def _load_class(module_name: str, class_name: str) -> type[Any]:
    module = importlib.import_module(module_name)
    return getattr(module, class_name)


def _is_tool_error(result: Any) -> bool:
    return (
        isinstance(result, dict)
        and result.get("ok") is False
        and "error" in result
    )


def _malformed_input_issue(
    dry_run_result: dict[str, Any] | None,
) -> str | None:
    if not dry_run_result:
        return None
    content = dry_run_result.get("content")
    inputfile = str(dry_run_result.get("inputfile", "")).lower()
    if not isinstance(content, str):
        return None
    if inputfile.endswith(".inp"):
        if _ORCA_ROUTE_RE.search(content) is None:
            return "ORCA route line missing or malformed"
        return None
    if _GAUSSIAN_ROUTE_RE.search(content) is None:
        return "Gaussian route line missing or malformed"
    return None


def _dedupe_strings(values: list[str]) -> list[str]:
    deduped: list[str] = []
    seen: set[str] = set()
    for value in values:
        if value not in seen:
            seen.add(value)
            deduped.append(value)
    return deduped
