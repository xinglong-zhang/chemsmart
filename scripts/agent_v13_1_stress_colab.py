"""Colab stress test for the v13.1 compact-SPEC adapter.

Inputs expected on the Colab VM:
- /content/chemsmart_worktree.tar.gz
- /content/eval_representative_v131.jsonl

Writes:
- /content/chemsmart_v131_stress_result.json
"""

from __future__ import annotations

import json
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT_TAR = Path("/content/chemsmart_worktree.tar.gz")
SRC = Path("/content/chemsmart_src")
EVAL = Path("/content/eval_representative_v131.jsonl")
OUT = Path("/content/chemsmart_v131_stress_result.json")


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", shlex.join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def setup_repo() -> None:
    if SRC.exists():
        shutil.rmtree(SRC)
    SRC.mkdir(parents=True)
    with tarfile.open(ROOT_TAR, "r:gz") as archive:
        archive.extractall(SRC)
    run([sys.executable, "-m", "pip", "install", "-q", "-e", "."], cwd=SRC)
    run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "-q",
            "--force-reinstall",
            "--no-deps",
            "numpy==1.26.4",
            "scipy==1.11.4",
            "ase==3.22.1",
        ]
    )


def targeted_cases() -> list[dict]:
    return [
        {
            "name": "gaussian_ts_canonical_route_dropped",
            "spec": {
                "intent": "workflow",
                "jobs": [
                    {
                        "id": 1,
                        "kind": "gaussian.ts",
                        "file": "ts.xyz",
                        "charge": 0,
                        "mult": 1,
                        "settings": {
                            "additional_opt_options_in_route": [
                                "ts",
                                "calcfc/noeigentest",
                            ]
                        },
                    }
                ],
            },
            "must_not_contain": ["[", "]", "calcfc/noeigentest"],
        },
        {
            "name": "gaussian_opt_freq_route_param",
            "spec": {
                "intent": "workflow",
                "jobs": [
                    {
                        "id": 1,
                        "kind": "gaussian.opt",
                        "file": "opt.xyz",
                        "charge": 0,
                        "mult": 1,
                        "settings": {"freq": True},
                    }
                ],
            },
            "must_contain": ["--additional-route-parameters freq"],
            "must_not_contain": ["--freq", "opt_freq"],
        },
        {
            "name": "dias_nested_fragment_first_only",
            "spec": {
                "intent": "workflow",
                "jobs": [
                    {
                        "id": 1,
                        "kind": "gaussian.dias",
                        "file": "dimer.xyz",
                        "charge": 0,
                        "mult": 1,
                        "settings": {
                            "fragment_indices": [[1, 2, 3], [4, 5, 6]]
                        },
                    }
                ],
            },
            "must_contain": ["--fragment-indices 1,2,3"],
            "must_not_contain": ["1 2 3", "4,5,6 dias"],
        },
        {
            "name": "qmmm_atom_ranges_comma_rendered",
            "spec": {
                "intent": "workflow",
                "jobs": [
                    {
                        "id": 1,
                        "kind": "gaussian.qmmm",
                        "file": "qm.xyz",
                        "charge": 0,
                        "mult": 1,
                        "settings": {
                            "high_level_atoms": [1, 2, 3],
                            "low_level_atoms": [4, 5, 6],
                        },
                    }
                ],
            },
            "must_contain": [
                "--high-level-atoms 1,2,3",
                "--low-level-atoms 4,5,6",
            ],
        },
        {
            "name": "gaussian_opt_to_orca_sp_chain",
            "spec": {
                "intent": "workflow",
                "jobs": [
                    {
                        "id": 1,
                        "kind": "gaussian.opt",
                        "file": "chain.xyz",
                        "charge": 0,
                        "mult": 1,
                    },
                    {
                        "id": 2,
                        "kind": "orca.sp",
                        "geom_from": 1,
                        "charge": 0,
                        "mult": 1,
                        "execution": "submit",
                        "server": "chemnode1",
                    },
                ],
            },
            "command_count": 2,
            "must_contain": ["chemsmart sub -s chemnode1 orca"],
        },
    ]


def run_targeted() -> list[dict]:
    from chemsmart.agent import v8_adapter
    from chemsmart.agent.synthesis import _normalize_v8_spec

    results = []
    for case in targeted_cases():
        out = v8_adapter.adapt(case["spec"], validate=True)
        commands = out.get("commands") or []
        joined = "\n".join(commands)
        ok = bool(out.get("valid")) and bool(commands)
        ok = ok and len(commands) == case.get("command_count", 1)
        for needle in case.get("must_contain", []):
            ok = ok and needle in joined
        for needle in case.get("must_not_contain", []):
            ok = ok and needle not in joined
        normalized = _normalize_v8_spec(case["spec"])
        ok = ok and normalized["status"] == "ready"
        results.append(
            {
                "name": case["name"],
                "ok": ok,
                "commands": commands,
                "errors": out.get("errors") or [],
            }
        )
    return results


def run_eval_stress() -> dict:
    from chemsmart.agent import v8_adapter

    total = 0
    workflow = 0
    nonworkflow = 0
    valid = 0
    failures: list[dict] = []
    for line in EVAL.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        total += 1
        record = json.loads(line)
        content = record["messages"][2]["content"]
        spec = json.loads(content)
        out = v8_adapter.adapt(spec, validate=True)
        if spec.get("intent") == "workflow":
            workflow += 1
            if out.get("valid") and out.get("commands"):
                valid += 1
            elif len(failures) < 20:
                failures.append(
                    {
                        "id": record.get("id"),
                        "intent": spec.get("intent"),
                        "errors": out.get("errors") or [],
                        "commands": out.get("commands") or [],
                    }
                )
        else:
            nonworkflow += 1
            if out.get("valid") is True and not out.get("commands"):
                valid += 1
            elif len(failures) < 20:
                failures.append(
                    {
                        "id": record.get("id"),
                        "intent": spec.get("intent"),
                        "errors": out.get("errors") or [],
                    }
                )
    return {
        "total": total,
        "workflow": workflow,
        "nonworkflow": nonworkflow,
        "valid": valid,
        "failures": failures,
    }


def main() -> int:
    setup_repo()
    os.chdir(SRC)
    result = {
        "targeted": run_targeted(),
        "eval_stress": run_eval_stress(),
    }
    target_ok = all(item["ok"] for item in result["targeted"])
    eval_ok = result["eval_stress"]["valid"] == result["eval_stress"]["total"]
    result["verdict"] = "PASS" if target_ok and eval_ok else "FAIL"
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2), flush=True)
    return 0 if result["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
