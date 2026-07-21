from __future__ import annotations

from chemsmart.agent.harness.runner import evaluate_harness


def test_evaluate_harness_accepts_dict_plan_for_gaussian_ts():
    plan = {
        "steps": [
            {"tool": "build_job", "args": {"kind": "gaussian.ts"}},
            {"tool": "dry_run_input", "args": {"job": "$step1"}},
        ]
    }
    dry_run_results = [
        {
            "inputfile": "ts.com",
            "command": "chemsmart run gaussian -f ts.xyz -c 0 -m 1 ts",
            "cli_grounded": True,
            "content": "%chk=ts.chk\n# b3lyp/6-31g* opt=(ts,calcfc,noeigentest)\n\n",
        }
    ]

    result = evaluate_harness(plan, dry_run_results)

    assert result.verdict == "ok"
    assert [r.rule_id for r in result.rule_results] == [
        "cli.grounding",
        "gaussian.ts.route",
    ]


def test_evaluate_harness_rejects_dry_run_without_cli_command():
    plan = {
        "steps": [
            {"tool": "build_job", "args": {"kind": "gaussian.opt"}},
            {"tool": "dry_run_input", "args": {"job": "$step1"}},
        ]
    }

    result = evaluate_harness(
        plan,
        [{"inputfile": "opt.com", "content": "# opt b3lyp/6-31g*\n"}],
    )

    assert result.verdict == "reject"
    assert result.failed_rule_ids == ["cli.grounding.missing"]
