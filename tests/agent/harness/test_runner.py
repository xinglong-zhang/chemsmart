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
            "content": "%chk=ts.chk\n# b3lyp/6-31g* opt=(ts,calcfc,noeigentest)\n\n",
        }
    ]

    result = evaluate_harness(plan, dry_run_results)

    assert result.verdict == "ok"
    assert len(result.rule_results) == 1
    assert result.rule_results[0].rule_id == "gaussian.ts.route"
