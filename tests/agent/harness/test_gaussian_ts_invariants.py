from __future__ import annotations

from chemsmart.agent.harness.invariants.gaussian_ts import (
    check_gaussian_ts_route,
)


def test_canonical_gaussian_ts_route_passes():
    result = check_gaussian_ts_route(
        "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest)"
    )

    assert result.verdict == "ok"
    assert result.issues == ()


def test_python_list_leak_rejects():
    result = check_gaussian_ts_route(
        "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest,['ts'])"
    )

    assert result.verdict == "reject"
    assert any(
        "Python list/string leak" in issue.message for issue in result.issues
    )


def test_duplicate_ts_keyword_rejects():
    result = check_gaussian_ts_route(
        "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest,ts)"
    )

    assert result.verdict == "reject"
    assert any("duplicate" in issue.message for issue in result.issues)


def test_calcfc_and_calcall_together_rejects():
    result = check_gaussian_ts_route(
        "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest,calcall)"
    )

    assert result.verdict == "reject"
    assert any(
        "calcfc and calcall" in issue.message for issue in result.issues
    )


def test_missing_noeigentest_rejects():
    result = check_gaussian_ts_route("# b3lyp/6-31g* opt=(ts,calcfc)")

    assert result.verdict == "reject"
    assert any("noeigentest" in issue.message for issue in result.issues)


def test_true_extra_route_option_passes():
    result = check_gaussian_ts_route(
        "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest,maxstep=8)"
    )

    assert result.verdict == "ok"


def test_calcall_variant_passes_without_calcfc():
    result = check_gaussian_ts_route(
        "# b3lyp/6-31g* opt=(ts,noeigentest,calcall)"
    )

    assert result.verdict == "ok"
