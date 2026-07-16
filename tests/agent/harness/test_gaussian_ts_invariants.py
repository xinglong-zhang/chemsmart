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


def test_mangled_ts_extra_rejects():
    result = check_gaussian_ts_route(
        "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest,calcfc/noeigentest)"
    )

    assert result.verdict == "reject"
    assert any("non-allowlisted" in issue.message for issue in result.issues)


def test_tssearch_type_extra_rejects():
    result = check_gaussian_ts_route(
        "# b3lyp/6-31g* opt=(ts,calcfc,noeigentest,tssearch_type=1)"
    )

    assert result.verdict == "reject"
    assert any("non-allowlisted" in issue.message for issue in result.issues)


def test_calcall_variant_passes_without_calcfc():
    result = check_gaussian_ts_route(
        "# b3lyp/6-31g* opt=(ts,noeigentest,calcall)"
    )

    assert result.verdict == "ok"


def test_duplicate_opt_block_via_route_params_rejects():
    # Leak: a second opt=(...) block from --additional-route-parameters.
    # The first block is canonical, so the old first-block-only check missed it.
    result = check_gaussian_ts_route(
        "# opt=(ts,calcfc,noeigentest) m062x 6-31G(d) "
        "opt=(ts,calcfc,noeigentest)"
    )

    assert result.verdict == "reject"
    assert any(
        "more than one opt=(...) block" in issue.message
        for issue in result.issues
    )


def test_second_parenthesized_opt_form_rejects():
    result = check_gaussian_ts_route(
        "# opt=(ts,calcfc,noeigentest) freq b3lyp/6-31g(d) "
        "opt(maxstep=8)"
    )

    assert result.verdict == "reject"
    assert any(
        "more than one Opt route form" in issue.message
        for issue in result.issues
    )


def test_bare_ts_token_outside_opt_block_rejects():
    # Leak: bare `ts` route token from --additional-route-parameters, outside
    # the runtime-owned opt=(...) block.
    result = check_gaussian_ts_route(
        "# opt=(ts,calcfc,noeigentest) m062x 6-31G(d) ts"
    )

    assert result.verdict == "reject"
    assert any(
        "outside the opt=(...) block" in issue.message
        for issue in result.issues
    )


def test_bare_calcfc_noeigentest_outside_opt_block_rejects():
    result = check_gaussian_ts_route(
        "# opt=(ts,calcfc,noeigentest) m062x 6-31G(d) calcfc noeigentest"
    )

    assert result.verdict == "reject"
    leaked = [
        issue
        for issue in result.issues
        if "outside the opt=(...) block" in issue.message
    ]
    assert leaked
    assert set(leaked[0].evidence["leaked_tokens"]) == {
        "calcfc",
        "noeigentest",
    }


def test_legit_extra_route_param_outside_opt_block_passes():
    # A genuine non-runtime route extra outside the opt block is not a leak.
    result = check_gaussian_ts_route(
        "# opt=(ts,calcfc,noeigentest) m062x 6-31G(d) maxcycles=100"
    )

    assert result.verdict == "ok"
