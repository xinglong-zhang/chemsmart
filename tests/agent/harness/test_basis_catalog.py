from __future__ import annotations

from chemsmart.agent.harness.basis_sets import (
    check_basis_intent,
    load_basis_catalog,
    resolve_basis_name,
    search_basis_sets,
)


def test_bse_catalog_has_program_splits_and_known_basis_names():
    catalog = load_basis_catalog()

    assert catalog["metadata"]["source_package"] == "basis_set_exchange"
    assert catalog["metadata"]["basis_set_count"] >= 700
    assert catalog["programs"]["gaussian"]["format"] == "gaussian94"
    assert catalog["programs"]["orca"]["format"] == "orca"
    assert "def2-SVP" in catalog["programs"]["gaussian"]["basis_names"]
    assert "def2-SVP" in catalog["programs"]["orca"]["basis_names"]


def test_resolve_basis_name_accepts_user_spelling_variants():
    result = resolve_basis_name("def2 TZVP", program="gaussian")

    assert result.verdict == "ok"
    assert result.canonical_name == "def2-TZVP"
    assert result.evidence
    assert result.evidence["family"] == "ahlrichs"


def test_check_basis_intent_distinguishes_concrete_from_qualitative():
    concrete = check_basis_intent(
        "Use M06-2X with def2-TZVP and SMD acetonitrile.",
        program="orca",
    )
    qualitative = check_basis_intent(
        "Use a good Karlsruhe triple-zeta basis for this system.",
        program="gaussian",
    )

    assert concrete.verdict == "ok"
    assert concrete.canonical_name == "def2-TZVP"
    assert qualitative.verdict == "ask_user"
    assert qualitative.candidates


def test_unknown_basis_name_is_rejected_with_candidates():
    result = resolve_basis_name("def2-not-real", program="gaussian")

    assert result.verdict == "reject"
    assert result.candidates == ()


def test_search_basis_sets_returns_top_k_only_for_qualitative_basis():
    result = search_basis_sets(
        "Karlsruhe triple zeta diffuse basis",
        program="gaussian",
        limit=4,
    )

    assert result["ok"] is True
    assert result["verdict"] in {"ask_user", "warn", "ok"}
    assert result["result_count"] <= 4
    assert (
        result["token_policy"] == "top_k_only; full catalog is never returned"
    )
    assert any(
        candidate["name"] == "def2-TZVPD" for candidate in result["candidates"]
    )


def test_search_basis_sets_preserves_ri_fit_role():
    result = search_basis_sets(
        "RI fit auxiliary basis for def2 TZVP",
        program="orca",
        limit=5,
    )

    assert result["ok"] is True
    assert result["requested_role"] == "rifit"
    assert result["result_count"] <= 5
    assert result["candidates"][0]["role"] == "rifit"
    assert any(
        candidate["name"] == "def2-TZVP-RIFIT"
        for candidate in result["candidates"]
    )


def test_search_basis_sets_handles_spoken_pople_style_query():
    result = search_basis_sets(
        "Pople split-valence basis with polarization, like six thirty one star",
        program="gaussian",
        limit=6,
    )

    assert result["ok"] is True
    assert result["result_count"] <= 6
    assert any(
        candidate["name"] in {"6-31G*", "6-31G(d)"}
        for candidate in result["candidates"]
    )
