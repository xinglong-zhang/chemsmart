"""Request-scoped prompt packs and hard synthesis budget."""

from __future__ import annotations

import pytest

from chemsmart.agent.cli_schema import build_chemsmart_cli_schema
from chemsmart.agent.prompts.synthesis import (
    SYNTHESIS_PROMPT_MAX_CHARS,
    build_synthesis_system_prompt,
)
from chemsmart.agent.schema_prune import prune_schema_for_request


@pytest.fixture(scope="module")
def schema():
    return build_chemsmart_cli_schema()


def _prompt(schema, request: str) -> str:
    return build_synthesis_system_prompt(
        prune_schema_for_request(schema, request)
    )


@pytest.mark.parametrize(
    "query",
    [
        "optimize h2o.xyz locally with Gaussian, charge 0 multiplicity 1",
        "submit an ORCA TS for ts.xyz to chemnode1",
        "Gaussian TDDFT 5 singlet states root 1 for pyridine.xyz",
        "scan Gaussian bond 1-2 and freeze 1-3",
        "ORCA NEB with reactant.xyz product.xyz and 8 images",
        "help me with my molecule",
    ],
)
def test_request_scoped_prompts_fit_hard_budget(schema, query: str) -> None:
    assert len(_prompt(schema, query)) <= SYNTHESIS_PROMPT_MAX_CHARS


def test_base_policy_keeps_routing_and_option_collisions(schema) -> None:
    prompt = _prompt(schema, "submit Gaussian opt of h2o.xyz on chemnode1")
    assert "chemsmart sub" in prompt
    assert "chemsmart run" in prompt
    assert "-p/--project" in prompt
    assert "-P/--pubchem" in prompt
    assert "--num-steps" in prompt and "--nstates" in prompt
    assert "debug" in prompt and "verbose" in prompt


def test_td_pack_is_present_only_for_td_request(schema) -> None:
    td = _prompt(schema, "Gaussian TDDFT 5 singlet states root 1")
    opt = _prompt(schema, "Gaussian optimization of h2o.xyz")
    for token in ("--states", "--nstates", "--root", "--eqsolv"):
        assert token in td
    assert "Do not invent `--singlets`" in td
    assert "ORCA has no td" in td
    assert "--eqsolv" not in opt


def test_scan_pack_preserves_scan_modred_disambiguation(schema) -> None:
    prompt = _prompt(schema, "Gaussian relaxed scan bond 1-2 freeze 1-3")
    assert "Gaussian relaxed scan" in prompt
    assert "gaussian|orca modred" in prompt
    assert "Vary/range -> scan" in prompt
    assert "--constrained-coordinates" in prompt


def test_nci_pack_distinguishes_raw_and_gaussian(schema) -> None:
    prompt = _prompt(schema, "run NCIPLOT on dimer.wfn locally")
    assert "nciplot -f FILE" in prompt
    assert "Gaussian workflow" in prompt


def test_xtb_pack_teaches_no_project_and_option_placement(schema) -> None:
    xtb = _prompt(schema, "GFN2-xTB optimization of water.xyz")
    assert "needs no project YAML" in xtb
    assert "hess" in xtb
    assert "BOTH `-sm` and `-si`" in xtb
    assert "before the leaf" in xtb
    # Explicit Gaussian prunes xtb away, and the pack must vanish with it.
    gaussian_only = _prompt(schema, "Gaussian optimization of h2o.xyz")
    assert "needs no project YAML" not in gaussian_only


def test_base_policy_orders_program_options_before_subcommand(
    schema,
) -> None:
    prompt = _prompt(schema, "Gaussian optimization of h2o.xyz")
    assert "between the program name and its job subcommand" in prompt


def test_xtb_prompts_fit_hard_budget(schema) -> None:
    for query in (
        "GFN-FF pre-optimization of ligand.xyz before a DFT TS search",
        "xTB hessian of water.xyz with ALPB water",
        "semiempirical single point of cluster.xyz",
    ):
        assert len(_prompt(schema, query)) <= SYNTHESIS_PROMPT_MAX_CHARS


def test_compact_signature_drops_verbose_schema_fields(schema) -> None:
    prompt = _prompt(schema, "Gaussian optimization of h2o.xyz")
    assert '"commands"' in prompt
    assert '"description"' not in prompt
    assert '"default"' not in prompt
    assert '"help"' not in prompt
    assert '"nargs"' not in prompt


def test_orca_qmmm_signature_fits_budget_without_repeated_arg_keys(
    schema,
) -> None:
    prompt = _prompt(
        schema,
        "Submit an ORCA QM/MM optimization for enzyme.xyz to hpc1 with "
        "QM atoms 1-12 and AMBER=HardFirst.",
    )

    assert len(prompt) <= SYNTHESIS_PROMPT_MAX_CHARS
    assert "--high-level-atoms" in prompt
    assert '"flag"' not in prompt
    assert '"type":"str"' not in prompt
