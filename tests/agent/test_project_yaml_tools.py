from __future__ import annotations

import yaml

from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
    ResolvedDecision,
)
from chemsmart.agent.project_yaml import (
    critic_project_yaml,
    extract_project_protocol,
    read_project_yaml,
    render_project_yaml,
    update_project_yaml,
    validate_project_yaml,
    write_project_yaml,
)
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.registry import ToolRegistry

PROFESSOR_CO2_PROMPT = """
I want to set up a project yaml for my project on CO2, you can name it as
co2.yaml using the following reported methods:

The model catalyst was conformationally sampled to locate the most stable
complex. The conformational sampling was carried out using Grimme's CREST
program, which used metadynamics (MTD) with genetic z-matrix crossing (GC)
performed at the GFN2-xTB extended semiempirical tight-binding level of theory
with opt=vtight option. Ten of the lowest energy GFN2-xTB optimized structures
from the CREST search were further optimized using density functional theory
(DFT), implemented in Gaussian16 rev. B.01 software, in the gas phase using the
B3LYP hybrid functional with Grimme's D3 dispersion correction with
Becke-Johnson damping (hereafter denoted B3LYP-D3BJ) and the def2-SVPD
Karlsruhe-family basis set for Br atom and def2-SVP basis set for all other
atoms. Minima and transition structures on the potential energy surface (PES)
were confirmed as such by harmonic frequency analysis, showing respectively
zero and one imaginary frequency.
"""


def test_professor_co2_prompt_renders_valid_gaussian_project_yaml():
    protocol = extract_project_protocol(
        PROFESSOR_CO2_PROMPT,
        project_name="co2.yaml",
        program="gaussian",
    )
    rendered = render_project_yaml(protocol)

    assert protocol["project_name"] == "co2"
    assert protocol["method"]["functional_route"] == (
        "b3lyp empiricaldispersion=gd3bj"
    )
    assert protocol["method"]["heavy_elements"] == ["Br"]
    assert protocol["method"]["heavy_elements_basis"] == "def2svpd"
    assert protocol["method"]["light_elements_basis"] == "def2svp"
    assert (
        "CREST/GFN2-xTB conformer sampling workflow"
        in protocol["unsupported_yaml_features"]
    )

    parsed = yaml.safe_load(rendered["yaml_text"])
    expected_block = {
        "functional": "b3lyp empiricaldispersion=gd3bj",
        "basis": "gen",
        "freq": True,
        "heavy_elements": ["Br"],
        "heavy_elements_basis": "def2svpd",
        "light_elements_basis": "def2svp",
    }
    expected_solv_block = dict(expected_block)
    expected_solv_block["freq"] = False
    assert parsed == {"gas": expected_block, "solv": expected_solv_block}

    validation = validate_project_yaml(
        rendered["yaml_text"],
        program="gaussian",
        project_name="co2",
    )
    assert validation["verdict"] == "ok"
    assert validation["runtime_summary"]["opt"]["functional"] == (
        "b3lyp empiricaldispersion=gd3bj"
    )
    assert validation["runtime_summary"]["opt"]["basis"] == "gen"

    critic = critic_project_yaml(
        rendered["yaml_text"],
        protocol=protocol,
        program="gaussian",
        project_name="co2",
    )
    assert critic["verdict"] == "warn"
    assert any(
        issue["rule_id"] == "protocol.unsupported_yaml_feature"
        for issue in critic["issues"]
    )


def test_project_yaml_handles_cited_basis_and_method_only_render_input():
    protocol = extract_project_protocol(
        "Gaussian16 gas phase B3LYP-D3BJ with def2-SVPD[12,13] "
        "Karlsruhe-family basis set for Br atomand def2-SVP[12,14] "
        "for all other atoms. Frequency analysis confirmed minima.",
        project_name="co2.yaml",
        program="gaussian",
    )

    rendered = render_project_yaml(
        protocol["method"],
        project_name="co2",
        program="gaussian",
    )
    parsed = yaml.safe_load(rendered["yaml_text"])

    assert protocol["method"]["heavy_elements"] == ["Br"]
    assert protocol["method"]["heavy_elements_basis"] == "def2svpd"
    assert parsed["gas"]["functional"] == "b3lyp empiricaldispersion=gd3bj"
    assert parsed["gas"]["basis"] == "gen"
    assert parsed["gas"]["heavy_elements"] == ["Br"]


def test_project_yaml_extracts_solvent_and_rejects_missing_method():
    protocol = extract_project_protocol(
        "Use Gaussian optimization with M06-2X-D3BJ/def2-TZVP and confirm "
        "the minimum by frequency analysis. Use SMD(acetonitrile) for the "
        "solvated single-point stage.",
        project_name="nitrile",
        program="gaussian",
    )
    rendered = render_project_yaml(protocol)
    parsed = yaml.safe_load(rendered["yaml_text"])

    assert parsed["gas"]["functional"] == "m062x empiricaldispersion=gd3bj"
    assert parsed["gas"]["basis"] == "def2tzvp"
    assert "solvent_model" not in parsed["gas"]
    assert parsed["solv"]["freq"] is False
    assert parsed["solv"]["solvent_model"] == "smd"
    assert parsed["solv"]["solvent_id"] == "acetonitrile"

    invalid = validate_project_yaml(
        "gas:\n  functional: null\n  basis: def2svp\nsolv:\n"
        "  functional: null\n  basis: def2svp\n",
        program="gaussian",
    )
    assert invalid["verdict"] == "reject"
    assert any(
        issue["rule_id"] == "yaml.method_missing_functional"
        for issue in invalid["issues"]
    )


def test_project_yaml_preserves_a_distinct_gaussian_td_method():
    protocol = extract_project_protocol(
        "Optimize in the gas phase using PBE0/def2-SVP. For TD-DFT, use "
        "CAM-B3LYP/def2-SVP to calculate the absorption spectrum.",
        project_name="photo",
        program="gaussian",
    )
    rendered = render_project_yaml(protocol)
    parsed = yaml.safe_load(rendered["yaml_text"])

    assert protocol["method"]["functional"] == "pbe0"
    assert protocol["td"]["functional"] == "camb3lyp"
    assert parsed["gas"]["functional"] == "pbe0"
    assert parsed["td"]["functional"] == "camb3lyp"
    assert parsed["td"]["basis"] == "def2svp"

    validation = validate_project_yaml(
        rendered["yaml_text"], program="gaussian", project_name="photo"
    )
    assert validation["verdict"] == "ok"
    assert validation["runtime_summary"]["td"]["functional"] == "camb3lyp"

    critic = critic_project_yaml(
        rendered["yaml_text"], protocol=protocol, program="gaussian"
    )
    assert critic["verdict"] == "ok"


def test_single_basis_protocol_renders_and_validates_without_mixed_basis_leak():
    # A plain single-basis protocol must NOT emit heavy/light basis sections;
    # doing so under a non-gen basis trips the chemsmart mixed-basis guard.
    protocol = extract_project_protocol(
        "All structures were optimized in water using the SMD implicit "
        "solvation model at the B3LYP-D3BJ/def2-SVP level of theory in "
        "Gaussian 16. Frequency analysis confirmed all minima.",
        project_name="h2o",
        program="gaussian",
    )
    assert protocol["method"]["heavy_elements"] == []
    assert protocol["method"]["light_elements_basis"] is None

    rendered = render_project_yaml(protocol, project_name="h2o")
    parsed = yaml.safe_load(rendered["yaml_text"])
    assert parsed["gas"]["basis"] == "def2svp"
    assert "light_elements_basis" not in parsed["gas"]
    assert "heavy_elements" not in parsed["gas"]

    validation = validate_project_yaml(
        rendered["yaml_text"], program="gaussian", project_name="h2o"
    )
    assert validation["verdict"] == "ok"
    assert validation["runtime_summary"]["opt"]["basis"] == "def2svp"

    critic = critic_project_yaml(
        rendered["yaml_text"],
        protocol=protocol,
        program="gaussian",
        project_name="h2o",
    )
    assert critic["verdict"] == "ok"


def test_project_yaml_canonicalizes_model_supplied_mixed_basis_method_dict():
    rendered = render_project_yaml(
        {
            "functional": "B3LYP-D3BJ",
            "basis": "def2-SVP",
            "freq": True,
            "heavy_elements": ["I"],
            "heavy_elements_basis": "def2-SVPD",
        },
        project_name="iodobenzene_mixed",
        program="gaussian",
    )
    parsed = yaml.safe_load(rendered["yaml_text"])

    assert parsed["gas"]["functional"] == "b3lyp empiricaldispersion=gd3bj"
    assert parsed["gas"]["basis"] == "gen"
    assert parsed["gas"]["heavy_elements"] == ["I"]
    assert parsed["gas"]["heavy_elements_basis"] == "def2svpd"
    assert parsed["gas"]["light_elements_basis"] == "def2svp"


def test_project_yaml_extracts_spoken_light_atom_basis_phrase():
    protocol = extract_project_protocol(
        "I want to create a project yaml named co2.yaml for a CO2 reduction "
        "project. Use Gaussian for DFT refinements. The method should be "
        "B3LYP-D3BJ with def2-SVP for light atoms and def2-SVPD for Br. "
        "Use gas phase optimization followed by harmonic frequency analysis "
        "to confirm minima and transition states. The conformer search was "
        "done externally with CREST at GFN2-xTB, so include that as protocol "
        "context but do not make it the Gaussian calculation method.",
        project_name="co2.yaml",
        program="gaussian",
    )
    rendered = render_project_yaml(protocol)
    parsed = yaml.safe_load(rendered["yaml_text"])

    assert protocol["method"]["basis"] == "gen"
    assert protocol["method"]["heavy_elements_basis"] == "def2svpd"
    assert protocol["method"]["light_elements_basis"] == "def2svp"
    assert parsed["gas"]["basis"] == "gen"
    assert parsed["gas"]["heavy_elements"] == ["Br"]
    assert parsed["gas"]["heavy_elements_basis"] == "def2svpd"
    assert parsed["gas"]["light_elements_basis"] == "def2svp"

    critic = critic_project_yaml(
        rendered["yaml_text"],
        protocol=protocol,
        program="gaussian",
        project_name="co2",
    )
    assert critic["verdict"] == "warn"
    assert all(issue["severity"] != "reject" for issue in critic["issues"])
    assert (
        validate_project_yaml(
            rendered["yaml_text"],
            program="gaussian",
            project_name="iodobenzene_mixed",
        )["verdict"]
        == "ok"
    )


def test_validate_accepts_render_result_dict_chained_from_model():
    # The tool-loop model often passes the whole render_project_yaml result as
    # yaml_text; the harness must unwrap the yaml_text field instead of erroring.
    protocol = extract_project_protocol(
        "Optimize in water with SMD at B3LYP-D3BJ/def2-SVP; freq confirms "
        "minima.",
        project_name="h2o",
        program="gaussian",
    )
    rendered = render_project_yaml(protocol, project_name="h2o")

    from_string = validate_project_yaml(
        rendered["yaml_text"], program="gaussian", project_name="h2o"
    )
    from_dict = validate_project_yaml(
        rendered, program="gaussian", project_name="h2o"
    )
    assert from_string["verdict"] == from_dict["verdict"] == "ok"

    critic = critic_project_yaml(
        rendered, protocol=protocol, program="gaussian", project_name="h2o"
    )
    assert critic["verdict"] == "ok"


def test_validate_dedups_identical_candidate(monkeypatch):
    # Re-validating an unchanged candidate must not repeat the runtime loader
    # (dedup guard against build-mode re-validation loops).
    import chemsmart.agent.project_yaml as pj

    pj._VALIDATION_CACHE.clear()
    yaml_text = (
        "gas:\n  functional: b3lyp\n  basis: def2svp\n  freq: true\n"
        "solv:\n  functional: b3lyp\n  basis: def2svp\n  freq: false\n"
    )

    calls = {"n": 0}
    real_loader = pj._load_project_yaml_via_runtime

    def counting_loader(**kwargs):
        calls["n"] += 1
        return real_loader(**kwargs)

    monkeypatch.setattr(pj, "_load_project_yaml_via_runtime", counting_loader)

    first = pj.validate_project_yaml(yaml_text, program="gaussian")
    second = pj.validate_project_yaml(yaml_text, program="gaussian")

    assert first["verdict"] == second["verdict"] == "ok"
    assert calls["n"] == 1  # runtime loader ran once, second hit the cache
    assert "revalidation_skipped" not in first
    assert second["revalidation_skipped"] is True


def test_project_yaml_validator_rejects_mixed_basis_without_gen():
    yaml_text = """
gas:
  functional: b3lyp empiricaldispersion=gd3bj
  basis: def2svp
  heavy_elements: [Br]
  heavy_elements_basis: def2svpd
  light_elements_basis: def2svp
"""

    result = validate_project_yaml(yaml_text, program="gaussian")

    assert result["verdict"] == "reject"
    assert any(
        issue["rule_id"] == "yaml.gaussian.mixed_basis_without_gen"
        for issue in result["issues"]
    )


def test_project_yaml_write_uses_workspace_config_dir(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    yaml_text = """
gas:
  functional: b3lyp empiricaldispersion=gd3bj
  basis: def2svp
  freq: true
solv:
  functional: b3lyp empiricaldispersion=gd3bj
  basis: def2svp
  freq: false
"""

    result = write_project_yaml("co2.yaml", yaml_text, program="gaussian")

    target = tmp_path / ".chemsmart" / "gaussian" / "co2.yaml"
    assert result["ok"] is True
    assert result["written_path"] == str(target)
    assert target.read_text(encoding="utf-8").endswith("\n")

    from chemsmart.settings.gaussian import GaussianProjectSettings

    settings = GaussianProjectSettings.from_project("co2")
    assert settings.opt_settings().functional == (
        "b3lyp empiricaldispersion=gd3bj"
    )


def test_project_yaml_read_and_update_use_workspace_config_dir(
    monkeypatch,
    tmp_path,
):
    monkeypatch.chdir(tmp_path)
    yaml_text = """
gas:
  functional: b3lyp
  basis: def2svp
  freq: true
solv:
  functional: b3lyp
  basis: def2svp
  freq: false
"""
    written = write_project_yaml("demo", yaml_text, program="gaussian")
    assert written["ok"] is True

    loaded = read_project_yaml()
    assert loaded["ok"] is True
    assert loaded["project_name"] == "demo"
    assert loaded["program"] == "gaussian"
    assert loaded["parsed"]["gas"]["functional"] == "b3lyp"

    updated = update_project_yaml(
        {"gas.functional": "m062x", "solv.functional": "m062x"},
        project_name="demo",
        program="gaussian",
    )
    assert updated["ok"] is True
    assert "-  functional: b3lyp" in updated["diff"]
    assert "+  functional: m062x" in updated["diff"]

    loaded_after = read_project_yaml("demo", "gaussian")
    assert loaded_after["parsed"]["gas"]["functional"] == "m062x"


def test_project_yaml_tools_are_registered_and_write_requires_approval():
    registry = ToolRegistry.default()
    names = {tool.name for tool in registry.list_tools()}

    assert {
        "extract_project_protocol",
        "render_project_yaml",
        "validate_project_yaml",
        "critic_project_yaml",
        "write_project_yaml",
        "read_project_yaml",
        "update_project_yaml",
        "synthesize_command",
        "repair_command",
        "execute_chemsmart_command",
    }.issubset(names)

    request = ToolRequest(
        request_id="req",
        provider="test",
        provider_call_id="call",
        name="write_project_yaml",
        arguments_json="{}",
        arguments={},
        raw={},
    )
    policy = PermissionPolicy(mode=PermissionMode.PERMISSION)
    resolved = policy.resolve(request)
    assert resolved.decision == ResolvedDecision.NEEDS_USER
    policy.record("write_project_yaml", ApprovalDecision.ALLOW_SESSION)
    assert "write_project_yaml" not in policy.session_allow


def test_plain_d3_dispersion_is_preserved_not_dropped():
    # "b3lyp-d3" (zero-damping D3, not BJ) must render the gd3 route keyword;
    # it was previously silently dropped because only d3bj was handled.
    for functional in ("b3lyp-d3", "B3LYP empiricaldispersion=gd3"):
        rendered = render_project_yaml(
            {
                "functional": functional,
                "basis": "genecp",
                "freq": True,
                "heavy_elements": ["Br"],
                "heavy_elements_basis": "SDD",
                "light_elements_basis": "6-31G**",
            },
            project_name="ch3br",
            program="gaussian",
        )
        parsed = yaml.safe_load(rendered["yaml_text"])
        assert parsed["gas"]["functional"] == (
            "b3lyp empiricaldispersion=gd3"
        ), functional

    # D3BJ still renders gd3bj (no regression).
    rendered = render_project_yaml(
        {"functional": "b3lyp-d3bj", "basis": "def2svp", "freq": True},
        project_name="x",
        program="gaussian",
    )
    parsed = yaml.safe_load(rendered["yaml_text"])
    assert parsed["gas"]["functional"] == "b3lyp empiricaldispersion=gd3bj"

    # ORCA plain D3 maps to the ORCA D3 keyword.
    rendered = render_project_yaml(
        {"functional": "pbe0", "dispersion": "d3", "basis": "def2tzvp"},
        project_name="cat",
        program="orca",
    )
    parsed = yaml.safe_load(rendered["yaml_text"])
    assert parsed["gas"]["dispersion"] == "D3"
