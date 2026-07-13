from chemsmart.agent.harness.sub_intent import build_sub_intent_assertions


def _failed_ids(rows):
    return {row["id"] for row in rows if row["status"] == "fail"}


def test_valid_sub_command_preserves_file_and_research_intent():
    command = (
        "chemsmart sub -s mock-pbs --fake gaussian -p mock "
        "-f examples/h2o.xyz -c 0 -m 1 -r freq opt"
    )
    rows = build_sub_intent_assertions(
        command,
        {
            "program": "gaussian",
            "job": "opt",
            "server": "mock-pbs",
            "filename": "examples/h2o.xyz",
            "charge": 0,
            "multiplicity": 1,
            "route_contains": "freq",
        },
    )
    assert _failed_ids(rows) == set()


def test_successful_submission_with_pubchem_substitution_is_rejected():
    command = (
        "chemsmart sub -s mock-pbs --fake gaussian -p mock -P water "
        "-c 0 -m 1 -r freq opt"
    )
    rows = build_sub_intent_assertions(
        command,
        {
            "program": "gaussian",
            "job": "opt",
            "server": "mock-pbs",
            "filename": "examples/h2o.xyz",
            "charge": 0,
            "multiplicity": 1,
        },
    )
    assert "sub.intent_filename" in _failed_ids(rows)


def test_submission_rejects_project_and_kind_drift():
    command = (
        "chemsmart sub -s mock-pbs --fake gaussian -p changed "
        "-f examples/h2o.xyz -c 0 -m 1 sp"
    )
    rows = build_sub_intent_assertions(
        command,
        {
            "program": "gaussian",
            "kind": "gaussian.opt",
            "project": "requested",
            "server": "mock-pbs",
            "filename": "examples/h2o.xyz",
        },
    )

    assert {"sub.intent_kind", "sub.intent_project"} <= _failed_ids(rows)


def test_neb_assertions_cover_endpoint_and_image_count():
    command = (
        "chemsmart sub -s mock-pbs --fake orca -p mock "
        "-f examples/reactant.xyz -c 0 -m 1 neb "
        "-e examples/product.xyz --nimages 8"
    )
    rows = build_sub_intent_assertions(
        command,
        {
            "program": "orca",
            "job": "neb",
            "server": "mock-pbs",
            "filename": "examples/reactant.xyz",
            "charge": 0,
            "multiplicity": 1,
            "ending_xyzfile": "examples/product.xyz",
            "nimages": 8,
        },
    )
    assert _failed_ids(rows) == set()


def test_qmmm_layer_state_uses_qmmm_charge_aliases():
    command = (
        "chemsmart sub -s mock-pbs --fake orca -p mock "
        "-f examples/enzyme.xyz -c 0 -m 1 opt qmmm "
        "-j QMMM -ct 0 -mt 1 -ha 1-4 -lm AMBER=HardFirst"
    )
    rows = build_sub_intent_assertions(
        command,
        {
            "program": "orca",
            "job": "opt",
            "server": "mock-pbs",
            "filename": "examples/enzyme.xyz",
            "charge": 0,
            "multiplicity": 1,
            "required_tokens": ["qmmm"],
        },
    )
    assert _failed_ids(rows) == set()
