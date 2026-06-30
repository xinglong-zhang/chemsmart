"""Tests for compact-SPEC postprocess and adapter rendering."""

from __future__ import annotations

import copy
import json

from chemsmart.agent import synthesis as S
from chemsmart.agent import v8_adapter
from chemsmart.agent.kind_disambiguator import disambiguate


def test_postprocess_strips_invalid_and_drops_auto_ts_route():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.nci",
                "file": "a.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"freq": True},
            },
            {
                "id": 2,
                "kind": "gaussian.ts",
                "file": "b.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"ts": "calcfc/noeigentest", "freq": True},
            },
            {
                "id": 3,
                "kind": "gaussian.ts",
                "file": "c.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {
                    "additional_opt_options_in_route": [
                        "ts",
                        "calcfc",
                        "noeigentest",
                        "maxstep=15",
                    ]
                },
            },
        ],
    }

    out = v8_adapter.postprocess(spec)

    assert "settings" not in out["jobs"][0]
    assert "additional_opt_options_in_route" not in out["jobs"][1]["settings"]
    assert out["jobs"][1]["settings"]["freq"] is True
    assert out["jobs"][2]["settings"] == {
        "additional_opt_options_in_route": ["maxstep=15"]
    }


def test_atom_index_settings_render_as_runtime_parseable_ranges():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.dias",
                "file": "dimer.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"fragment_indices": [[1, 2, 3], [4, 5, 6]]},
            },
            {
                "id": 2,
                "kind": "gaussian.qmmm",
                "file": "qm.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {
                    "high_level_atoms": [1, 2, 3],
                    "low_level_atoms": [4, 5, 6],
                },
            },
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec))

    assert out["valid"], out["errors"]
    assert "--fragment-indices 1,2,3" in out["commands"][0]
    assert "--high-level-atoms 1,2,3" in out["commands"][1]
    assert "--low-level-atoms 4,5,6" in out["commands"][1]


def test_freq_renders_via_route_param_not_fake_flag():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "orca.opt",
                "file": "m.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"freq": True},
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec))

    assert out["valid"], out["errors"]
    assert "--additional-route-parameters freq" in out["commands"][0]
    assert out["commands"][0].endswith(" opt")
    assert "opt_freq" not in out["commands"][0]
    assert "--freq" not in out["commands"][0]
    assert "--jobtype" not in out["commands"][0]


def test_gaussian_sp_renders_real_sp_subcommand():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.sp",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -p test -f water.xyz -c 0 -m 1 sp"
    ]
    assert "singlepoint" not in out["commands"][0]


def test_orca_sp_renders_real_sp_subcommand():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "orca.sp",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run orca -p test -f water.xyz -c 0 -m 1 sp"
    ]
    assert "singlepoint" not in out["commands"][0]


def test_gaussian_freq_renders_real_opt_route_freq():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.freq",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -p test --additional-route-parameters freq "
        "-f water.xyz -c 0 -m 1 opt"
    ]
    assert "--jobtype" not in out["commands"][0]


def test_orca_program_options_render_before_job_subcommand():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "orca.opt",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
                "label": "water_orca_opt",
                "settings": {"additional_route_parameters": "GridX7"},
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec))

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run orca --additional-route-parameters GridX7 -f water.xyz "
        "-c 0 -m 1 -l water_orca_opt opt"
    ]


def test_gaussian_tddft_options_render_after_td_subcommand():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.tddft",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"nstates": 10, "root": 2, "eqsolv": True},
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -p test -f water.xyz -c 0 -m 1 td "
        "--nstates 10 --root 2 --eqsolv true"
    ]


def test_gaussian_tddft_drops_invalid_numeric_states_choice():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.tddft",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"nstates": 5, "states": 2, "root": 2},
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -p test -f water.xyz -c 0 -m 1 td "
        "--nstates 5 --root 2"
    ]


def test_default_project_renders_as_runtime_owned_program_option():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.opt",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -p test -f water.xyz -c 0 -m 1 opt"
    ]


def test_gaussian_dias_fragment_indices_render_after_subcommand():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.dias",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"fragment_indices": [1, 2]},
                "label": "water_dias",
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec))

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -f water.xyz -c 0 -m 1 "
        "-l water_dias dias --fragment-indices 1,2"
    ]


def test_qmmm_atom_indices_render_after_subcommand():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.qmmm",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"high_level_atoms": [1], "low_level_atoms": [2, 3]},
                "label": "water_qmmm",
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec))

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -f water.xyz -c 0 -m 1 "
        "-l water_qmmm qmmm --high-level-atoms 1 --low-level-atoms 2,3"
    ]


def test_orca_neb_options_render_after_subcommand():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "orca.neb",
                "file": "reactant.xyz",
                "product_file": "product.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"nimages": 6, "joboption": "NEB-TS"},
                "label": "neb_job",
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec))

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run orca -f reactant.xyz -c 0 -m 1 -l neb_job "
        "neb --nimages 6 --joboption NEB-TS -e product.xyz"
    ]


def test_adapt_chain_renders_ordered_commands():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.opt",
                "file": "m.xyz",
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
    }

    out = v8_adapter.adapt(json.dumps(spec))

    assert len(out["commands"]) == 2
    assert out["valid"], out["errors"]
    assert out["commands"][1].startswith("chemsmart sub -s chemnode1 orca")
    assert "<gaussian.opt-output>" in out["commands"][1]


def test_synthesis_routes_compact_spec():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.opt",
                "file": "a.xyz",
                "charge": 0,
                "mult": 1,
            }
        ],
    }

    assert S._is_v8_spec(spec) is True
    result = S._normalize_v8_spec(spec)

    assert result["status"] == "ready"
    assert result["command"].startswith("chemsmart ")


def test_synthesis_compact_decline_is_infeasible():
    result = S._normalize_v8_spec(
        {"intent": "decline", "message": "missing fragment indices"}
    )

    assert result["status"] == "infeasible"
    assert result["command"] == ""
    assert "fragment" in result["explanation"]


def test_kind_disambiguator_fixes_high_confidence_confusions():
    dias = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.dias",
                "file": "m.xyz",
                "charge": 0,
                "mult": 1,
            }
        ],
    }
    fixed, changed = disambiguate("Wiberg bond index analysis", copy.deepcopy(dias))
    assert changed is True
    assert fixed["jobs"][0]["kind"] == "gaussian.wbi"

    opt = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.opt",
                "file": "m.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"freeze_atoms": [3, 4]},
            }
        ],
    }
    fixed, changed = disambiguate(
        "constrained opt freezing the bond between 3 and 4",
        copy.deepcopy(opt),
    )
    assert changed is True
    assert fixed["jobs"][0]["kind"] == "gaussian.modred"
    assert fixed["jobs"][0]["settings"]["modred"] == [[3, 4]]

    orca_opt = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "orca.opt",
                "file": "m.xyz",
                "charge": 0,
                "mult": 1,
            }
        ],
    }
    fixed, changed = disambiguate(
        "Run ORCA OptTS locally for m.xyz, neutral singlet.",
        copy.deepcopy(orca_opt),
    )
    assert changed is True
    assert fixed["jobs"][0]["kind"] == "orca.ts"

    plain_opt, changed = disambiguate(
        "Run ORCA opt for m.xyz, neutral singlet.",
        copy.deepcopy(orca_opt),
    )
    assert changed is False
    assert plain_opt["jobs"][0]["kind"] == "orca.opt"


def test_gaussian_modred_renders_python_list_literal_coordinates():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.modred",
                "file": "examples/h2o.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"modred": [[1, 2]]},
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -p test -f examples/h2o.xyz -c 0 -m 1 "
        "modred --coordinates '[[1,2]]'"
    ]
