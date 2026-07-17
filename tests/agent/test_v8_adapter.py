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
                    "parent_job": "opt",
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
        "--nstates 10 --root 2 --eqsolv eqsolv"
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


def test_db_record_selector_renders_as_program_option():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.sp",
                "file": "results.db",
                "record_index": 1,
                "structure_index": "2",
                "charge": 0,
                "mult": 1,
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -p test --record-index 1 "
        "--structure-index 2 -f results.db -c 0 -m 1 sp"
    ]


def test_postprocess_promotes_db_selectors_from_settings():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "orca.sp",
                "file": "results.db",
                "charge": 0,
                "mult": 1,
                "settings": {"record_id": "abc123", "structure_index": "1"},
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="orca")

    assert out["valid"], out["errors"]
    assert out["spec"]["jobs"][0]["record_id"] == "abc123"
    assert out["spec"]["jobs"][0]["structure_index"] == "1"
    assert "settings" not in out["spec"]["jobs"][0]
    assert (
        "chemsmart run orca -p orca --record-id abc123 --structure-index 1 "
        "-f results.db -c 0 -m 1 sp"
    ) == out["commands"][0]


def test_gaussian_scan_definition_renders_runtime_scan_flags():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.scan",
                "file": "anisole.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {
                    "scan_definition": "B 1 2 S 12 0.08",
                    "additional_route_parameters": "scf=tight",
                },
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -p test --additional-route-parameters "
        "scf=tight -f anisole.xyz -c 0 -m 1 scan --coordinates "
        "'[[1,2]]' --num-steps 12 --step-size 0.08"
    ]


def test_gaussian_scan_definition_with_constraints_renders_runtime_flags():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.scan",
                "file": "scan.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {
                    "scan_definition": (
                        "B 1 2 S 12 0.08\n" "A 3 4 5 S 8 2.0\n" "D 1 2 3 4 F"
                    ),
                },
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -p test -f scan.xyz -c 0 -m 1 scan "
        "--coordinates '[[1,2],[3,4,5]]' --num-steps '[12,8]' "
        "--step-size '[0.08,2.0]' --constrained-coordinates "
        "'[[1,2,3,4]]'"
    ]


def test_invalid_gaussian_scan_definition_returns_adapter_error():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.scan",
                "file": "scan.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"scan_definition": "B 1 2 F"},
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec), default_project="test")

    assert out["valid"] is False
    assert out["commands"] == []
    assert out["errors"] == [
        "adapter render failed: scan_definition must include at least one "
        "S scan entry"
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


def test_qmmm_renders_under_explicit_parent_with_total_state():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.qmmm",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {
                    "parent_job": "sp",
                    "high_level_atoms": [1],
                    "medium_level_atoms": [2],
                    "low_level_atoms": [2, 3],
                },
                "label": "water_qmmm",
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec))

    assert out["valid"], out["errors"]
    assert out["commands"] == [
        "chemsmart run gaussian -f water.xyz -c 0 -m 1 "
        "-l water_qmmm sp qmmm --high-level-atoms 1 "
        "--medium-level-atoms 2 --low-level-atoms 2,3 "
        "--charge-total 0 --mult-total 1"
    ]


def test_qmmm_renders_explicit_layer_state_and_orca_mm_method():
    specs = [
        {
            "intent": "workflow",
            "jobs": [
                {
                    "id": 1,
                    "kind": "gaussian.qmmm",
                    "file": "enzyme.xyz",
                    "charge": 0,
                    "mult": 1,
                    "settings": {
                        "parent_job": "opt",
                        "high_level_atoms": [1, 2, 3],
                        "medium_level_atoms": [4, 5],
                        "low_level_atoms": [6, 7],
                        "charge_total": 0,
                        "mult_total": 1,
                        "charge_intermediate": 0,
                        "mult_intermediate": 1,
                        "charge_high": 0,
                        "mult_high": 1,
                    },
                }
            ],
        },
        {
            "intent": "workflow",
            "jobs": [
                {
                    "id": 1,
                    "kind": "orca.qmmm",
                    "file": "enzyme.xyz",
                    "charge": -1,
                    "mult": 1,
                    "settings": {
                        "parent_job": "opt",
                        "high_level_atoms": [1, 2, 3],
                        "jobtype": "QMMM",
                        "low_level_method": "AMBER=HardFirst",
                        "charge_total": -1,
                        "mult_total": 1,
                    },
                }
            ],
        },
    ]

    gaussian = v8_adapter.adapt(json.dumps(specs[0]), default_project="demo")
    orca = v8_adapter.adapt(json.dumps(specs[1]), default_project="demo")

    assert gaussian["valid"], gaussian["errors"]
    assert orca["valid"], orca["errors"]
    assert (
        "--charge-intermediate 0 --mult-intermediate 1"
        in gaussian["commands"][0]
    )
    assert (
        "--jobtype QMMM --low-level-method AMBER=HardFirst"
        in orca["commands"][0]
    )


def test_qmmm_without_parent_is_not_rendered_as_top_level_command():
    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "orca.qmmm",
                "file": "water.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"high_level_atoms": [1, 2]},
            }
        ],
    }

    out = v8_adapter.adapt(json.dumps(spec))

    assert out["valid"] is False
    assert out["commands"] == []
    assert "requires settings.parent_job" in out["errors"][0]


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
    fixed, changed = disambiguate(
        "Wiberg bond index analysis", copy.deepcopy(dias)
    )
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
