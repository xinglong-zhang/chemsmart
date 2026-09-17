"""ORCA's per-atom spin populations are declared, read by position, and
checked against 2S.

A session asked how the spin splits between nickel and two sulfurs; ORCA
printed it in the block the reader already parsed for charges, and the
reader kept the charge column and dropped the spin column in the same
loop (NOVEL-3 ino3, 2026-09-05). The proxy the session was left with
inverted the chemistry.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.io.orca.output import ORCAOutput

pytestmark = pytest.mark.capability(
    "selector:orca:sp:mulliken_atomic_spin_populations"
)

_QUARTET = Path("tests/data/ORCATests/outputs/fe3_quartet.out")
_CLOSED = Path("tests/data/ORCATests/outputs/CO2.out")


def test_the_reader_keeps_the_spin_column():
    output = ORCAOutput(str(_QUARTET))
    mulliken = output.mulliken_atomic_spin_populations
    loewdin = output.loewdin_atomic_spin_populations
    assert mulliken["Fe1"] == pytest.approx(2.950529, abs=1e-6)
    assert sum(mulliken.values()) == pytest.approx(3.0, abs=0.01)
    assert sum(loewdin.values()) == pytest.approx(3.0, abs=0.01)
    assert ORCAOutput(str(_CLOSED)).mulliken_atomic_spin_populations is None


def test_the_selector_is_declared_and_refused_on_a_closed_shell():
    from chemsmart.analysis.result_readers import reader_for

    reader = reader_for("orca")
    for jobtype in ("sp", "opt", "ts"):
        declared = reader.selectors_for_jobtype(jobtype)
        assert "mulliken_atomic_spin_populations" in declared
        assert "loewdin_atomic_spin_populations" in declared
    served = reader.read(
        reader.open_output(_QUARTET), "mulliken_atomic_spin_populations"
    )
    values, unit = served
    assert unit == "1"
    assert values[0] == pytest.approx(2.950529, abs=1e-6)
    assert sum(values) == pytest.approx(3.0, abs=0.01)
    with pytest.raises(Exception, match="printed no spin populations"):
        reader.read(
            reader.open_output(_CLOSED), "mulliken_atomic_spin_populations"
        )


@pytest.mark.capability("tool:inspect_run")
@pytest.mark.capability("tool:extract_result_quantities")
def test_orca_spin_populations_reach_the_host_evidence_plane_by_scheme(
    tmp_path,
):
    """The Agent sees two named partitions, never one generic spin density."""

    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="orca-spin-populations"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    host.artifacts["fe3-quartet"] = TrustedArtifactRefV1(
        artifact_id="fe3-quartet",
        kind="orca_output",
        sha256=file_sha256(_QUARTET),
        size_bytes=_QUARTET.stat().st_size,
        path=str(_QUARTET.resolve()),
        cli_value=str(_QUARTET.resolve()),
    )

    inspected = host._inspect_run(
        "t1", {"program": "orca", "artifact_id": "fe3-quartet"}
    )
    metadata = inspected["atom_resolved_metadata"]
    assert metadata["mulliken_atomic_spin_populations"] == {
        "semantic_quantity": "atomic_spin_population",
        "population_scheme": "Mulliken",
        "atom_order": "zero-based molecular atom order",
    }
    assert metadata["loewdin_atomic_spin_populations"] == {
        "semantic_quantity": "atomic_spin_population",
        "population_scheme": "Loewdin",
        "atom_order": "zero-based molecular atom order",
    }

    receipt = host._extract_result_quantities(
        "t2",
        {
            "program": "orca",
            "artifact_id": "fe3-quartet",
            "selectors": [
                {
                    "quantity_id": "mulliken_spin",
                    "selector": "mulliken_atomic_spin_populations",
                },
                {
                    "quantity_id": "loewdin_spin",
                    "selector": "loewdin_atomic_spin_populations",
                },
            ],
        },
    )
    values = {item.quantity_id: item.value for item in receipt.quantities}
    assert values["mulliken_spin"][0] == pytest.approx(2.950529)
    assert sum(values["mulliken_spin"]) == pytest.approx(3.0, abs=0.01)
    assert sum(values["loewdin_spin"]) == pytest.approx(3.0, abs=0.01)
    assert dict(receipt.selector_bindings) == {
        "mulliken_spin": "mulliken_atomic_spin_populations",
        "loewdin_spin": "loewdin_atomic_spin_populations",
    }
