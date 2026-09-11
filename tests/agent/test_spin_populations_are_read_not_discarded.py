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

from chemsmart.io.orca.output import ORCAOutput

pytestmark = pytest.mark.capability(
    "selector:mulliken_atomic_spin_populations"
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
