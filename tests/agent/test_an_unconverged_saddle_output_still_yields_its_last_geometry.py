"""Reading a result is linear in its size, converged or not.

REACH-1 po3 (2026-09-06): a saddle search that hit ORCA's iteration
cap terminated normally with 90 coordinate blocks and no rotational
spectrum. The reader's tail decoration indexed an empty list, the
structure reader caught the IndexError and fell to a fallback whose
``else`` sat on the ``if`` inside ``for line in self.contents``, so a
full-file regex scan ran once per line: 37 minutes per read, paid twice
by the executor and charged to the engine wall, then once more at
every wake's bootstrap.
"""

from __future__ import annotations

import time

import pytest

from chemsmart.io.orca.output import ORCAOutput

_BLOCK = """CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O      0.000000    0.000000    {oz:.6f}
  H      0.000000    0.757200   -0.469200
  H      0.000000   -0.757200   -0.469200

----------------------------
CARTESIAN COORDINATES (A.U.)
----------------------------
   NO LB      ZA    FRAG     MASS         X           Y           Z
   0 O     8.0000    0    15.999    0.000000    0.000000    {oz_au:.6f}
   1 H     1.0000    0     1.008    0.000000    1.430900   -0.886660
   2 H     1.0000    0     1.008    0.000000   -1.430900   -0.886660

FINAL SINGLE POINT ENERGY     {energy:.12f}

------------------
CARTESIAN GRADIENT
------------------

   1   O   :    0.000000000    0.000000000    {gz:.9f}
   2   H   :    0.000000000    0.001000000   -0.000500000
   3   H   :    0.000000000   -0.001000000   -0.000500000

"""

_THERMO = """
--------------------------
THERMOCHEMISTRY AT 353.15K
--------------------------

ENTHALPY
--------
Point Group:  C1, Symmetry Number:   1
"""


def _unconverged_saddle_output(path):
    text = (
        "                                 * O   R   C   A *\n"
        "|  1> ! OptTS Freq wb97x-d3bj def2-svp RIJCOSX\n"
        "|  2> * xyz 0 1\n"
        "|  3> O 0.0 0.0 0.1173\n"
        "|  4> H 0.0 0.7572 -0.4692\n"
        "|  5> H 0.0 -0.7572 -0.4692\n"
        "|  6> *\n"
        "****END OF INPUT****\n"
        + _BLOCK.format(oz=0.117300, oz_au=0.221665, energy=-76.4, gz=0.012)
        + _THERMO
        + _BLOCK.format(oz=0.140000, oz_au=0.264562, energy=-76.41, gz=0.004)
        + "       The optimization did not converge but reached the maximum \n"
        "       number of optimization cycles\n"
        "                             ****ORCA TERMINATED NORMALLY****\n"
        "TOTAL RUN TIME: 0 days 0 hours 0 minutes 5 seconds 0 msec\n"
    )
    path.write_text(text)
    return path


@pytest.mark.capability("selector:*")
def test_the_last_geometry_is_read_once_and_in_linear_time(
    tmp_path, monkeypatch
):
    out = _unconverged_saddle_output(tmp_path / "ts_optts.out")
    calls = {"fallback": 0}
    original = ORCAOutput._get_input_structure_in_output

    def counted(self):
        calls["fallback"] += 1
        return original(self)

    monkeypatch.setattr(ORCAOutput, "_get_input_structure_in_output", counted)
    output = ORCAOutput(filename=str(out))
    assert output.normal_termination is True
    assert output.converged is False
    assert output.rotational_constants_in_MHz is None
    assert output.rotational_constants_in_Hz is None

    started = time.perf_counter()
    molecule = output.molecule
    elapsed = time.perf_counter() - started
    assert elapsed < 2.0
    assert (
        calls["fallback"] == 0
    ), "the fallback is for outputs without structures"
    assert [round(float(z), 4) for z in molecule.positions[:, 2]] == [
        0.14,
        -0.4692,
        -0.4692,
    ], "the reached structure is the last block, not the input"
    assert molecule.is_optimized_structure is False
    assert output.molecule is output.molecule, "one read, shared"


def test_the_fallback_runs_once_when_no_structure_can_be_paired(
    tmp_path, monkeypatch
):
    """A single point echoing its input and printing no gradient has no
    structure to pair; the fallback reads the echoed input exactly once."""

    text = (
        "|  1> ! sp wb97x-d3bj def2-svp\n"
        "|  2> * xyz 0 1\n"
        "|  3> O 0.0 0.0 0.1173\n"
        "|  4> H 0.0 0.7572 -0.4692\n"
        "|  5> H 0.0 -0.7572 -0.4692\n"
        "|  6> *\n"
        "****END OF INPUT****\n"
        "FINAL SINGLE POINT ENERGY     -76.400000000000\n"
        "                             ****ORCA TERMINATED NORMALLY****\n"
    )
    out = tmp_path / "sp.out"
    out.write_text(text)
    calls = {"fallback": 0}
    original = ORCAOutput._get_input_structure_in_output

    def counted(self):
        calls["fallback"] += 1
        return original(self)

    monkeypatch.setattr(ORCAOutput, "_get_input_structure_in_output", counted)
    output = ORCAOutput(filename=str(out))
    molecule = output.molecule
    assert calls["fallback"] == 1
    assert len(molecule.chemical_symbols) == 3
    assert output.molecule is molecule
