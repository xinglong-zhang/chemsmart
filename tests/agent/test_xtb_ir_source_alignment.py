"""A paired xTB IR table must also own the returned frequency ordering."""

from pathlib import Path
from types import SimpleNamespace

from chemsmart.analysis.result_readers import _xtb_native_evidence_paths
from chemsmart.io.xtb.output import XTBOutput


class _ModeTable:
    def __init__(self, frequencies, intensities, path):
        self.vibrational_frequencies = frequencies
        self.ir_intensities = intensities
        self.filepath = Path(path)


class _OutputWithModeTables(XTBOutput):
    """Minimal output double that drives the public source-selection API."""

    def __init__(self, *, main, vibspectrum):
        self._main = main
        self._vibspectrum = vibspectrum

    @property
    def molecule(self):
        return SimpleNamespace(is_monoatomic=False)

    @property
    def main_out(self):
        return self._main

    @property
    def g98_file(self):
        return None

    @property
    def vibspectrum_file(self):
        return self._vibspectrum


def test_ir_mode_pair_uses_one_native_mode_table_not_equal_length_lists():
    """A matching length alone never establishes normal-mode identity."""

    main = _ModeTable([100.0, 200.0, 300.0], [1.0, 2.0, 3.0], "/tmp/xtb.out")
    vibspectrum = _ModeTable(
        [101.0, 201.0, 301.0], [10.0, 20.0, 30.0], "/tmp/vibspectrum"
    )
    output = _OutputWithModeTables(main=main, vibspectrum=vibspectrum)

    assert output.ir_spectrum_source is vibspectrum
    assert output.vibrational_frequencies == [101.0, 201.0, 301.0]
    assert output.ir_intensities == [10.0, 20.0, 30.0]
    assert _xtb_native_evidence_paths(output, "vibrational_frequencies") == (
        Path("/tmp/vibspectrum"),
    )
    assert _xtb_native_evidence_paths(output, "ir_intensities") == (
        Path("/tmp/vibspectrum"),
    )
