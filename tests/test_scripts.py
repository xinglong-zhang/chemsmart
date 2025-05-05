import pytest

from chemsmart.scripts.fukui import entry_point as fukui_entry_point


class TestFukuiScripts:
    def test_missing_radical_files(self, gaussian_opt_inputfile):
        """Test input validation: missing both radical files."""
        with pytest.raises(
            ValueError,
            match="At least one of radical cation or radical anion files must be provided",
        ):
            fukui_entry_point.callback(
                neutral_filename=gaussian_opt_inputfile,
                radical_cation_filename=None,
                radical_anion_filename=None,
                mode="mulliken",
            )

    def test_invalid_file_type(self):
        """Test invalid file type handling."""
        with pytest.raises(
            TypeError, match="File neutral.txt is of unknown filetype"
        ):
            fukui_entry_point.callback(
                neutral_filename="neutral.txt",
                radical_cation_filename="cation.log",
                radical_anion_filename=None,
                mode="mulliken",
            )

    def test_invalid_mode(
        self, gaussian_hirshfeld_outfile, gaussian_rc_hirshfeld_outfile
    ):
        """Test invalid charge mode handling."""
        with pytest.raises(
            ValueError,
            match="Unknown mode invalid. Supported modes are: mulliken, nbo, hirshfeld, cm5.",
        ):
            fukui_entry_point.callback(
                neutral_filename=gaussian_hirshfeld_outfile,
                radical_cation_filename=gaussian_rc_hirshfeld_outfile,
                radical_anion_filename=None,
                mode="invalid",
            )

    def test_gaussian_mulliken(
        self,
        gaussian_hirshfeld_outfile,
        gaussian_rc_hirshfeld_outfile,
        capture_log,
    ):
        """Test Gaussian file handling with Mulliken charges."""
        fukui_entry_point.callback(
            neutral_filename=gaussian_hirshfeld_outfile,
            radical_cation_filename=gaussian_rc_hirshfeld_outfile,
            radical_anion_filename=None,
            mode="mulliken",
        )
        print("********************")
        print([record.message for record in capture_log.records])
        print(capture_log.text)
        assert (
            "Using Mulliken Charges for computing Fukui Reactivity Indices."
            in capture_log.text
        )
        assert "Neutral System Charges:" in capture_log.text
        assert "O1     :   -0.360" in capture_log.text
        assert "O2     :   -0.317" in capture_log.text
        assert "C3     :   -0.090" in capture_log.text
        assert "H33    :    0.183" in capture_log.text
        assert "Radical Cationic System Charges:" in capture_log.text
        assert "O1     :    0.020" in capture_log.text
        assert "O2     :   -0.317" in capture_log.text
        assert "C3     :   -0.088" in capture_log.text
        assert "H33    :    0.184" in capture_log.text
        assert "Fukui Minus (f-)" in capture_log.text
        assert "O1        0.380" in capture_log.text
        assert "O2       -0.000" in capture_log.text
        assert "C3        0.003" in capture_log.text
        assert "H33       0.000" in capture_log.text
        assert (
            "Fukui Plus(f+)    Fukui Zero(f0)    Fukui Dual Descriptor(f(2))"
            in capture_log.text
        )
        assert "O1        0.000     0.000     0.000" in capture_log.text
        assert "O2        0.000     0.000     0.000" in capture_log.text
        assert "C3        0.000     0.000     0.000" in capture_log.text
        assert "H33       0.000     0.000     0.000" in capture_log.text

    def test_gaussian_hirshfeld(
        self,
        gaussian_hirshfeld_outfile,
        gaussian_rc_hirshfeld_outfile,
        capture_log,
    ):
        """Test Gaussian file handling with Hirshfeld charges."""
        fukui_entry_point.callback(
            neutral_filename=gaussian_hirshfeld_outfile,
            radical_cation_filename=gaussian_rc_hirshfeld_outfile,
            radical_anion_filename=None,
            mode="hirshfeld",
        )
        assert (
            "Using Hirshfeld Charges for computing Fukui Reactivity Indices"
            in capture_log.text
        )
        assert "Neutral System Charges:" in capture_log.text
        assert "O1     :   -0.222" in capture_log.text
        assert "O2     :   -0.176" in capture_log.text
        assert "C3     :   -0.030" in capture_log.text
        assert "H33    :    0.050" in capture_log.text
        assert "Radical Cationic System Charges:" in capture_log.text
        assert "O1     :    0.100" in capture_log.text
        assert "O2     :   -0.169" in capture_log.text
        assert "C3     :   -0.001" in capture_log.text
        assert "H33    :    0.051" in capture_log.text
        assert "Fukui Minus (f-)" in capture_log.text
        assert "O1        0.322" in capture_log.text
        assert "O2        0.006" in capture_log.text
        assert "C3        0.030" in capture_log.text
        assert "H33       0.000" in capture_log.text
        assert "O1        0.000     0.000     0.000" in capture_log.text
        assert "O2        0.000     0.000     0.000" in capture_log.text
        assert "C3        0.000     0.000     0.000" in capture_log.text
        assert "H33       0.000     0.000     0.000" in capture_log.text

    def test_gaussian_cm5(
        self,
        gaussian_hirshfeld_outfile,
        gaussian_rc_hirshfeld_outfile,
        capture_log,
    ):
        """Test Gaussian file handling with CM5 charges."""
        fukui_entry_point.callback(
            neutral_filename=gaussian_hirshfeld_outfile,
            radical_cation_filename=gaussian_rc_hirshfeld_outfile,
            radical_anion_filename=None,
            mode="cm5",
        )
        assert (
            "Using CM5 Charges for computing Fukui Reactivity Indices"
            in capture_log.text
        )
        assert "Neutral System Charges:" in capture_log.text
        assert "O1     :   -0.310" in capture_log.text
        assert "O2     :   -0.279" in capture_log.text
        assert "C3     :   -0.090" in capture_log.text
        assert "Fukui Minus (f-)" in capture_log.text
        assert "O1        " in capture_log.text
        assert "O2        " in capture_log.text
        assert "C3        " in capture_log.text
        assert "O1        0.000     0.000     0.000" in capture_log.text
        assert "O2        0.000     0.000     0.000" in capture_log.text
        assert "C3        0.000     0.000     0.000" in capture_log.text

    def test_gaussian_nbo(
        self,
        gaussian_hirshfeld_outfile,
        gaussian_rc_hirshfeld_outfile,
        capture_log,
    ):
        """Test Gaussian file handling with NBO charges (skipped)."""
        fukui_entry_point.callback(
            neutral_filename=gaussian_hirshfeld_outfile,
            radical_cation_filename=gaussian_rc_hirshfeld_outfile,
            radical_anion_filename=None,
            mode="nbo",
        )
        assert (
            "Using NBO Charges for computing Fukui Reactivity Indices"
            in capture_log.text
        )
        # Add more asserts here if NBO charges available
