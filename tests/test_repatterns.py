import re

from chemsmart.utils.repattern import (
    gaussian_dias_filename_point_with_fragment1,
    gaussian_dias_filename_point_with_fragment2,
    gaussian_dias_filename_point_without_fragment,
    gaussian_dias_filename_with_reactant,
    gaussian_freq_keywords_pattern,
    gaussian_opt_keywords_pattern,
    multiple_spaces_pattern,
)


def test_gaussian_opt_keywords_pattern():
    """Test pattern for matching optimization keywords in Gaussian route strings."""
    pattern = re.compile(gaussian_opt_keywords_pattern, re.IGNORECASE)

    # Test basic opt keyword
    assert pattern.search("# b3lyp opt freq") is not None, "Should match: opt"
    # Test opt with parameters
    assert (
        pattern.search("# b3lyp OPT=tight freq") is not None
    ), "Should match: OPT=tight"
    assert (
        pattern.search("# b3lyp opt = tight freq") is not None
    ), "Should match: opt = tight"
    # Test opt with parentheses
    assert (
        pattern.search("# b3lyp opt=(calcfc,tight) freq") is not None
    ), "Should match: opt=(calcfc,tight)"
    assert (
        pattern.search("# b3lyp opt = (calcfc,tight) freq") is not None
    ), "Should match: opt = (calcfc,tight)"


def test_gaussian_freq_keywords_pattern():
    """Test pattern for matching frequency keywords in Gaussian route strings."""
    pattern = re.compile(gaussian_freq_keywords_pattern, re.IGNORECASE)

    # Test basic freq keyword
    assert pattern.search("# b3lyp opt freq") is not None, "Should match: freq"
    # Test freq with parameters
    assert (
        pattern.search("# b3lyp opt FREQ=numer") is not None
    ), "Should match: freq=numer"
    assert (
        pattern.search("# b3lyp opt freq = analytical") is not None
    ), "Should match: freq = analytical"
    # Test that it doesn't match partial words
    assert (
        pattern.search("# b3lyp opt frequency") is not None
    ), "Should not match: frequency (partial word)"


def test_multiple_spaces_pattern():
    """Test pattern for matching multiple consecutive spaces."""
    pattern = re.compile(multiple_spaces_pattern)

    # Test multiple spaces
    assert (
        pattern.search("word1  word2") is not None
    ), "Should match: double space"
    assert (
        pattern.search("word1   word2") is not None
    ), "Should match: triple space"
    assert (
        pattern.search("word1    word2") is not None
    ), "Should match: four spaces"
    # Test tabs and mixed whitespace
    assert (
        pattern.search("word1\t\tword2") is not None
    ), "Should match: double tab"
    assert (
        pattern.search("word1 \t word2") is not None
    ), "Should match: mixed space/tab"
    assert (
        pattern.search("word1\n\nword2") is not None
    ), "Should match: double newline"
    # Test single space (should still match as it's whitespace)
    assert (
        pattern.search("word1 word2") is not None
    ), "Should match: single space"
    # Test no spaces
    assert pattern.search("word1word2") is None, "Should not match: no spaces"
    # Test whitespace at beginning and end
    assert pattern.search("  word") is not None, "Should match: leading spaces"
    assert (
        pattern.search("word  ") is not None
    ), "Should match: trailing spaces"


def test_gaussian_route_string_cleaning_patterns():
    """Test comprehensive route string cleaning with all three patterns."""
    test_route = "# b3lyp/6-31g(d)  opt=(calcfc,tight,maxstep=5)   freq=numer   scrf=(smd)  scf=qc"

    # Test removing opt keywords
    opt_pattern = re.compile(gaussian_opt_keywords_pattern, re.IGNORECASE)
    route_no_opt = opt_pattern.sub("", test_route)
    assert "opt" not in route_no_opt.lower(), "Should remove opt keywords"

    # Test removing freq keywords
    freq_pattern = re.compile(gaussian_freq_keywords_pattern, re.IGNORECASE)
    route_no_freq = freq_pattern.sub("", test_route)
    assert "freq" not in route_no_freq.lower(), "Should remove freq keywords"

    # Test cleaning multiple spaces
    spaces_pattern = re.compile(multiple_spaces_pattern)
    cleaned_route = spaces_pattern.sub(" ", test_route)
    assert "  " not in cleaned_route, "Should clean multiple spaces"

    # Test full cleaning pipeline
    full_clean = test_route
    full_clean = opt_pattern.sub("", full_clean)
    full_clean = freq_pattern.sub("", full_clean)
    full_clean = spaces_pattern.sub(" ", full_clean)
    full_clean = full_clean.strip()

    expected = "# b3lyp/6-31g(d) scrf=(smd) scf=qc"
    assert (
        full_clean == expected
    ), f"Full cleaning should produce: '{expected}', got: '{full_clean}'"


class TestGaussianDiasRegexPatterns:
    def test_point_without_fragment_matches(self):
        pattern = re.compile(gaussian_dias_filename_point_without_fragment)
        assert (
            pattern.match("test_p123_data.log") is not None
        ), "Should match: test_p123_data.log"
        assert pattern.match("test_p123_data.log").groups() == (
            None,
            None,
            "123",
            "data",
        ), "Correct groups for test_p123_data.log"
        assert (
            pattern.match("test_p123.log") is not None
        ), "Should match: test_p123.log"
        assert pattern.match("test_p123.log").groups() == (
            None,
            None,
            "123",
            None,
        ), "Correct groups for test_p123.log"

    def test_point_without_fragment_non_matches(self):
        pattern = re.compile(gaussian_dias_filename_point_without_fragment)
        assert (
            pattern.match("test_p123_f1_data.log") is None
        ), "Should not match: test_p123_f1_data.log"
        assert (
            pattern.match("test_p123_.log") is None
        ), "Should not match: test_p123_.log"
        assert (
            pattern.match("NiRRBenz_c2_scan_p16_ts_dias_p1_f1.log") is None
        ), "Should not match: NiRRBenz_c2_scan_p16_ts_dias_p1_f1.log"
        assert (
            pattern.match("NiRRBenz_c2_scan_p16_ts_dias_p1_f2.log") is None
        ), "Should not match: NiRRBenz_c2_scan_p16_ts_dias_p1_f2.log"

    def test_point_with_fragment1_matches(self):
        pattern = re.compile(gaussian_dias_filename_point_with_fragment1)
        assert (
            pattern.match("test_p123_f1_data.log") is not None
        ), "Should match: test_p123_f1_data.log"
        assert pattern.match("test_p123_f1_data.log").groups() == (
            "123",
            "f1",
            "data",
        ), "Correct groups for test_p123_f1_data.log"
        assert (
            pattern.match("test_p123_f1.log") is not None
        ), "Should match: test_p123_f1.log"
        assert pattern.match("test_p123_f1.log").groups() == (
            "123",
            "f1",
            None,
        ), "Correct groups for test_p123_f1.log"

    def test_point_with_fragment1_non_matches(self):
        pattern = re.compile(gaussian_dias_filename_point_with_fragment1)
        assert (
            pattern.match("test_p123.log") is None
        ), "Should not match: test_p123.log"
        assert (
            pattern.match("test_p123_.log") is None
        ), "Should not match: test_p123_.log"
        assert (
            pattern.match("test_p123_f2_data.log") is None
        ), "Should not match: test_p123_f2_data.log"
        assert (
            pattern.match("test_p123_data.log") is None
        ), "Should not match: test_p123_data.log"

    def test_point_with_fragment2_matches(self):
        pattern = re.compile(gaussian_dias_filename_point_with_fragment2)
        assert (
            pattern.match("test_p123_f2_data.log") is not None
        ), "Should match: test_p123_f2_data.log"
        assert pattern.match("test_p123_f2_data.log").groups() == (
            "123",
            "f2",
            "data",
        ), "Correct groups for test_p123_f2_data.log"
        assert (
            pattern.match("test_p123_f2.log") is not None
        ), "Should match: test_p123_f2.log"
        assert pattern.match("test_p123_f2.log").groups() == (
            "123",
            "f2",
            None,
        ), "Correct groups for test_p123_f2.log"

    def test_point_with_fragment2_non_matches(self):
        pattern = re.compile(gaussian_dias_filename_point_with_fragment2)
        assert (
            pattern.match("test_p123.log") is None
        ), "Should not match: test_p123.log"
        assert (
            pattern.match("test_p123_f1_data.log") is None
        ), "Should not match: test_p123_f1_data.log"
        assert (
            pattern.match("test_p123_data.log") is None
        ), "Should not match: test_p123_data.log"

    def test_with_reactant_matches(self):
        pattern = re.compile(gaussian_dias_filename_with_reactant)
        assert (
            pattern.match("test_r1_data.log") is not None
        ), "Should match: test_r1_data.log"
        assert pattern.match("test_r1_data.log").groups() == (
            "1",
            "data",
        ), "Correct groups for test_r1_data.log"
        assert (
            pattern.match("test_r2.log") is not None
        ), "Should match: test_r2.log"
        assert pattern.match("test_r2.log").groups() == (
            "2",
            None,
        ), "Correct groups for test_r2.log"

    def test_with_reactant_non_matches(self):
        pattern = re.compile(gaussian_dias_filename_with_reactant)
        assert (
            pattern.match("test_r3_data.log") is None
        ), "Should not match: test_r3_data.log"
        assert (
            pattern.match("test_p123.log") is None
        ), "Should not match: test_p123.log"
