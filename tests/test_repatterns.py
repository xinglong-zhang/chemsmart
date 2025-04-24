import re

from chemsmart.utils.repattern import (
    gaussian_dias_filename_point_with_fragment1,
    gaussian_dias_filename_point_with_fragment2,
    gaussian_dias_filename_point_without_fragment,
    gaussian_dias_filename_with_reactant,
)


class TestGaussianDiasRegexPatterns:
    def test_point_without_fragment_matches(self):
        pattern = re.compile(gaussian_dias_filename_point_without_fragment)
        assert (
            pattern.match("test_p123_data.log") is not None
        ), "Should match: test_p123_data.log"
        assert pattern.match("test_p123_data.log").groups() == (
            "123",
            "data",
        ), "Correct groups for test_p123_data.log"
        assert (
            pattern.match("test_p123.log") is not None
        ), "Should match: test_p123.log"
        assert pattern.match("test_p123.log").groups() == (
            "123",
            None,
        ), "Correct groups for test_p123.log"

    def test_point_without_fragment_non_matches(self):
        pattern = re.compile(gaussian_dias_filename_point_without_fragment)
        assert (
            pattern.match("test_p123_fdata.log") is None
        ), "Should not match: test_p123_fdata.log"
        assert (
            pattern.match("test_p123_f1_data.log") is None
        ), "Should not match: test_p123_f1_data.log"
        assert (
            pattern.match("test_p123_.log") is None
        ), "Should not match: test_p123_.log"

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
