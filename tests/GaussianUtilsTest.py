from chemsmart.utils.utils import (
    get_list_from_string_range,
    str_indices_to_list,
)


class TestGetListFromStringRange:
    def test_get_list_from_string_range(self):
        s1 = "[1-3,28-31,34-41]"
        s1_list = get_list_from_string_range(string_of_range=s1)
        assert s1_list == [1, 2, 3, 28, 29, 30, 31, 34, 35, 36, 37, 38, 39, 40, 41]

        s2 = "1-3,28-31,34-41"
        s2_list = get_list_from_string_range(string_of_range=s2)
        assert s2_list == [1, 2, 3, 28, 29, 30, 31, 34, 35, 36, 37, 38, 39, 40, 41]

        s3 = "1,3,33,37,42,43,44,45"
        s3_list = get_list_from_string_range(string_of_range=s3)
        assert s3_list == [1, 3, 33, 37, 42, 43, 44, 45]


class TestGetListFromStringZeroIndexed:
    def test_get_list_from_string(self):
        s1 = "1:9"
        s1_list = str_indices_to_list(str_indices=s1)
        assert s1_list == [1, 2, 3, 4, 5, 6, 7, 8]

        s2 = "1,2,4"
        s2_list = str_indices_to_list(str_indices=s2)
        assert s2_list == [1, 2, 4]

        s3 = "1-9"
        s3_list = str_indices_to_list(str_indices=s3)
        assert s3_list == [1, 2, 3, 4, 5, 6, 7, 8]

        s4 = "[1-9]"
        s4_list = str_indices_to_list(str_indices=s4)
        assert s4_list == [1, 2, 3, 4, 5, 6, 7, 8]
