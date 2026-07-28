import pytest

from chemsmart.utils.cluster import normalize_scheduler_job_label


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("pka_scale_pka_batch", "pka_scale_pka_batch"),
        ("array_pka_scale_pka_batch", "pka_scale_pka_batch"),
        ("pka_scale_pka_batch_array", "pka_scale_pka_batch"),
        ("water_sp", "water_sp"),
    ],
)
def test_normalize_scheduler_job_label(raw, expected):
    assert normalize_scheduler_job_label(raw) == expected
