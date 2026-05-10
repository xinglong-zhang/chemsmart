from chemsmart.agent.wizard.normalize import (
    FALLBACK_NUM_HOURS,
    WALLTIME_SANITY_CEILING_HOURS,
    choose_queue,
    normalize_resources,
    normalize_walltime,
)
from chemsmart.agent.wizard.parsers import QueueFacts


def _queue(name, **kwargs):
    defaults = dict(
        default=False,
        max_walltime_hours=None,
        default_walltime_hours=None,
        default_mem_gb=None,
        default_cores=None,
        gpus_per_node=None,
        enabled=True,
        started=True,
    )
    defaults.update(kwargs)
    return QueueFacts(name=name, **defaults)


def test_choose_queue_prefers_declared_default():
    assert (
        choose_queue([_queue("long"), _queue("debug", default=True)])
        == "debug"
    )


def test_choose_queue_falls_back_to_shortest_enabled_started_non_gpu():
    assert (
        choose_queue([_queue("gpu", gpus_per_node=1), _queue("std")]) == "std"
    )


def test_choose_queue_returns_none_without_non_gpu_candidate():
    assert (
        choose_queue([_queue("gpu", gpus_per_node=1, default=False)]) is None
    )


def test_normalize_walltime_caps_max_to_sanity_ceiling():
    assert (
        normalize_walltime(_queue("long", max_walltime_hours=240))
        == WALLTIME_SANITY_CEILING_HOURS
    )


def test_normalize_walltime_prefers_queue_default_without_cap():
    assert (
        normalize_walltime(
            _queue("long", default_walltime_hours=48, max_walltime_hours=72)
        )
        == 48
    )


def test_normalize_walltime_uses_queue_max_when_default_missing():
    assert normalize_walltime(_queue("long", max_walltime_hours=72)) == 72


def test_normalize_walltime_falls_back_when_queue_lacks_limits():
    assert normalize_walltime(_queue("long")) == FALLBACK_NUM_HOURS


def test_normalize_walltime_caps_large_queue_max_to_sanity_ceiling():
    assert (
        normalize_walltime(_queue("long", max_walltime_hours=720))
        == WALLTIME_SANITY_CEILING_HOURS
    )


def test_normalize_resources_defaults_gpus_to_zero():
    assert normalize_resources(
        _queue("std", default_mem_gb=64, default_cores=16)
    ) == {
        "mem_gb": 64,
        "cores": 16,
        "gpus": 0,
        "threads": 16,
    }
