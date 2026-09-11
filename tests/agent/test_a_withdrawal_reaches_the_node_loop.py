"""The human's withdrawal reaches the executor, not only the cycle.

The charter promises cancellation at any node boundary, and the
executor checks for it before launching each node. The driver's own
execute hook dropped the stop file, so that check could never fire:
`GoalDriver._stopped()` read the file between cycles and nothing read
it between nodes. Observed twice on real windows -- REACH-1 po3 ran
2.5 hours of nodes after its STOP, and OPEN-1 po3 launched a
transition-state search fourteen minutes after one (2026-09-07).
"""

from __future__ import annotations

import inspect

import pytest

from chemsmart.agent import driver as driver_module


@pytest.mark.capability("rule:execution.cancelled.human")
def test_the_default_execute_hook_carries_the_stop_file():
    parameters = inspect.signature(driver_module._default_execute).parameters
    assert "stop_file" in parameters
    source = inspect.getsource(driver_module._default_execute)
    assert "stop_file=stop_file" in source


@pytest.mark.capability("rule:execution.cancelled.human")
def test_the_driver_offers_the_withdrawal_to_a_hook_that_takes_it(tmp_path):
    seen: dict[str, object] = {}

    def hook(*, approval_file, workspace, run_directory, stop_file=None):
        seen["stop_file"] = stop_file
        return None

    def older_hook(*, approval_file, workspace, run_directory):
        seen["called"] = True
        return None

    assert driver_module._execute_hook_takes_stop_file(hook)
    assert not driver_module._execute_hook_takes_stop_file(older_hook)
    assert driver_module._execute_hook_takes_stop_file(lambda **kw: None)


@pytest.mark.capability("rule:execution.cancelled.human")
def test_the_executor_stops_before_the_next_node():
    from chemsmart.agent import executor as executor_module

    source = inspect.getsource(executor_module.ApprovedWorkflowExecutor)
    assert "self.should_stop is not None and self.should_stop()" in source
    assert "execution.cancelled.human" in source
