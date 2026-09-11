"""The label-valued extremum operations are on the stem.

A spin-state session declared the ground state's multiplicity as a
required output and was refused for want of a producer while
coordinate_at_minimum sat behind the scan guide, which its task could
not open (NOVEL-2 ino1, 2026-09-04). Which of several values is lowest
is every comparison question.
"""

from __future__ import annotations

import pytest

from .test_a_guide_opens_when_something_asks import _host, _operations

pytestmark = pytest.mark.capability("operation:coordinate_at_minimum")


def test_the_extremum_operations_are_on_the_bare_stem(tmp_path):
    host = _host(tmp_path)
    operations = _operations(host.surface)
    assert "coordinate_at_minimum" in operations
    assert "coordinate_at_maximum" in operations
    assert "scan" not in host.active_guides
