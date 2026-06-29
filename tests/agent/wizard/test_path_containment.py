from __future__ import annotations

import pytest

from chemsmart.agent.wizard.cache import cache_path
from chemsmart.agent.wizard.verify import verify_server_yaml


def test_cache_path_rejects_path_traversal():
    with pytest.raises(ValueError):
        cache_path("../foo")


def test_verify_server_yaml_reports_invalid_server_name():
    result = verify_server_yaml("../../tmp/x")

    assert result.errors
    assert "Invalid server name" in result.errors[0]
