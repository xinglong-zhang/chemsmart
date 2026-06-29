from __future__ import annotations

import json
from pathlib import Path

from chemsmart.agent.harness.models import HarnessResult


def write_harness_result(session_dir: Path, result: HarnessResult) -> None:
    (session_dir / "harness_result.json").write_text(
        json.dumps(result.to_dict(), indent=2, sort_keys=True),
        encoding="utf-8",
    )
