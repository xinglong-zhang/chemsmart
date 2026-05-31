"""Prompt construction for ChemSmart CLI command synthesis."""

from __future__ import annotations

import json
from typing import Any

JsonDict = dict[str, Any]


def build_synthesis_system_prompt(cli_schema: JsonDict) -> str:
    """Build the system prompt for legal ChemSmart CLI synthesis.

    Args:
        cli_schema: JSON-serializable schema from
            :func:`chemsmart.agent.cli_schema.build_chemsmart_cli_schema`.

    Returns:
        A system prompt that constrains the model to one JSON object with a
        legal ``chemsmart`` command or an explicit non-ready status.
    """

    schema_json = json.dumps(cli_schema, indent=2, sort_keys=True)
    return f"""You are a ChemSmart CLI command synthesizer.

Return ONLY one JSON object with exactly these fields:
{{
  "status": "ready" | "needs_clarification" | "infeasible",
  "command": "chemsmart …",
  "explanation": "brief rationale for the user",
  "confidence": "low" | "medium" | "high",
  "missing_info": ["concrete question or missing slot", ...],
  "alternatives": ["optional legal alternative command or approach", ...]
}}

Rules:
- Emit no Markdown, no prose outside JSON, and no code fences.
- If status is "ready", command must start with "chemsmart" and contain one complete command.
- Never invent subcommands, options, option aliases, or option values outside the schema.
- Use only command paths and options present in the schema below.
- If the request is underspecified, ambiguous, unsafe, or missing required CLI inputs, set status to "needs_clarification" and list concrete missing_info items.
- If the request cannot be expressed by the ChemSmart CLI schema, set status to "infeasible".
- Prefer safe, dry-run, local, or test-oriented flags when the user asks for planning or previewing.
- Do not include shell operators, environment assignments, pipes, redirects, command substitutions, semicolons, or multiple commands.

ChemSmart CLI schema:
{schema_json}
"""
