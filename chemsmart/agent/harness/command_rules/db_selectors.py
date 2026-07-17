"""Database selector-flag rules for Gaussian/ORCA job submission."""

from __future__ import annotations

from chemsmart.agent.harness.command_rules.models import (
    CommandContractIssue,
    reject,
)
from chemsmart.agent.harness.command_rules.tokens import (
    flag_option_present,
    flag_option_value,
)

DB_SELECTOR_FLAGS = {
    "record_index": ("--ri", "--record-index"),
    "record_id": ("--rid", "--record-id"),
    "structure_index": ("--si", "--structure-index"),
    "structure_id": ("--sid", "--structure-id"),
    "molecule_id": ("--mid", "--molecule-id"),
}
DB_RECORD_SELECTORS = {"record_index", "record_id"}
DB_TOP_LEVEL_SELECTORS = {"record_index", "record_id", "structure_id"}


def db_selector_issues(
    program: str,
    program_tokens: list[str],
) -> list[CommandContractIssue]:
    present = {
        key: flag_option_value(program_tokens, aliases)
        for key, aliases in DB_SELECTOR_FLAGS.items()
        if flag_option_present(program_tokens, aliases)
    }
    source = flag_option_value(program_tokens, ("-f", "--filename"))
    is_db = bool(source and str(source).lower().endswith(".db"))
    if not present and not is_db:
        return []
    issues: list[CommandContractIssue] = []

    if "molecule_id" in present:
        issues.append(
            reject(
                "cmd.semantic.db_molecule_id_job",
                (
                    "--mid/--molecule-id is not supported for Gaussian/ORCA "
                    "job submission; select a structure with --sid or "
                    "record plus structure index"
                ),
                {
                    "program": program,
                    "filename": source,
                    "selectors": sorted(present),
                },
            )
        )

    if not is_db:
        issues.append(
            reject(
                "cmd.semantic.db_selector_without_db",
                "database selectors require a .db input file",
                {
                    "program": program,
                    "filename": source,
                    "selectors": sorted(present),
                },
            )
        )
        return issues

    top_level = [key for key in DB_TOP_LEVEL_SELECTORS if key in present]
    if len(top_level) != 1:
        issues.append(
            reject(
                "cmd.semantic.db_selector_cardinality",
                (
                    ".db job input must select exactly one of "
                    "--ri/--record-index, --rid/--record-id, or "
                    "--sid/--structure-id"
                ),
                {
                    "program": program,
                    "filename": source,
                    "selectors": sorted(present),
                },
            )
        )

    if "structure_index" in present and not (
        DB_RECORD_SELECTORS & set(present)
    ):
        issues.append(
            reject(
                "cmd.semantic.db_structure_index_requires_record",
                (
                    "--si/--structure-index can only be used together with "
                    "--ri/--record-index or --rid/--record-id"
                ),
                {
                    "program": program,
                    "filename": source,
                    "selectors": sorted(present),
                },
            )
        )
    return issues


__all__ = [
    "DB_RECORD_SELECTORS",
    "DB_SELECTOR_FLAGS",
    "DB_TOP_LEVEL_SELECTORS",
    "db_selector_issues",
]
