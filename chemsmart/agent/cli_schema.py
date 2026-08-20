"""Deterministic, read-only projection of the live Click command tree."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import click

from chemsmart.agent._contracts import ContractError, canonical_sha256


@dataclass(frozen=True)
class ClickOptionV1:
    parameter_name: str
    primary_flag: str
    flags: tuple[str, ...]
    secondary_flags: tuple[str, ...]
    required: bool
    multiple: bool
    nargs: int
    is_flag: bool
    type_name: str


@dataclass(frozen=True)
class ClickCommandV1:
    path: tuple[str, ...]
    options: tuple[ClickOptionV1, ...]
    child_names: tuple[str, ...]

    def option(self, parameter_name: str) -> ClickOptionV1 | None:
        normalized = str(parameter_name or "").strip().replace("-", "_")
        return next(
            (
                item
                for item in self.options
                if item.parameter_name == normalized
            ),
            None,
        )


@dataclass(frozen=True)
class LiveClickSchemaV1:
    schema_version: str
    commands: tuple[ClickCommandV1, ...]
    schema_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.live-click-schema.v1":
            raise ContractError("unsupported live Click schema version")
        paths = tuple(item.path for item in self.commands)
        if paths != tuple(sorted(set(paths))):
            raise ContractError(
                "Click command paths must be sorted and unique"
            )
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "commands": self.commands,
            }
        )
        if self.schema_sha256 != expected:
            raise ContractError("live Click schema digest mismatch")

    def command(self, path: Iterable[str]) -> ClickCommandV1 | None:
        normalized = tuple(str(item) for item in path)
        return next(
            (item for item in self.commands if item.path == normalized), None
        )

    def has_path(self, path: Iterable[str]) -> bool:
        return self.command(path) is not None


def build_live_click_schema(
    root: click.Command | None = None,
) -> LiveClickSchemaV1:
    """Walk Click command objects without invoking callbacks or engines."""

    if root is None:
        from chemsmart.cli.main import entry_point

        root = entry_point
    records: list[ClickCommandV1] = []

    def walk(command: click.Command, path: tuple[str, ...]) -> None:
        options = []
        for parameter in command.params:
            if not isinstance(parameter, click.Option):
                continue
            long_flags = tuple(
                sorted(
                    flag for flag in parameter.opts if flag.startswith("--")
                )
            )
            secondary = tuple(sorted(parameter.secondary_opts))
            preferred = (
                max(long_flags, key=lambda item: (len(item), item))
                if long_flags
                else sorted(parameter.opts)[0]
            )
            options.append(
                ClickOptionV1(
                    parameter_name=str(parameter.name),
                    primary_flag=preferred,
                    flags=tuple(sorted(parameter.opts)),
                    secondary_flags=secondary,
                    required=bool(parameter.required),
                    multiple=bool(parameter.multiple),
                    nargs=int(parameter.nargs),
                    is_flag=bool(parameter.is_flag),
                    type_name=type(parameter.type).__name__,
                )
            )
        child_map = getattr(command, "commands", None)
        children = child_map if isinstance(child_map, dict) else {}
        records.append(
            ClickCommandV1(
                path=path,
                options=tuple(
                    sorted(options, key=lambda item: item.parameter_name)
                ),
                child_names=tuple(sorted(str(name) for name in children)),
            )
        )
        for name, child in sorted(
            children.items(), key=lambda item: str(item[0])
        ):
            walk(child, path + (str(name),))

    walk(root, ())
    ordered = tuple(sorted(records, key=lambda item: item.path))
    body = {
        "schema_version": "chemsmart.live-click-schema.v1",
        "commands": ordered,
    }
    return LiveClickSchemaV1(**body, schema_sha256=canonical_sha256(body))


__all__ = [
    "ClickCommandV1",
    "ClickOptionV1",
    "LiveClickSchemaV1",
    "build_live_click_schema",
]
