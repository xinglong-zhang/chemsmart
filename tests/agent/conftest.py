from __future__ import annotations

from types import MappingProxyType, SimpleNamespace

import click
import pytest

from chemsmart.agent.capabilities import (
    build_command_compiled_preview_overlay,
    build_program_component_conformance_receipt,
    load_program_capabilities,
)
from chemsmart.agent.cli_schema import build_live_click_schema


def _execution_group(name: str, *, include_test: bool):
    decorators = [
        click.option("--fake/--no-fake", default=False),
        click.option("--scratch/--no-scratch", default=True),
    ]
    if include_test:
        decorators.append(click.option("--test/--no-test", default=False))

    def callback(**_kwargs):
        return None

    for decorator in decorators:
        callback = decorator(callback)
    return click.group(name=name)(callback)


def _program_group(name: str, *, project_option: bool = True):
    options = [
        click.option("--filename", type=str),
        click.option("--charge", type=int),
        click.option("--multiplicity", type=int),
        click.option("--gpu/--no-gpu", default=None),
    ]
    if project_option:
        options.append(click.option("--project", type=str))

    def callback(**_kwargs):
        return None

    for decorator in options:
        callback = decorator(callback)
    group = click.group(name=name)(callback)
    group.add_command(click.Command(name="sp", callback=lambda: None))
    return group


@pytest.fixture
def fake_click_root():
    root = click.Group(name="chemsmart")
    for target, include_test in (("run", False), ("sub", True)):
        group = _execution_group(target, include_test=include_test)
        group.add_command(_program_group("demo"))
        group.add_command(_program_group("optional"))
        group.add_command(_program_group("pyscf"))
        root.add_command(group)
    return root


@pytest.fixture
def fake_click_schema(fake_click_root):
    return build_live_click_schema(fake_click_root)


@pytest.fixture
def fake_capability_registry():
    mapping = MappingProxyType(
        {
            "demo": SimpleNamespace(
                program="demo",
                requires_project_configuration=True,
                supports_project_configuration=True,
                jobtypes=("sp",),
                project_owned_parameters=("basis", "functional"),
                engines=("cpu",),
            ),
            "optional": SimpleNamespace(
                program="optional",
                requires_project_configuration=False,
                supports_project_configuration=True,
                jobtypes=("sp",),
                project_owned_parameters=("basis",),
                engines=("cpu",),
            ),
            "pyscf": SimpleNamespace(
                program="pyscf",
                requires_project_configuration=True,
                supports_project_configuration=True,
                jobtypes=("sp",),
                project_owned_parameters=("basis", "functional"),
                engines=("cpu", "gpu"),
            ),
        }
    )
    module = SimpleNamespace(PROGRAM_CAPABILITIES=mapping)
    return load_program_capabilities(module)


@pytest.fixture
def fake_support_overlay(fake_capability_registry, fake_click_schema):
    receipts = []
    for ordinal, capability in enumerate(
        fake_capability_registry.programs, start=1
    ):
        base = ordinal * 5
        receipts.append(
            build_program_component_conformance_receipt(
                program=capability.program,
                registry_sha256=fake_capability_registry.registry_sha256,
                live_cli_schema_sha256=fake_click_schema.schema_sha256,
                fixture_bundle_sha256=f"{base:064x}",
                covered_jobtypes=capability.jobtypes,
                covered_engines=capability.engines,
                compiler_receipt_sha256=f"{base + 1:064x}",
                preview_receipt_sha256=f"{base + 2:064x}",
                preflight_receipt_sha256=f"{base + 3:064x}",
                verifier_receipt_sha256=f"{base + 4:064x}",
                compiler_status="passed",
                preview_status="passed",
                preflight_status="passed",
                verifier_status="passed",
            )
        )
    return build_command_compiled_preview_overlay(
        fake_capability_registry,
        conformance_receipts=tuple(receipts),
        live_schema=fake_click_schema,
    )
