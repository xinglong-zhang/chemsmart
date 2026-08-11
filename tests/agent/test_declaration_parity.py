"""What ChemSmart supports and what it declares must not drift apart.

`settings/capabilities.py` is the registry a caller reads to learn what a
project YAML may set.  ORCA's relativistic, reference-determinant, frozen-core
and element-basis settings existed in `jobs/orca/settings.py` and emitted
correct `%basis`/`%method`/`%scf` blocks, yet were absent from that registry --
so an agent reading the declaration could not ask for any of them, and a live
run reported them as capability gaps that were really declaration gaps.

These tests fail when a settings field becomes settable without being declared.
"""

import pytest

from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

#: Fields every program's settings object carries that a *project* never owns.
#: Identity comes from the bound geometry, the job type from the CLI leaf, and
#: the rest are rendering scratch rather than scientific configuration.
_NEVER_PROJECT_OWNED = frozenset(
    {
        "charge",
        "multiplicity",
        "jobtype",
        "title",
        "input_string",
        "route_to_be_written",
        "gen_genecp_file",
        "gbw",
        "modred",
        "invert_constraints",
        "chk",
        "heavy_elements",
        "light_elements_basis",
    }
)


def _settings_fields(settings_class):
    return {
        name
        for name in vars(settings_class.default())
        if not name.startswith("_")
    }


@pytest.mark.parametrize(
    "program,import_path",
    [
        ("orca", "chemsmart.jobs.orca.settings:ORCAJobSettings"),
        ("gaussian", "chemsmart.jobs.gaussian.settings:GaussianJobSettings"),
    ],
)
def test_settable_scientific_fields_are_declared(program, import_path):
    """A field a project YAML can set must appear in the registry."""

    import importlib

    module_name, class_name = import_path.split(":")
    settings_class = getattr(importlib.import_module(module_name), class_name)
    declared = set(PROGRAM_CAPABILITIES[program].project_owned_parameters)
    settable = _settings_fields(settings_class) - _NEVER_PROJECT_OWNED
    undeclared = sorted(settable - declared)
    assert not undeclared, (
        f"{program} settings expose {undeclared} but "
        "settings/capabilities.py does not declare them; a caller reading the "
        "registry cannot request them. Declare them or add them to "
        "_NEVER_PROJECT_OWNED with a reason."
    )


@pytest.mark.xfail(
    reason=(
        "Pre-existing drift, not introduced here and not yet verified: the "
        "shared union declares names such as scf_tol, aux_basis and "
        "solvent_options that no attribute of the settings object carries. "
        "They may be phantom promises, or accepted by from_dict under other "
        "attribute names. Recorded so the reverse direction is visible; "
        "resolving it needs a pass over each program's from_dict."
    ),
    strict=False,
)
@pytest.mark.parametrize("program", sorted(PROGRAM_CAPABILITIES))
def test_a_declared_parameter_is_actually_settable(program):
    """The reverse drift: declaring a parameter no settings object accepts."""

    import importlib

    modules = {
        "orca": "chemsmart.jobs.orca.settings:ORCAJobSettings",
        "gaussian": "chemsmart.jobs.gaussian.settings:GaussianJobSettings",
    }
    if program not in modules:
        pytest.skip(f"{program} settings are not dict-shaped")
    module_name, class_name = modules[program].split(":")
    settings_class = getattr(importlib.import_module(module_name), class_name)
    declared = set(PROGRAM_CAPABILITIES[program].project_owned_parameters)
    phantom = sorted(declared - _settings_fields(settings_class))
    assert not phantom, (
        f"{program} declares {phantom} which its settings object has no field "
        "for; the declaration promises a control that does not exist."
    )


def test_the_relativistic_correlated_settings_are_reachable():
    """The exact gaps a live run reported are closed."""

    declared = set(PROGRAM_CAPABILITIES["orca"].project_owned_parameters)
    assert {
        "relativistic",
        "reference",
        "frozen_core",
        "frozen_core_electrons",
        "ri_approximation",
        "heavy_elements_basis",
        "mdci_cutoff",
    } <= declared


@pytest.mark.parametrize("program", sorted(PROGRAM_CAPABILITIES))
def test_a_project_taking_program_declares_its_section_names(program):
    """Without this, a wrong section shape is only caught deep in a loader."""

    capability = PROGRAM_CAPABILITIES[program]
    if not capability.supports_project_configuration:
        return
    assert (
        capability.project_section_names
    ), f"{program} accepts a project YAML but declares no section names"
