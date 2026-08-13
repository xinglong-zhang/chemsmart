"""Gaussian owns the native spelling of generic orbital-basis aliases."""

from io import StringIO
from types import SimpleNamespace

import pytest

from chemsmart.agent.program_verifiers import _settings_match
from chemsmart.agent.projects import (
    project_document,
    project_section_application_observation,
)
from chemsmart.io.gaussian.route import GaussianRoute
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianLinkJobSettings,
    GaussianQMMMJobSettings,
    gaussian_native_basis_token,
)
from chemsmart.jobs.gaussian.writer import GaussianInputWriter


@pytest.mark.parametrize(
    ("generic", "native"),
    (
        ("def2-SVP", "def2svp"),
        ("def2-SVPD", "def2svpd"),
        ("def2-TZVP", "def2tzvp"),
        ("def2-TZVPPD", "def2tzvppd"),
        ("def2-QZVP", "def2qzvp"),
        ("def2-QZVPPD", "def2qzvppd"),
    ),
)
def test_known_generic_def2_aliases_have_native_tokens(generic, native):
    assert gaussian_native_basis_token(generic) == native
    assert gaussian_native_basis_token(generic.lower()) == native


@pytest.mark.parametrize("basis", ("aug-cc-pVDZ", "6-31G(d)", "lanl2dz"))
def test_unrelated_punctuation_and_unknown_names_remain_explicit(basis):
    assert gaussian_native_basis_token(basis) == basis


@pytest.mark.parametrize(
    ("generic", "native"),
    (
        ("def2-SVP", "def2svp"),
        ("def2-TZVP", "def2tzvp"),
        ("def2-QZVP", "def2qzvp"),
        ("aug-cc-pVDZ", "aug-cc-pVDZ"),
        ("6-31G(d)", "6-31G(d)"),
    ),
)
def test_writer_parser_and_preview_share_one_basis_meaning(generic, native):
    settings = GaussianJobSettings(
        functional="B3LYP",
        basis=generic,
        charge=0,
        multiplicity=1,
        jobtype="sp",
        freq=False,
    )
    job = SimpleNamespace(
        settings=settings,
        jobrunner=SimpleNamespace(num_cores=1, mem_gb=1),
        molecule=None,
    )
    route_buffer = StringIO()

    GaussianInputWriter(job)._write_route_section(route_buffer)
    route = route_buffer.getvalue().splitlines()[0]
    parsed_basis = GaussianRoute(route).basis

    assert settings.basis == generic
    assert parsed_basis.casefold() == native.casefold()
    assert not _settings_match(settings, {"basis": generic})


@pytest.mark.parametrize(
    ("generic", "native"),
    (("def2-SVP", "def2svp"), ("aug-cc-pVDZ", "aug-cc-pVDZ")),
)
def test_genecp_fallback_uses_native_rules_without_stripping_punctuation(
    generic, native
):
    settings = GaussianJobSettings(
        functional="B3LYP",
        basis="genecp",
        heavy_elements=("I",),
        heavy_elements_basis="def2-TZVP",
        light_elements_basis=generic,
        charge=0,
        multiplicity=1,
        jobtype="sp",
        freq=False,
    )
    job = SimpleNamespace(
        settings=settings,
        jobrunner=SimpleNamespace(num_cores=1, mem_gb=1),
        molecule=SimpleNamespace(chemical_symbols=("O", "H", "H")),
    )
    route_buffer = StringIO()

    GaussianInputWriter(job)._write_route_section(route_buffer)
    route = route_buffer.getvalue().splitlines()[0]

    assert native in route
    assert "augccpvdz" not in route.casefold()


def test_qmmm_typed_bases_use_the_same_native_vocabulary():
    settings = GaussianQMMMJobSettings(
        parent_jobtype="sp",
        high_level_functional="B3LYP",
        high_level_basis="def2-TZVP",
        low_level_force_field="UFF",
        charge_total=0,
        mult_total=1,
    )

    assert settings.high_level_basis == "def2-TZVP"
    assert "oniom(B3LYP/def2tzvp:UFF)" in settings.route_string


def test_preview_accepts_generic_alias_against_native_parser_value():
    parsed = GaussianJobSettings(basis="def2tzvp")

    assert not _settings_match(parsed, {"basis": "def2-TZVP"})
    assert _settings_match(parsed, {"basis": "cc-pVTZ"})


def test_loader_application_does_not_misreport_native_spelling_as_override():
    document = project_document(
        program="gaussian",
        sections={"gas": {"functional": "B3LYP", "basis": "def2-TZVP"}},
    )

    observation = project_section_application_observation(
        document,
        jobtype="opt",
        applied_settings={"functional": "B3LYP", "basis": "def2tzvp"},
    )

    assert observation["status"] == "effective_settings_present"
    assert "overridden_settings" not in observation


@pytest.mark.parametrize(
    ("settings_basis", "link_route"),
    (
        ("def2-TZVP", "# opt B3LYP/def2-TZVP"),
        ("def2-TZVP", "# opt b3lyp/def2tzvp"),
        ("def2tzvp", "# opt B3LYP/def2-TZVP"),
    ),
)
def test_custom_link_route_recognizes_generic_and_native_basis_aliases(
    settings_basis, link_route,
):
    settings = GaussianLinkJobSettings(
        functional="B3LYP",
        basis=settings_basis,
        charge=0,
        multiplicity=1,
        jobtype="opt",
        link_route=link_route,
    )

    route = settings.link_route_string

    assert route.casefold().count("b3lyp") == 1
    assert route.casefold().count("def2") == 1
    assert "def2-TZVP" not in route
    assert "def2tzvp" in route.casefold()
    assert "geom=check" in route
    assert "guess=read" in route


def test_custom_link_route_materializes_missing_basis_in_native_spelling():
    settings = GaussianLinkJobSettings(
        functional="B3LYP",
        basis="def2-TZVP",
        charge=0,
        multiplicity=1,
        jobtype="opt",
        link_route="# opt",
    )

    route = settings.link_route_string

    assert "B3LYP def2tzvp" in route
    assert "def2-TZVP" not in route


def test_custom_link_route_materializes_ab_initio_method_without_none():
    settings = GaussianLinkJobSettings(
        functional=None,
        ab_initio="MP2",
        basis="def2-TZVP",
        charge=0,
        multiplicity=1,
        jobtype="opt",
        link_route="# opt",
    )

    route = settings.link_route_string

    assert "MP2 def2tzvp" in route
    assert "None" not in route


def test_custom_link_route_materializes_semiempirical_method_without_basis():
    settings = GaussianLinkJobSettings(
        functional=None,
        semiempirical="PM6",
        basis=None,
        charge=0,
        multiplicity=1,
        jobtype="opt",
        link_route="# opt",
    )

    route = settings.link_route_string

    assert "PM6" in route
    assert "None" not in route
