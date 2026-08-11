"""Generated-input tests for ORCA electronic-structure settings.

The probe system is a relativistic open-shell nickel complex, deliberately not
the system of any benchmark case: these constructs are required for any
relativistic, open-shell, or benchmark-quality correlated calculation.
"""

import io

import pytest

from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.writer import ORCAInputWriter


def _writer(settings):
    writer = ORCAInputWriter.__new__(ORCAInputWriter)
    writer.settings = settings
    return writer


def _blocks(settings):
    writer = _writer(settings)
    buffer = io.StringIO()
    writer._write_basis_block(buffer)
    writer._write_method_block(buffer)
    writer._write_scf_block(buffer)
    return buffer.getvalue()


def _relativistic_open_shell_settings(**overrides):
    values = dict(
        jobtype="sp",
        ab_initio="DLPNO-CCSD(T1)",
        basis="cc-pVQZ-DK",
        aux_basis="AutoAux",
        relativistic="dkh2",
        reference="rohf",
        frozen_core="fc_electrons",
        frozen_core_electrons=10,
        ri_approximation="none",
        defgrid="grid6",
        scf_tol="verytight",
        heavy_elements=["Ni"],
        heavy_elements_basis="cc-pwCVQZ-DK",
        charge=2,
        multiplicity=3,
    )
    values.update(overrides)
    return ORCAJobSettings(**values)


def test_route_carries_relativistic_hamiltonian_and_ri_choice():
    route = _relativistic_open_shell_settings().route_string
    assert route.startswith("!")
    assert "DKH2" in route
    assert "NoRI" in route
    assert "DLPNO-CCSD(T1)" in route
    assert "cc-pVQZ-DK" in route
    assert "AutoAux" in route


def test_basis_block_assigns_a_per_element_basis():
    assert '  NewGTO Ni "cc-pwCVQZ-DK" end\n' in _blocks(
        _relativistic_open_shell_settings()
    )


def test_method_block_carries_frozen_core_policy():
    rendered = _blocks(_relativistic_open_shell_settings())
    assert "%method\n" in rendered
    assert "  FrozenCore FC_ELECTRONS\n" in rendered
    assert "  NCore 10\n" in rendered


def test_scf_block_carries_the_reference_determinant():
    assert "  HFTyp ROHF\n" in _blocks(_relativistic_open_shell_settings())


def test_scf_block_opens_for_a_reference_without_convergence_settings():
    settings = ORCAJobSettings(
        jobtype="sp",
        functional="bp86",
        basis="def2-svp",
        reference="uhf",
        charge=0,
        multiplicity=3,
    )
    rendered = _blocks(settings)
    assert rendered.startswith("%scf\n")
    assert "  HFTyp UHF\n" in rendered


def test_blocks_are_absent_when_nothing_is_requested():
    settings = ORCAJobSettings(
        jobtype="sp", functional="bp86", basis="def2-svp", charge=0,
        multiplicity=1,
    )
    rendered = _blocks(settings)
    assert "%basis" not in rendered
    assert "%method" not in rendered


def test_restricted_closed_shell_reference_is_refused_for_an_open_shell():
    with pytest.raises(ValueError, match="cannot represent multiplicity"):
        ORCAJobSettings(
            jobtype="sp", functional="bp86", basis="def2-svp",
            reference="rhf", charge=0, multiplicity=3,
        )


def test_relativistic_hamiltonian_requires_a_recontracted_basis():
    with pytest.raises(ValueError, match="recontracted basis"):
        ORCAJobSettings(
            jobtype="sp", functional="bp86", basis="def2-TZVPP",
            relativistic="dkh2", charge=0, multiplicity=1,
        )


@pytest.mark.parametrize(
    "basis", ["DKH-def2-TZVPP", "SARC-DKH-TZVP", "ZORA-def2-SVP", "cc-pVQZ-DK"]
)
def test_recognised_relativistic_basis_families_are_accepted(basis):
    settings = ORCAJobSettings(
        jobtype="sp", functional="bp86", basis=basis, relativistic="dkh2",
        charge=0, multiplicity=1,
    )
    assert "DKH2" in settings.route_string


def test_frozen_core_electron_count_requires_the_matching_policy():
    with pytest.raises(ValueError, match="fc_electrons"):
        ORCAJobSettings(
            jobtype="sp", functional="bp86", basis="def2-svp",
            frozen_core_electrons=10, charge=0, multiplicity=1,
        )


@pytest.mark.parametrize(
    "field,value",
    [
        ("relativistic", "sr-zora-but-not-really"),
        ("reference", "casscf"),
        ("frozen_core", "freeze_everything"),
        ("ri_approximation", "maybe"),
    ],
)
def test_unknown_keyword_values_are_refused_not_passed_through(field, value):
    """The route must not become a free-text channel."""

    with pytest.raises(ValueError, match="Unsupported"):
        ORCAJobSettings(
            jobtype="sp", functional="bp86", basis="def2-svp", charge=0,
            multiplicity=1, **{field: value},
        )
