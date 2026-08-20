"""Generated-input tests for ORCA electronic-structure settings.

The probe system is a relativistic open-shell nickel complex, deliberately not
the system of any benchmark case: these constructs are required for any
relativistic, open-shell, or benchmark-quality correlated calculation.
"""

import io

import pytest

from chemsmart.io.orca.route import ORCARoute
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.writer import ORCAInputWriter
from chemsmart.settings.orca import ORCAProjectSettings


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
        jobtype="sp",
        functional="bp86",
        basis="def2-svp",
        charge=0,
        multiplicity=1,
    )
    rendered = _blocks(settings)
    assert "%basis" not in rendered
    assert "%method" not in rendered


def test_restricted_closed_shell_reference_is_refused_for_an_open_shell():
    with pytest.raises(ValueError, match="cannot represent multiplicity"):
        ORCAJobSettings(
            jobtype="sp",
            functional="bp86",
            basis="def2-svp",
            reference="rhf",
            charge=0,
            multiplicity=3,
        )


def test_relativistic_hamiltonian_requires_a_recontracted_basis():
    with pytest.raises(ValueError, match="recontracted basis"):
        ORCAJobSettings(
            jobtype="sp",
            functional="bp86",
            basis="def2-TZVPP",
            relativistic="dkh2",
            charge=0,
            multiplicity=1,
        )


@pytest.mark.parametrize(
    "basis", ["DKH-def2-TZVPP", "SARC-DKH-TZVP", "ZORA-def2-SVP", "cc-pVQZ-DK"]
)
def test_recognised_relativistic_basis_families_are_accepted(basis):
    settings = ORCAJobSettings(
        jobtype="sp",
        functional="bp86",
        basis=basis,
        relativistic="dkh2",
        charge=0,
        multiplicity=1,
    )
    assert "DKH2" in settings.route_string


def test_frozen_core_electron_count_requires_the_matching_policy():
    with pytest.raises(ValueError, match="fc_electrons"):
        ORCAJobSettings(
            jobtype="sp",
            functional="bp86",
            basis="def2-svp",
            frozen_core_electrons=10,
            charge=0,
            multiplicity=1,
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
            jobtype="sp",
            functional="bp86",
            basis="def2-svp",
            charge=0,
            multiplicity=1,
            **{field: value},
        )


@pytest.mark.parametrize("aux_basis", [None, "def2/J", "def2/JK"])
def test_dlpno_coupled_cluster_requires_correlation_fitting_basis(aux_basis):
    with pytest.raises(ValueError, match="correlation-fitting AuxC"):
        ORCAJobSettings.from_dict(
            {
                "jobtype": "sp",
                "ab_initio": "DLPNO-CCSD(T)",
                "basis": "def2-TZVP",
                "aux_basis": aux_basis,
                "charge": 0,
                "multiplicity": 1,
            }
        )


def test_dlpno_auxiliary_basis_error_refuses_silent_autoaux_injection():
    with pytest.raises(ValueError, match="will not silently add AutoAux"):
        ORCAJobSettings(
            jobtype="sp",
            ab_initio="DLPNO-CCSD(T)",
            basis="def2-TZVP",
            charge=0,
            multiplicity=1,
        )


def test_dlpno_project_yaml_fails_before_native_execution_without_auxc(
    tmp_path,
):
    project = tmp_path / "dlpno-missing-auxc.yaml"
    project.write_text(
        "sp:\n"
        "  ab_initio: DLPNO-CCSD(T)\n"
        "  basis: def2-TZVP\n"
        "  charge: 0\n"
        "  multiplicity: 1\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="aux_basis is missing"):
        ORCAProjectSettings.from_project(str(project)).sp_settings()


def test_dlpno_explicit_auxc_must_match_orbital_basis():
    with pytest.raises(ValueError, match="same family and zeta level"):
        ORCAJobSettings(
            jobtype="sp",
            ab_initio="DLPNO-CCSD(T)",
            basis="def2-TZVP",
            aux_basis="def2-SVP/C",
            charge=0,
            multiplicity=1,
        )


def test_dlpno_auxc_compatibility_comes_from_exact_registered_pairing():
    """A plausible suffix must not turn an unknown name into an AuxC role."""

    with pytest.raises(ValueError, match="cannot be qualified"):
        ORCAJobSettings(
            jobtype="sp",
            ab_initio="DLPNO-CCSD(T)",
            basis="invented-TZVP",
            aux_basis="invented-TZVP/C",
            charge=0,
            multiplicity=1,
        )


@pytest.mark.parametrize(
    ("basis", "aux_basis", "role"),
    [
        ("def2-TZVP", "def2-TZVP/C", "correlation"),
        ("cc-pVTZ-F12", "cc-pVTZ-F12-MP2fit", "correlation"),
        ("def2-TZVP", "AutoAux", "autoaux"),
    ],
)
def test_dlpno_correlation_auxiliary_basis_round_trips(basis, aux_basis, role):
    settings = ORCAJobSettings(
        jobtype="sp",
        ab_initio="DLPNO-CCSD(T1)",
        basis=basis,
        aux_basis=aux_basis,
        charge=0,
        multiplicity=1,
    )
    parsed = ORCARoute(settings.route_string)

    assert parsed.ab_initio == "dlpno-ccsd(t1)"
    assert parsed.auxiliary_basis == aux_basis.casefold()
    assert parsed.auxiliary_basis_role == role
    assert parsed.correlation_auxiliary_basis == aux_basis.casefold()
