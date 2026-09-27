"""A flip of named centres is a typed request, written in ORCA's numbering.

``broken_symmetry: true`` mixes one pair of orbitals on a singlet, which
reaches the antiferromagnetic state of one-electron sites (R10 Q18: H2,
twisted ethylene, p-benzyne). On ino2's Ni(II)2, two S = 1 centres, it
landed 33.40 mEh above that state with the spin on the O and Cl bridges,
while ORCA's FlipSpin of one nickel from the high-spin quintet reached it
(R11 truth-2, CUHK 2157086: -3890.968544659585 Eh, <S**2> 1.998044). Every
FlipSpin word was refused and sent to the singlet mixing guess, and the
archived ``FlipSpin 1,2`` had flipped one nickel and the bridging oxygen,
because ORCA counts atoms from 0.

``site_spin_flip`` names the centres by the host's own 1-based atom
numbers and the Ms asked for; the bound multiplicity is the high-spin
state the flip starts from. The host writes HFTyp UHF, FlipSpin in ORCA's
numbering and FinalMs, and reads the written input back to the request.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from .gaussian_fake_preview import fake_preview
from .test_a_native_word_has_its_typed_setting import _refusal

pytestmark = pytest.mark.capability("setting:orca:site_spin_flip")

#: ino2's hydroxo- and chloro-bridged dinickel core, truth-2's oracle
#: geometry (dinickel-oh-cl.xyz, sha256 be1a5c68...): atoms 1 and 2 are Ni.
NI2 = """29
dinickel-oh-cl
Ni -1.54500000 0.00000000 0.00000000
Ni 1.54500000 0.00000000 0.00000000
O 0.00000000 1.32890965 0.00000000
H 0.00000000 2.22953204 0.36024896
Cl 0.00000000 -1.76483928 0.00000000
N -3.66500000 0.00000000 0.00000000
H -4.00548300 0.96149432 0.00000000
H -4.00548300 -0.48074716 -0.83267851
H -4.00548300 -0.48074716 0.83267851
N -1.36022983 0.00000000 -2.11193276
H -1.80947255 -0.83267851 -2.49301999
H -1.80947255 0.83267851 -2.49301999
H -0.37271923 0.00000000 -2.36732036
N -1.36022983 -0.00000000 2.11193276
H -1.80947255 0.83267851 2.49301999
H -1.80947255 -0.83267851 2.49301999
H -0.37271923 -0.00000000 2.36732036
N 3.66500000 0.00000000 0.00000000
H 4.00548300 -0.96149432 0.00000000
H 4.00548300 0.48074716 -0.83267851
H 4.00548300 0.48074716 0.83267851
N 1.36022983 -0.00000000 -2.11193276
H 2.28839032 0.00000000 -2.53491987
H 0.85163700 -0.83267851 -2.40922024
H 0.85163700 0.83267851 -2.40922024
N 1.36022983 0.00000000 2.11193276
H 0.85163700 -0.83267851 2.40922024
H 2.28839032 0.00000000 2.53491987
H 0.85163700 0.83267851 2.40922024
"""

LEVEL = {"functional": "b3lyp", "basis": "def2-svp"}
FLIP_SECOND_NICKEL = {"atoms": [2], "final_ms": 0}


def _scf_block(written):
    lines = [line.strip() for line in written.splitlines()]
    start = lines.index("%scf")
    return lines[start + 1 : lines.index("end", start)]


def test_the_host_writes_the_flip_in_orcas_numbering_and_reads_it_back(
    tmp_path,
):
    """Through the planning session's own preview chain: the project is
    validated, the public ``run --fake`` writes the input, and the
    preview verifier reads it back against the request."""

    from chemsmart.settings.capabilities import PROJECT_OWNED_PARAMETERS

    assert "site_spin_flip" in PROJECT_OWNED_PARAMETERS["orca"]

    receipt, written = fake_preview(
        tmp_path,
        "orca",
        {"gas": {**LEVEL, "site_spin_flip": FLIP_SECOND_NICKEL}},
        NI2,
        (2, 5),
        "sp",
    )

    assert receipt.status == "valid", receipt
    scf = _scf_block(written)
    assert "HFTyp UHF" in scf
    # Host atom 2, the second nickel, is ORCA's atom 1.
    assert "FlipSpin 1" in scf and "FinalMs 0.0" in scf
    assert not any("guessmix" in line.lower() for line in scf)
    assert "* xyz 2 5" in [line.strip() for line in written.splitlines()]


def _run(tmp_path, program, sections, state):
    """``chemsmart run --fake`` of one sp stage on the Ni2 geometry."""

    from chemsmart.agent.live_session import _preview_server_profile
    from chemsmart.cli.main import entry_point

    xyz = tmp_path / "ni2.xyz"
    xyz.write_text(NI2, encoding="utf-8")
    project = tmp_path / "project.yaml"
    project.write_text(yaml.safe_dump(sections), encoding="utf-8")
    server = tmp_path / "server.yaml"
    server.write_text(_preview_server_profile(), encoding="utf-8")
    charge, multiplicity = state
    with CliRunner().isolated_filesystem(temp_dir=tmp_path) as cwd:
        result = CliRunner().invoke(
            entry_point,
            [
                "run",
                "--server",
                str(server),
                "--fake",
                "--no-scratch",
                program,
                "--project",
                str(project),
                "--filename",
                str(xyz),
                "--charge",
                str(charge),
                "--multiplicity",
                str(multiplicity),
                "sp",
            ],
        )
        written = [path.read_text() for path in Path(cwd).rglob("*.inp")]
    return result, written


@pytest.mark.parametrize(
    "flip, multiplicity",
    [
        # A singlet has no high-spin state to flip from.
        (FLIP_SECOND_NICKEL, 1),
        # The quintet reaches Ms 2, 1, 0, -1 and -2: 2 flips nothing,
        # 3 and 0.5 are not among them.
        ({"atoms": [2], "final_ms": 2}, 5),
        ({"atoms": [2], "final_ms": 3}, 5),
        ({"atoms": [2], "final_ms": 0.5}, 5),
        # The geometry has 29 atoms.
        ({"atoms": [30], "final_ms": 0}, 5),
    ],
)
def test_a_flip_the_state_cannot_carry_is_refused_by_name(
    tmp_path, flip, multiplicity
):
    result, written = _run(
        tmp_path,
        "orca",
        {"gas": {**LEVEL, "site_spin_flip": flip}},
        (2, multiplicity),
    )
    assert result.exit_code != 0
    _refused_by_the_request(result)
    assert not written


def _refused_by_the_request(result):
    """The request's own check refused it, not a loader that lacks it."""

    message = str(result.exception)
    assert "site_spin_flip" in message
    assert "is not in list of keywords" not in message


def test_one_state_gets_one_broken_symmetry_mechanism(tmp_path):
    result, written = _run(
        tmp_path,
        "orca",
        {
            "gas": {
                **LEVEL,
                "broken_symmetry": True,
                "site_spin_flip": FLIP_SECOND_NICKEL,
            }
        },
        (2, 5),
    )
    assert result.exit_code != 0
    _refused_by_the_request(result)
    assert not written


@pytest.mark.parametrize(
    "program, section",
    [("gaussian", "gas"), ("pyscf", "sp")],
)
def test_a_program_that_writes_no_site_flip_refuses_it(
    tmp_path, program, section
):
    """Gaussian and PySCF have no named-site flip ChemSmart writes; the
    loader refuses the key rather than run the unflipped state."""

    level = (
        {"functional": "b3lyp", "basis": "def2svp"}
        if program == "gaussian"
        else dict(LEVEL)
    )
    result, _written = _run(
        tmp_path,
        program,
        {section: {**level, "site_spin_flip": FLIP_SECOND_NICKEL}},
        (2, 5),
    )
    assert result.exit_code != 0
    assert "site_spin_flip" in str(result.exception)


def test_a_native_flip_word_names_the_typed_request():
    report = _refusal(
        "orca",
        {"gas": {**LEVEL, "additional_route_parameters": "FlipSpin 1,2"}},
    )
    assert "site_spin_flip" in report["route"]
