"""An optimiser cap ORCA honours and the hub could not say.

ORCA stops a geometry optimisation at max(3N, 60) steps of its own
accord. A live goal met that cap three times on a floppy 21-atom Fe(II)
ammine, was told by the repair menu to "consider the optimiser settings
the project exposes", found none, and recorded the conclusion in its own
decision record: "the project surface exposes no geometry-convergence or
maxiter parameter (the cycle-2 attempt to pass MAXITER 500 was rejected
as an unrecognized input keyword), so no legal lever exists."

That attempt is the whole argument for closing the hole. MaxIter is a
``%geom`` setting, not a route-line keyword, so passing it through
``additional_route_parameters`` produced native input ORCA refused, and
two engine calls died in one second. The model did not then invent
native input -- the hub invariant held -- it re-planned around the gap
and paid for it in budget. A hub is only defensible if it can express
what the science needs.

What is pinned: the two fields reach the program, come back out, and
never collide with ORCA's other MaxIter.
"""

from __future__ import annotations

import io
import tempfile
from pathlib import Path

import pytest

from chemsmart.io.orca.input import ORCAInput
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.writer import ORCAInputWriter


class _Molecule:
    frozen_atoms = None
    is_monoatomic = False


class _Job:
    molecule = _Molecule()


def _written(**overrides):
    settings = ORCAJobSettings.default().copy()
    settings.jobtype = "opt"
    settings.charge = 0
    settings.multiplicity = 1
    settings.functional = "B3LYP"
    settings.basis = "def2-svp"
    for name, value in overrides.items():
        setattr(settings, name, value)

    writer = ORCAInputWriter.__new__(ORCAInputWriter)
    writer.settings = settings
    writer.job = _Job()
    buffer = io.StringIO()
    buffer.write(settings.route_string + "\n")
    writer._write_geom_block(buffer)
    writer._write_scf_block(buffer)
    buffer.write("* xyz 0 1\nO 0.0 0.0 0.0\n*\n")
    path = Path(tempfile.mkdtemp()) / "written.inp"
    path.write_text(buffer.getvalue())
    return buffer.getvalue(), ORCAInput(str(path))


@pytest.mark.capability("setting:orca:geom_maxiter")
def test_the_iteration_cap_is_written_into_geom_and_read_back():
    text, parsed = _written(geom_maxiter=300)
    assert "%geom" in text
    assert "MaxIter 300" in text
    assert parsed.geom_maxiter == 300
    # The route line is where the failed attempt put it, and ORCA
    # refused that input; the block is the only place it belongs.
    assert "maxiter" not in text.splitlines()[0].lower()


@pytest.mark.capability("setting:orca:geom_maxiter")
def test_orcas_two_maxiters_never_answer_for_each_other():
    """One caps the SCF, the other the optimiser; the reader separates them.

    Before the reader was block-scoped it searched the whole file for
    the word, so writing a second MaxIter anywhere would have made an
    optimiser cap answer to ``scf_maxiter``.
    """

    _text, parsed = _written(geom_maxiter=300, scf_maxiter=250)
    assert parsed.geom_maxiter == 300
    assert parsed.scf_maxiter == 250

    _text, only_scf = _written(scf_maxiter=250)
    assert only_scf.geom_maxiter is None
    assert only_scf.scf_maxiter == 250

    _text, only_geom = _written(geom_maxiter=300)
    assert only_geom.geom_maxiter == 300
    assert only_geom.scf_maxiter is None


@pytest.mark.capability("setting:orca:opt_convergence")
def test_the_convergence_preset_is_orcas_own_route_word():
    text, _parsed = _written(opt_convergence="tight")
    assert "TightOpt" in text.splitlines()[0]
    text, _parsed = _written(opt_convergence="loose")
    assert "LooseOpt" in text.splitlines()[0]
    # ORCA's default preset has no keyword, so stating it changes
    # nothing about the input rather than writing a word ORCA rejects.
    text, _parsed = _written(opt_convergence="normal")
    assert "opt" in text.splitlines()[0].lower()
    assert "normalopt" not in text.splitlines()[0].lower()


@pytest.mark.capability("setting:orca:opt_convergence")
def test_an_unknown_preset_is_refused_where_it_is_written():
    with pytest.raises(ValueError, match="opt_convergence"):
        ORCAJobSettings.default().copy().__class__(opt_convergence="verytight")
    with pytest.raises(ValueError, match="geom_maxiter"):
        ORCAJobSettings.default().copy().__class__(geom_maxiter=0)


@pytest.mark.capability("setting:orca:geom_maxiter")
def test_only_one_geom_block_is_ever_opened():
    """ORCA reads one %geom block, so the cap rides inside whichever opens.

    A second block would be exactly the class of malformed native input
    this change exists to stop the model from having to hand-write.
    """

    settings = ORCAJobSettings.default().copy()
    settings.jobtype = "modred"
    settings.charge = 0
    settings.multiplicity = 1
    settings.functional = "B3LYP"
    settings.basis = "def2-svp"
    settings.geom_maxiter = 300
    settings.modred = [[1, 2]]

    writer = ORCAInputWriter.__new__(ORCAInputWriter)
    writer.settings = settings
    writer.job = _Job()
    buffer = io.StringIO()
    writer._write_geom_block(buffer)
    writer._write_modred_block(buffer)
    text = buffer.getvalue()

    assert text.count("%geom") == 1, text
    assert "MaxIter 300" in text


@pytest.mark.capability("setting:orca:geom_maxiter")
def test_the_cap_survives_the_verifiers_own_round_trip(tmp_path):
    """The preview verifier reads a written input back through
    ORCAJobSettings.from_inpfile; that path dropped geom_maxiter for
    every stage (REACH-1 po3: a ts project validated with the cap,
    refused as a semantic mismatch, and the session deleted the cap it
    needed)."""

    from chemsmart.jobs.orca.settings import ORCATSJobSettings

    text, _parsed = _written(geom_maxiter=300)
    written = tmp_path / "opt.inp"
    written.write_text(text)
    assert ORCAJobSettings.from_inpfile(str(written)).geom_maxiter == 300

    saddle = tmp_path / "ts.inp"
    saddle.write_text(
        "! OptTS Freq B3LYP def2-svp\n"
        "%geom\n  Calc_Hess True\n  MaxIter 150\nend\n"
        "* xyz 0 1\nO 0.0 0.0 0.0\n*\n"
    )
    read_back = ORCAJobSettings.from_inpfile(str(saddle))
    assert isinstance(read_back, ORCATSJobSettings)
    assert read_back.geom_maxiter == 150
    assert read_back.jobtype == "ts"
