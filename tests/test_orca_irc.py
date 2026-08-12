from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.agent.program_verifiers import _settings_match
from chemsmart.analysis.result_readers import reader_for
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.jobs.orca.settings import ORCAIRCJobSettings, ORCAJobSettings
from chemsmart.jobs.orca.writer import ORCAInputWriter
from chemsmart.settings.orca import ORCAProjectSettings
from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES


def test_orca_irc_project_materializes_and_round_trips_direction(tmp_path):
    project_path = tmp_path / "reaction-path.yaml"
    project_path.write_text(
        "irc:\n"
        "  functional: b3lyp\n"
        "  basis: def2-SVP\n"
        "  direction: both\n"
        "  inithess: calc_anfreq\n"
        "  maxiter: 25\n",
        encoding="utf-8",
    )

    settings = ORCAProjectSettings.from_project(project_path).irc_settings()

    assert isinstance(settings, ORCAIRCJobSettings)
    assert settings.direction == "both"
    assert settings.inithess == "calc_anfreq"
    assert settings.maxiter == 25
    # The CLI supplies molecular state from the input structure/job options;
    # this unit probe materializes the same post-merge settings directly.
    settings.charge = 0
    settings.multiplicity = 1

    molecule = Molecule(
        symbols=["H", "H"],
        positions=np.array([[0.0, 0.0, -0.37], [0.0, 0.0, 0.37]]),
    )
    job = SimpleNamespace(
        settings=settings,
        molecule=molecule,
        jobrunner=SimpleNamespace(num_cores=2, mem_gb=2),
        folder=str(tmp_path),
        label="reaction_path",
    )
    ORCAInputWriter(job).write()

    input_path = tmp_path / "reaction_path.inp"
    rendered = input_path.read_text(encoding="utf-8")
    assert "%irc" in rendered
    assert "direction both" in rendered.casefold()
    assert "inithess calc_anfreq" in rendered.casefold()
    assert "maxiter 25" in rendered.casefold()

    observed = ORCAJobSettings.from_filepath(str(input_path))
    assert isinstance(observed, ORCAIRCJobSettings)
    assert observed.direction == "both"
    assert (
        _settings_match(
            observed, {"direction": "both"}, native_input=input_path
        )
        == []
    )


def test_orca_irc_direction_has_one_explicit_native_vocabulary(tmp_path):
    capability = PROGRAM_CAPABILITIES["orca"]

    assert "direction" in capability.project_owned_parameters
    assert dict(capability.project_parameter_domains)["direction"] == (
        "backward",
        "both",
        "down",
        "forward",
    )
    assert ORCAIRCJobSettings(direction=" FORWARD ").direction == "forward"
    assert ORCAIRCJobSettings(direction=None).direction is None
    with pytest.raises(ValueError, match="ORCA IRC direction"):
        ORCAIRCJobSettings(direction="reverse")

    invalid_project = tmp_path / "reverse-is-not-orca.yaml"
    invalid_project.write_text(
        "irc:\n  functional: b3lyp\n  basis: def2-SVP\n  direction: reverse\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="ORCA IRC direction"):
        ORCAProjectSettings.from_project(invalid_project)


def test_orca_output_exposes_explicit_irc_direction_to_typed_reader(tmp_path):
    output_path = tmp_path / "reaction_path.out"
    output_path.write_text(
        "|  1> ! IRC B3LYP def2-SVP\n"
        "|  2> %irc\n"
        "|  3>   Direction both\n"
        "|  4> end\n"
        "****ORCA TERMINATED NORMALLY****\n",
        encoding="utf-8",
    )

    output = ORCAOutput(str(output_path))

    assert output.irc_direction == "both"
    assert reader_for("orca").read(output, "irc_direction") == ("both", "")
