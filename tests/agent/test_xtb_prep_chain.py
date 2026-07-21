"""xTB pre-optimization -> heavy DFT job preparation chain."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from chemsmart.agent.tools import (
    build_gaussian_settings,
    build_job,
    build_molecule,
    build_xtb_settings,
    dry_run_input,
    extract_optimized_geometry,
    run_local,
)
from chemsmart.agent.tools_fs import save_geometry

_XTB_FIXTURES = Path(__file__).parents[1] / "data" / "XTBTests" / "outputs"


def _xtb_job_over_fixture(tmp_path: Path, fixture_name: str):
    """Build an XTBJob whose folder is a copy of a REAL xtb output folder."""

    source = _XTB_FIXTURES / fixture_name
    folder = tmp_path / fixture_name
    shutil.copytree(source, folder)
    geometry = next(
        path
        for path in folder.iterdir()
        if path.suffix == ".xyz" and not path.name.startswith("xtbopt")
    )
    job = build_job(
        "xtb.opt",
        molecule=build_molecule(str(geometry)),
        settings=build_xtb_settings(charge=0, multiplicity=1),
        label=fixture_name,
    )
    job.set_folder(str(folder))
    return job


class TestExtractXtbGeometry:
    def test_extracts_real_xtbopt_geometry_with_grounding(
        self, tmp_path
    ) -> None:
        # water_ohess carries a genuine xtbopt.xyz written by real xtb 6.7.1;
        # synthetic stand-ins hid two parser bugs before, so the fixture must
        # stay real output.
        job = _xtb_job_over_fixture(tmp_path, "water_ohess")

        molecule = extract_optimized_geometry(job)

        assert list(molecule.symbols) == ["O", "H", "H"]
        assert molecule.charge == 0
        assert molecule.multiplicity == 1
        source = getattr(molecule, "_agent_source_filepath")
        assert source.endswith("xtbopt.xyz")
        assert Path(source).is_file()

    def test_missing_xtbopt_raises_actionable_error(self, tmp_path) -> None:
        job = _xtb_job_over_fixture(tmp_path, "water_ohess")
        (Path(job.folder) / "xtbopt.xyz").unlink()

        with pytest.raises(FileNotFoundError, match="xtbopt"):
            extract_optimized_geometry(job)


class TestSummarizeXtbLocalOutput:
    def test_real_hessian_output_reports_energy_and_frequencies(
        self, tmp_path
    ) -> None:
        # run_local's summary previously handled only Gaussian/ORCA, so a
        # genuine xtb run returned an empty summary and the agent could not
        # report the energy the prep chain depends on. Verify against real
        # xtb 6.7.1 output, not a synthetic string.
        from chemsmart.agent.services.local_runtime import (
            _summarize_xtb_output,
        )

        job = _xtb_job_over_fixture(tmp_path, "water_ohess")
        # The runner names the output <label>.out; point the job's label at
        # the fixture's real output file.
        job.label = "water_ohess"

        summary = _summarize_xtb_output(job)

        assert summary["converged"] is True
        assert summary["energy"] is not None
        assert summary["energy"] < 0
        assert summary["frequency_count"] == 3
        assert summary["imag_freqs"] == []


class TestXtbPrepChain:
    def test_fake_xtb_preopt_feeds_a_grounded_gaussian_job(
        self,
        monkeypatch,
        tmp_path,
        single_molecule_xyz_file,
        xtb_jobrunner_no_scratch,
        gaussian_jobrunner_no_scratch,
    ) -> None:
        """build_molecule -> xtb.opt (fake run) -> extract -> save_geometry
        -> gaussian.opt dry run, with the refined file grounding `-f`."""

        monkeypatch.chdir(tmp_path)
        molecule = build_molecule(single_molecule_xyz_file)
        xtb_job = build_job(
            "xtb.opt",
            molecule=molecule,
            settings=build_xtb_settings(charge=0, multiplicity=1),
            label="prep",
            jobrunner=xtb_jobrunner_no_scratch,
        )
        xtb_job.set_folder(str(tmp_path / "prep"))

        run_result = run_local(xtb_job)

        assert run_result["ok"] is True
        assert (tmp_path / "prep" / "xtbopt.xyz").is_file()

        refined = extract_optimized_geometry(xtb_job)
        saved = save_geometry(refined, "prep_refined.xyz")
        assert "error" not in saved

        gaussian_job = build_job(
            "gaussian.opt",
            molecule=refined,
            settings=build_gaussian_settings(
                "B3LYP", "6-31G*", charge=0, multiplicity=1
            ),
            label="dft_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        gaussian_job.set_folder(str(tmp_path / "dft"))

        dry_run = dry_run_input(gaussian_job)

        assert dry_run["cli_grounded"] is True
        assert "prep_refined.xyz" in dry_run["command"]
        assert "-f" in dry_run["command"]
