import io
import os
import re
import shutil
from filecmp import cmp
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.orca import (
    ORCAModredJob,
    ORCAOptJob,
    ORCAScanJob,
    ORCASinglePointJob,
    ORCATSJob,
)
from chemsmart.jobs.orca.neb import ORCANEBJob
from chemsmart.jobs.orca.settings import (
    ORCAIRCJobSettings,
    ORCAJobSettings,
    ORCANEBJobSettings,
    ORCAQMMMJobSettings,
    ORCATSJobSettings,
)
from chemsmart.jobs.orca.writer import ORCAInputWriter
from chemsmart.settings.orca import ORCAProjectSettings


def _writer(settings, molecule=None):
    job = SimpleNamespace(
        settings=settings,
        molecule=(
            molecule if molecule is not None else MagicMock(frozen_atoms=None)
        ),
        jobrunner=MagicMock(num_cores=4, mem_gb=8),
        label="test_job",
        folder="/tmp",
    )
    return ORCAInputWriter(job=job)


class TestORCAInputWriter:
    def test_write_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
        orca_written_opt_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )

        assert isinstance(job, ORCAOptJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_opt.inp")
        assert os.path.isfile(orca_file)
        assert cmp(
            orca_file, orca_written_opt_file, shallow=False
        )  # writes input file as expected

        # job run will result in the job being run and
        # the output file copied back to run folder
        # job.run()
        # assert job.is_complete()

    def test_write_opt_job_with_route(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_opt_file_with_route,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.route_to_be_written = "pbepbe/6-31g(d) opt"
        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_opt_with_route",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAOptJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_opt_with_route.inp")
        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_opt_file_with_route, shallow=False)

    def test_write_modred_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_modred_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.modred_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = [[1, 2], [3, 4, 5]]
        job = ORCAModredJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_modred",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAModredJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_modred.inp")
        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_modred_file, shallow=False)

    def test_write_scan_single_degree_of_freedom_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_scan_single_degree_of_freedom_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [[1, 2]],
            "num_steps": [10],
            "dist_start": [1.5],
            "dist_end": [3.5],
        }
        job = ORCAScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_scan",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAScanJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_scan.inp")
        assert os.path.isfile(orca_file)
        assert cmp(
            orca_file,
            orca_written_scan_single_degree_of_freedom_file,
            shallow=False,
        )

    def test_write_scan_multiple_degrees_of_freedom_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_scan_multiple_degrees_of_freedom_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [[1, 2], [3, 4, 5], [6, 7, 8, 9]],
            "num_steps": [10, 15, 20],
            "dist_start": [1.5, 70.0, 80.0],
            "dist_end": [3.5, 85.0, 60.0],
        }
        job = ORCAScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_scan",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAScanJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_scan.inp")
        assert os.path.isfile(orca_file)
        assert cmp(
            orca_file,
            orca_written_scan_multiple_degrees_of_freedom_file,
            shallow=False,
        )

    def test_write_scan_multiple_degrees_of_freedom_with_constraints_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_scan_multiple_degrees_of_freedom_with_constraints_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [[1, 2], [3, 4, 5]],
            "num_steps": [10, 15],
            "dist_start": [1.5, 70.0],
            "dist_end": [3.5, 85.0],
            "constrained_coordinates": [[5, 6], [7, 8, 9]],
        }
        job = ORCAScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_scan",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAScanJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_scan.inp")
        assert os.path.isfile(orca_file)
        assert cmp(
            orca_file,
            orca_written_scan_multiple_degrees_of_freedom_with_constraints_file,
            shallow=False,
        )

    def test_write_ts_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_ts_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.ts_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = ORCATSJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_ts",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCATSJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_ts.inp")
        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_ts_file, shallow=False)

    def test_write_opt_input_from_logfile(
        self,
        tmpdir,
        orca_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        orca_jobrunner_no_scratch,
        orca_written_ts_from_nhc_singlet_log_file,
    ):
        """Taking the Gaussian nhc_neutral_singlet.log output
        and write orca input file."""
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        job_settings = ORCAJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        # also merge the title keywords
        keywords = ("charge", "multiplicity", "title")
        ts_settings = ts_settings.merge(job_settings, keywords=keywords)
        job = ORCATSJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=ts_settings,
            label="orca_ts_from_log",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCATSJob)
        orca_writer = ORCAInputWriter(job=job)
        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_ts_from_log.inp")
        assert os.path.isfile(orca_file)
        assert cmp(
            orca_file,
            orca_written_ts_from_nhc_singlet_log_file,
            shallow=False,
        )

    def test_write_sp_input_with_solvation_from_logfile(
        self,
        tmpdir,
        orca_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        orca_jobrunner_no_scratch,
        orca_written_sp_from_nhc_singlet_log_with_solvent_file,
    ):
        """Test writing simple .inp input file using settings from .out file,
        including solvation."""

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        sp_settings = project_settings.sp_settings()
        job_settings = ORCAJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        sp_settings = sp_settings.merge(job_settings)
        job = ORCASinglePointJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=sp_settings,
            label="orca_sp_from_log_with_solvent",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCASinglePointJob)
        orca_writer = ORCAInputWriter(job=job)
        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_sp_from_log_with_solvent.inp")
        assert os.path.isfile(orca_file)
        assert cmp(
            orca_file, orca_written_sp_from_nhc_singlet_log_with_solvent_file
        )

    def test_write_opt_input_on_monoatomic_species(
        self,
        tmpdir,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_he_monoatomic_opt_file,
    ):

        helium = Molecule(
            symbols=["He"],
            positions=[(0.0, 0.0, 0.0)],
        )
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        opt_settings = project_settings.opt_settings()
        opt_settings.charge = 0
        opt_settings.multiplicity = 1
        job = ORCAOptJob.from_molecule(
            molecule=helium,
            settings=opt_settings,
            label="orca_he_monoatomic_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAOptJob)
        orca_writer = ORCAInputWriter(job=job)
        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_he_monoatomic_opt.inp")
        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_he_monoatomic_opt_file)

    def test_write_neb_input(
        self,
        tmpdir,
        orca_input_nebts_reactant_xyz_file,
        orca_input_nebts_product_xyz_file,
        orca_input_nebts_ts_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_neb_file,
    ):

        reactant_xyz_file = tmpdir.join("reactant.xyz")
        product_xyz_file = tmpdir.join("product.xyz")
        ts_xyz_file = tmpdir.join("ts.xyz")
        shutil.copy(orca_input_nebts_reactant_xyz_file, reactant_xyz_file)
        product_xyz_file = os.path.basename(
            shutil.copy(orca_input_nebts_product_xyz_file, product_xyz_file)
        )
        ts_xyz_file = os.path.basename(
            shutil.copy(orca_input_nebts_ts_xyz_file, ts_xyz_file)
        )
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.neb_settings()
        settings.semiempirical = "GFN2-xTB"
        settings.joboption = "NEB-TS"
        settings.jobtype = "neb"
        settings.nimages = 5
        settings.charge = 0
        settings.multiplicity = 1
        settings.ending_xyzfile = product_xyz_file
        settings.intermediate_xyzfile = ts_xyz_file
        settings.preopt_ends = True
        job = ORCANEBJob.from_filename(
            filename=reactant_xyz_file,
            settings=settings,
            label="orca_neb",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCANEBJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_neb.inp")

        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_neb_file, shallow=False)

    def test_smd_solvent_uses_smd_in_route_no_cpcm_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """SMD solvation must write SMD(solvent) in route; no %cpcm block needed."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "smd"
        settings.solvent_id = "water"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_smd_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_smd_opt.inp")
        assert os.path.isfile(orca_file)

        content = open(orca_file).read()
        # Route must use SMD(solvent) — canonical ORCA 6.0 form
        assert "SMD(water)" in content
        assert "CPCM" not in content
        # No %cpcm block needed for pure SMD (route handles activation)
        assert "%cpcm" not in content
        assert "SMD true" not in content
        assert "SMDsolvent" not in content

    def test_cpcm_solvent_uses_cpcm_in_route_no_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Pure CPCM solvation must write CPCM(solvent) in route only, no %cpcm block."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cpcm"
        settings.solvent_id = "toluene"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cpcm_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cpcm_opt.inp")
        assert os.path.isfile(orca_file)

        content = open(orca_file).read()
        # Route must have CPCM(toluene)
        assert "CPCM(toluene)" in content
        # No %cpcm block needed for pure CPCM
        assert "%cpcm" not in content

    def test_smd_with_additional_solvent_options(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """SMD + additional_solvent_options: route SMD(solvent), options in %cpcm block."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "smd"
        settings.solvent_id = "water"
        settings.additional_solvent_options = "SurfaceType gepol_ses"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_smd_extra_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_smd_extra_opt.inp")
        content = open(orca_file).read()

        # Route uses SMD(solvent) — canonical ORCA 6.0 form
        assert "SMD(water)" in content
        assert "CPCM" not in content
        # %cpcm block only for the extra surface-type option
        assert "%cpcm" in content
        assert "SurfaceType gepol_ses" in content
        # SMD activation is handled by the route line, NOT the block
        assert "SMD true" not in content
        assert "SMDsolvent" not in content

    def test_custom_epsilon_cpcm_no_solvent_name(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Custom CPCM with Epsilon/Refrac: route must have bare CPCM (no parens)."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cpcm"
        settings.solvent_id = None
        settings.additional_solvent_options = "Epsilon 78.36\nRefrac 1.33"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_custom_epsilon_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_custom_epsilon_opt.inp")
        assert os.path.isfile(orca_file)

        content = open(orca_file).read()
        # Route must have bare CPCM (no solvent in parentheses) for custom epsilon
        assert re.search(r"\bCPCM\b", content)
        assert "CPCM(" not in content
        # %cpcm block must contain Epsilon and Refrac on separate lines
        assert "%cpcm" in content
        assert "Epsilon 78.36" in content
        assert "Refrac 1.33" in content
        # SMD should not be present
        assert "SMD true" not in content

    def test_custom_epsilon_with_rsolv_and_surface_type(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Multi-line additional_solvent_options: each option on its own line."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cpcm"
        settings.solvent_id = None
        settings.additional_solvent_options = (
            "Epsilon 78.36\nRefrac 1.33\nSurfaceType gepol_ses\nRsolv 1.30"
        )

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_multiline_solvent_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_multiline_solvent_opt.inp")
        content = open(orca_file).read()

        assert "%cpcm" in content
        assert "Epsilon 78.36" in content
        assert "Refrac 1.33" in content
        assert "SurfaceType gepol_ses" in content
        assert "Rsolv 1.30" in content
        # Verify each option is on its own properly indented line
        lines = content.splitlines()
        cpcm_lines = [
            line
            for line in lines
            if line.strip()
            in {
                "Epsilon 78.36",
                "Refrac 1.33",
                "SurfaceType gepol_ses",
                "Rsolv 1.30",
            }
        ]
        assert len(cpcm_lines) == 4
        for line in cpcm_lines:
            assert line.startswith("  ")  # all indented by 2 spaces

    def test_smd_with_surface_type_and_rsolv(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """SMD + SurfaceType + Rsolv: route SMD(solvent), options in %cpcm block."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "smd"
        settings.solvent_id = "water"
        settings.additional_solvent_options = (
            "SurfaceType gepol_ses\nRsolv 1.30"
        )

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_smd_surface_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_smd_surface_opt.inp")
        content = open(orca_file).read()

        # Route uses SMD(solvent) — canonical ORCA 6.0 form
        assert "SMD(water)" in content
        assert "CPCM" not in content
        assert "%cpcm" in content
        # SMD activation is handled by the route line, NOT the block
        assert "SMD true" not in content
        assert "SMDsolvent" not in content
        assert "SurfaceType gepol_ses" in content
        assert "Rsolv 1.30" in content

    def test_custom_solvent_from_yaml_project(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_custom_solv_project_name,
        orca_jobrunner_no_scratch,
    ):
        """custom_solvent from YAML project is written into %cpcm block.

        The ``custom_solv.yaml`` project sets::

            custom_solvent: |
              Epsilon 16.7
              Refrac 1.275

        The route must use bare ``CPCM`` (no named solvent in parens) and the
        ``%cpcm`` block must contain the Epsilon and Refrac lines.
        """
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_custom_solv_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_custom_solv_sp",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_custom_solv_sp.inp")
        assert os.path.isfile(orca_file)

        content = open(orca_file).read()
        # Route must have bare CPCM (no named solvent)
        assert re.search(r"\bCPCM\b", content)
        assert "CPCM(" not in content
        # %cpcm block must contain Epsilon and Refrac
        assert "%cpcm" in content
        assert "Epsilon 16.7" in content
        assert "Refrac 1.275" in content
        # SMD should not be activated
        assert "SMD true" not in content

    def test_custom_solvent_programmatic(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """custom_solvent set programmatically is written into %cpcm block."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cpcm"
        settings.solvent_id = None
        settings.custom_solvent = "Epsilon 16.7\nRefrac 1.275"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_custom_solv_prog",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_custom_solv_prog.inp")
        content = open(orca_file).read()

        assert re.search(r"\bCPCM\b", content)
        assert "CPCM(" not in content
        assert "%cpcm" in content
        assert "Epsilon 16.7" in content
        assert "Refrac 1.275" in content
        assert "SMD true" not in content

    def test_custom_solvent_with_smd_and_additional_options(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """custom_solvent + SMD + additional_solvent_options all appear in %cpcm block.

        With ORCA 6.0: route uses ``SMD(water)``, and the %cpcm block
        contains ``custom_solvent`` lines first, then ``additional_solvent_options``.
        SMD activation is handled by the route line alone.
        """
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "smd"
        settings.solvent_id = "water"
        settings.custom_solvent = "Epsilon 78.36\nRefrac 1.33"
        settings.additional_solvent_options = "SurfaceType gepol_ses"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_custom_smd_combo",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_custom_smd_combo.inp")
        content = open(orca_file).read()

        # Route uses SMD(water) — canonical ORCA 6.0 form
        assert "SMD(water)" in content
        assert "CPCM" not in content
        assert "%cpcm" in content
        # SMD activation is handled by the route line, NOT the block
        assert "SMD true" not in content
        assert "SMDsolvent" not in content
        # custom_solvent lines in block
        assert "Epsilon 78.36" in content
        assert "Refrac 1.33" in content
        # additional_solvent_options in block
        assert "SurfaceType gepol_ses" in content

        # Verify order: custom_solvent before additional_solvent_options
        epsilon_pos = content.index("Epsilon 78.36")
        surface_pos = content.index("SurfaceType gepol_ses")
        assert epsilon_pos < surface_pos

    def test_cpcmc_solvent_uses_cpcmc_in_route_no_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """CPCMC solvation writes CPCMC(solvent) in route and no %cpcm block.

        ``cpcmc`` applies C-PCM with the COSMO epsilon function
        (``!CPCMC(solvent)``).  This is the replacement for the legacy COSMO
        model which was removed in ORCA 4.0.
        """
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cpcmc"
        settings.solvent_id = "water"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cpcmc_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cpcmc_opt.inp")
        content = open(orca_file).read()

        assert "CPCMC(water)" in content
        # No standalone CPCM or old COSMO keywords
        assert "COSMO" not in content
        # No block needed for pure CPCMC with named solvent
        assert "%cpcm" not in content

    def test_cpcmc_custom_epsilon_no_solvent_name(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """CPCMC with custom Epsilon (no solvent_id) writes bare CPCMC + %cpcm block."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cpcmc"
        settings.additional_solvent_options = "Epsilon 16.7\nRefrac 1.275"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cpcmc_custom",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cpcmc_custom.inp")
        content = open(orca_file).read()

        # Bare CPCMC in route (no parens) — use regex to match whole word
        assert re.search(r"\bCPCMC\b(?!\()", content)
        assert "CPCMC(" not in content
        assert "COSMO" not in content
        # %cpcm block with custom parameters
        assert "%cpcm" in content
        assert "Epsilon 16.7" in content
        assert "Refrac 1.275" in content
        assert "end" in content

    def test_cpcmc_custom_solvent_yaml(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """CPCMC with custom_solvent YAML goes into %cpcm block."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cpcmc"
        settings.custom_solvent = "Epsilon 32.6\nRefrac 1.34"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cpcmc_yaml",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cpcmc_yaml.inp")
        content = open(orca_file).read()

        assert "%cpcm" in content
        assert "Epsilon 32.6" in content
        assert "Refrac 1.34" in content
        assert "COSMO" not in content

    def test_cosmors_uses_cosmors_in_route_and_cosmors_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """COSMO-RS writes COSMORS(solvent) in route and %cosmors block."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solvent_id = "water"
        settings.additional_solvent_options = "Temperature 298.15"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cosmors_opt.inp")
        content = open(orca_file).read()

        # Route must use COSMORS(solvent) — canonical ORCA 6.0 openCOSMO-RS form
        assert "COSMORS(water)" in content
        assert "COSMO(" not in content
        assert "CPCM" not in content
        assert "%cosmors" in content
        assert "Temperature 298.15" in content
        assert "%cpcm" not in content

    def test_cosmors_custom_solvent_in_cosmors_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """COSMO-RS with custom_solvent YAML goes into %cosmors block."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.custom_solvent = "Temperature 298.15\nDensity 1.0"

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_custom",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cosmors_custom.inp")
        content = open(orca_file).read()

        assert "%cosmors" in content
        assert "Temperature 298.15" in content
        assert "Density 1.0" in content
        assert "%cpcm" not in content
        assert "COSMO(" not in content

    def test_cosmors_full_params_from_yaml_project(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_custom_solv_cosmors_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Full %cosmors parameter set written from custom_solvent in project YAML.

        The ``custom_solv_cosmors.yaml`` project sets ``solvent_model: cosmors``
        and provides a multi-parameter ``custom_solvent`` block containing
        COSMO-RS parameters (``aeff``, ``lnalpha``, ``lnchb``, ``temp``,
        ``dftfunc``, ``dftbas``, ``solvent``).  All lines must be written into
        the ``%cosmors … end`` block unchanged.
        """
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_custom_solv_cosmors_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_full_params",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cosmors_full_params.inp")
        assert os.path.isfile(orca_file)

        content = open(orca_file).read()

        # Route must use COSMORS(water) — solvent_id from YAML
        assert "COSMORS(water)" in content
        assert "CPCM" not in content
        assert "COSMO(" not in content

        # All custom_solvent lines (except filtered 'solvent') appear in %cosmors block
        assert "%cosmors" in content
        assert "aeff" in content
        assert "lnalpha" in content
        assert "lnchb" in content
        assert "temp" in content
        assert "dftfunc" in content
        assert "dftbas" in content
        assert "end" in content

        # ORCA 6.1 guard: 'solvent "water"' must NOT appear in the block since
        # solvent_id is already encoded in the route as COSMORS(water).
        # Verify the block lines do not contain a bare 'solvent "..."' keyword.
        lines = content.splitlines()
        cosmors_start = next(
            i for i, ln in enumerate(lines) if ln.strip() == "%cosmors"
        )
        cosmors_end = next(
            i
            for i, ln in enumerate(lines)
            if ln.strip() == "end" and i > cosmors_start
        )
        block_lines = lines[cosmors_start + 1 : cosmors_end]
        for line in block_lines:
            stripped_lower = line.strip().lower()
            assert not stripped_lower.startswith("solvent "), (
                f"'solvent name' keyword must not appear in %cosmors block "
                f"when solvent_id is already in the route: {line!r}"
            )
            # All non-empty lines inside the block must be indented
            if line.strip():
                assert line.startswith("  "), f"Line not indented: {line!r}"

        # No %cpcm block should appear
        assert "%cpcm" not in content

    def test_cosmors_programmatic_full_params(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Full %cosmors parameter set set programmatically via custom_solvent.

        Verifies that all ORCA COSMO-RS block parameters can be set via
        ``settings.custom_solvent`` and that ``additional_solvent_options``
        appends further lines in the same ``%cosmors`` block.
        """
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solvent_id = "water"
        # All known %cosmors parameters as custom_solvent
        settings.custom_solvent = (
            "aeff              5.925\n"
            "lnalpha           0.202\n"
            "lnchb             0.166\n"
            "chbt              1.50\n"
            "sigmahb           9.61e-3\n"
            "rav               0.50\n"
            "fcorr             2.40\n"
            "ravcorr           1.00\n"
            "astd              41.624\n"
            "zcoord            10.0\n"
            "dgsolv_eta       -4.4480\n"
            "dgsolv_omegaring  0.2630\n"
            "temp              298.15\n"
            'dftfunc           "BP86"\n'
            'dftbas            "def2-TZVPD"\n'
            # solvent and solventfilename are independent: solvent names the
            # internal database entry; solventfilename points to the .cosmorsxyz
            # file — they can differ (e.g. internal "THF" with a custom file)
            'solvent           "THF"\n'
            'solventfilename   "water"\n'
            "orbs_vac          false"
        )
        settings.additional_solvent_options = 'dftfunc "PBE0"'

        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_all_params",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cosmors_all_params.inp")
        content = open(orca_file).read()

        # Route must use COSMORS(water)
        assert "COSMORS(water)" in content
        assert "%cpcm" not in content

        # Extract the %cosmors block lines for verification
        lines = content.splitlines()
        block_lines = []
        in_block = False
        for line in lines:
            if line.strip() == "%cosmors":
                in_block = True
                continue
            if in_block and line.strip() == "end":
                break
            if in_block:
                block_lines.append(line)
        block_content = "\n".join(block_lines)

        # All non-solvent parameters must appear inside the %cosmors block
        assert "aeff" in block_content
        assert "lnalpha" in block_content
        assert "chbt" in block_content
        assert "sigmahb" in block_content
        assert "rav" in block_content
        assert "fcorr" in block_content
        assert "dgsolv_eta" in block_content
        assert "orbs_vac" in block_content

        # ORCA 6.1 guard: 'solvent "THF"' must be filtered (solvent_id="water" in route)
        assert '  solvent           "THF"' not in content
        # 'solventfilename' is a distinct keyword — must appear inside the block
        assert 'solventfilename   "water"' in block_content

        # additional_solvent_options appended after custom_solvent
        custom_pos = content.index("orbs_vac")
        additional_pos = content.index('dftfunc "PBE0"')
        assert custom_pos < additional_pos

    def test_cosmors_solvent_name_filtered_when_in_route(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """ORCA 6.1 guard: 'solvent "name"' filtered from %cosmors block when route has COSMORS(name).

        When solvent_id is set, the route already encodes the solvent as
        COSMORS(solvent_id).  Writing 'solvent "name"' in the %cosmors block
        causes an ORCA 6.1 INPUT ERROR (DUPLICATED KEYWORD).  The writer must
        filter out those lines automatically — case-insensitively.

        'solventfilename "..."' is a distinct keyword (points to a .cosmorsxyz
        file) and must NOT be filtered even when solvent_id is set.
        """
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solvent_id = "water"
        # Include 'solvent' in all three case variants (all should be filtered),
        # 'solventfilename' (preserved), and regular parameters (preserved).
        settings.custom_solvent = (
            "temp              298.15\n"
            'solvent           "water"\n'  # lowercase — filtered
            'Solvent           "acetone"\n'  # mixed case — also filtered
            'SOLVENT           "toluene"\n'  # uppercase — also filtered
            'solventfilename   "custom_thf"\n'  # different keyword — NOT filtered
            'dftfunc           "BP86"\n'
        )

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_filter_test",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cosmors_filter_test.inp")
        content = open(orca_file).read()

        # Route must use COSMORS(water)
        assert "COSMORS(water)" in content

        # Extract the %cosmors block lines for verification
        lines = content.splitlines()
        block_lines = []
        in_block = False
        for line in lines:
            if line.strip() == "%cosmors":
                in_block = True
                continue
            if in_block and line.strip() == "end":
                break
            if in_block:
                block_lines.append(line)

        # %cosmors block is still written (has non-solvent content)
        assert block_lines, "Expected %cosmors block to be present"
        block_content = "\n".join(block_lines)
        assert "temp" in block_content
        assert "dftfunc" in block_content

        # All three case variants of 'solvent "..."' must be filtered
        for line in block_lines:
            stripped_lower = line.strip().lower()
            assert not stripped_lower.startswith(
                "solvent "
            ), f"'solvent name' (any case) must not appear in %cosmors block: {line!r}"

        # 'solventfilename' starts with 'solventf', not 'solvent ' — preserved
        assert 'solventfilename   "custom_thf"' in block_content

    def test_cosmors_solvent_name_preserved_without_solvent_id(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """When solvent_id is NOT set, 'solvent "name"' in custom_solvent is preserved.

        Option 2 workflow: user omits solvent_id and instead puts
        'solvent "name"' directly in custom_solvent → route gets bare
        'COSMORS' keyword, and the %cosmors block contains the solvent name.
        """
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        # Explicitly clear solvent_id so the base project YAML's solvent_id
        # (if any) does not carry over — solvent name goes into block only
        settings.solvent_id = None
        settings.custom_solvent = (
            'solvent           "THF"\n' "temp              298.15\n"
        )

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_no_id",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_cosmors_no_id.inp")
        content = open(orca_file).read()

        # Bare COSMORS in route (no parentheses)
        assert "COSMORS" in content
        assert "COSMORS(" not in content
        # 'solvent "THF"' must appear in %cosmors block (no filtering)
        assert "%cosmors" in content
        assert 'solvent           "THF"' in content
        assert "temp" in content

    # ------------------------------------------------------------------
    # Tests for solventfilename (-sf CLI option)
    # ------------------------------------------------------------------

    def test_solventfilename_written_in_cosmors_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """solventfilename setting writes solventfilename line in %cosmors block."""
        import os

        # create a dummy .cosmorsxyz file so the writer can copy it
        sf_path = os.path.join(str(tmpdir), "mysolvent.cosmorsxyz")
        open(sf_path, "w").close()

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solventfilename = sf_path

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_sf",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        out_dir = str(tmpdir.join("out"))
        orca_writer.write(target_directory=out_dir)
        orca_file = os.path.join(out_dir, "orca_cosmors_sf.inp")
        content = open(orca_file).read()

        # %cosmors block must be written
        assert "%cosmors" in content
        # solventfilename line (without .cosmorsxyz extension)
        assert 'solventfilename "mysolvent"' in content
        # Route uses bare COSMORS (no solvent_id given)
        assert "COSMORS" in content

    def test_solventfilename_with_solvent_id_in_route(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """solventfilename + solvent_id: route has COSMORS(id), block has solventfilename."""
        import os

        sf_path = os.path.join(str(tmpdir), "water.cosmorsxyz")
        open(sf_path, "w").close()

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solvent_id = "water"
        settings.solventfilename = sf_path

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_sf_id",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        out_dir = str(tmpdir.join("out2"))
        orca_writer.write(target_directory=out_dir)
        orca_file = os.path.join(out_dir, "orca_cosmors_sf_id.inp")
        content = open(orca_file).read()

        # Route uses COSMORS(water)
        assert "COSMORS(water)" in content
        # %cosmors block with solventfilename
        assert "%cosmors" in content
        assert 'solventfilename "water"' in content

    def test_solventfilename_file_copied_to_target_directory(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """The .cosmorsxyz file is copied to the target directory during write."""
        import os

        sf_path = os.path.join(str(tmpdir), "custom.cosmorsxyz")
        open(sf_path, "w").close()

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solventfilename = sf_path

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_sf_copy",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        out_dir = str(tmpdir.join("out3"))
        orca_writer.write(target_directory=out_dir)

        # The .cosmorsxyz file must be copied alongside the .inp file
        copied = os.path.join(out_dir, "custom.cosmorsxyz")
        assert os.path.isfile(copied)

    def test_solventfilename_derived_cosmorsxyz_missing_skips_copy(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """solventfilename points to a non-.cosmorsxyz file, but the
        model isn't cosmors, so _write_solvent_block never creates the
        derived .cosmorsxyz file -- the copy step must silently skip
        it instead of erroring."""
        import os

        sf_path = os.path.join(str(tmpdir), "water.xyz")
        open(sf_path, "w").close()

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cpcm"
        settings.solventfilename = sf_path

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_sf_missing_derived",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        out_dir = str(tmpdir.join("out5"))
        orca_writer.write(target_directory=out_dir)

        assert os.path.isfile(
            os.path.join(out_dir, "orca_sf_missing_derived.inp")
        )
        assert not os.path.isfile(os.path.join(out_dir, "water.cosmorsxyz"))

    def test_solventfilename_same_source_and_dest_skips_copy(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Writing directly into the same directory that already holds
        the .cosmorsxyz file must not try to copy it onto itself."""
        import os

        sf_path = os.path.join(str(tmpdir), "same_dir.cosmorsxyz")
        open(sf_path, "w").close()

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solventfilename = sf_path

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_sf_same_dir",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        # Write into the same directory that already contains sf_path.
        orca_writer.write(target_directory=str(tmpdir))

        assert os.path.isfile(
            os.path.join(str(tmpdir), "orca_sf_same_dir.inp")
        )
        # Still exactly one copy of the file -- no self-copy occurred.
        assert os.path.isfile(sf_path)

    def test_solventfilename_not_written_for_non_cosmors_model(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """solventfilename is silently ignored for non-cosmors models (cpcm, smd)."""
        import os

        sf_path = os.path.join(str(tmpdir), "water.cosmorsxyz")
        open(sf_path, "w").close()

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cpcm"
        settings.solvent_id = "water"
        settings.solventfilename = sf_path

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_sf_cpcm",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        out_dir = str(tmpdir.join("out4"))
        orca_writer.write(target_directory=out_dir)
        orca_file = os.path.join(out_dir, "orca_sf_cpcm.inp")
        content = open(orca_file).read()

        # No %cpcm block should be written (CPCM with named solvent, no extras)
        assert "%cpcm" not in content
        # solventfilename line must NOT appear for cpcm model
        assert "solventfilename" not in content

    def test_solventfilename_orca_out_converted_to_cosmorsxyz(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
        orca_he_output_freq,
    ):
        """Non-.cosmorsxyz solvent file (ORCA .out) is converted to .cosmorsxyz."""
        import os
        import shutil

        # Copy the ORCA output file to a temp location so it has a clean basename
        sf_path = os.path.join(str(tmpdir), "He_solvent.out")
        shutil.copy2(orca_he_output_freq, sf_path)

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solventfilename = sf_path

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_out",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        out_dir = str(tmpdir.join("out_conv"))
        orca_writer.write(target_directory=out_dir)
        orca_file = os.path.join(out_dir, "orca_cosmors_out.inp")
        content = open(orca_file).read()

        # %cosmors block must be written with solventfilename (without extension)
        assert "%cosmors" in content
        assert 'solventfilename "He_solvent"' in content

        # The converted .cosmorsxyz file must be present in the target directory
        copied = os.path.join(out_dir, "He_solvent.cosmorsxyz")
        assert os.path.isfile(
            copied
        ), f"Expected converted .cosmorsxyz file at {copied}"
        # The original .out file must NOT be copied (only the .cosmorsxyz)
        assert not os.path.isfile(os.path.join(out_dir, "He_solvent.out"))

    def test_solventfilename_xyz_converted_to_cosmorsxyz(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Non-.cosmorsxyz solvent file (.xyz) is converted to .cosmorsxyz."""
        import os
        import shutil

        # Copy the xyz file to a temp location so it has a clean basename
        sf_path = os.path.join(str(tmpdir), "mysolvent.xyz")
        shutil.copy2(single_molecule_xyz_file, sf_path)

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solventfilename = sf_path

        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_cosmors_xyz",
            jobrunner=orca_jobrunner_no_scratch,
        )
        orca_writer = ORCAInputWriter(job=job)
        out_dir = str(tmpdir.join("out_xyz"))
        orca_writer.write(target_directory=out_dir)
        orca_file = os.path.join(out_dir, "orca_cosmors_xyz.inp")
        content = open(orca_file).read()

        # %cosmors block must be written with solventfilename (without extension)
        assert "%cosmors" in content
        assert 'solventfilename "mysolvent"' in content

        # The converted .cosmorsxyz file must be present in the target directory
        copied = os.path.join(out_dir, "mysolvent.cosmorsxyz")
        assert os.path.isfile(
            copied
        ), f"Expected converted .cosmorsxyz file at {copied}"

    # ------------------------------------------------------------------
    # Tests for route deduplication and improved solvent detection
    # ------------------------------------------------------------------

    def test_no_duplicate_solvent_keyword_when_already_in_route(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Route deduplication guard: COSMORS(water) must appear exactly once."""
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.solvent_model = "cosmors"
        settings.solvent_id = "water"
        # Simulate a misconfigured defgrid that accidentally contains the
        # solvent keyword (triggering the duplication scenario).
        settings.defgrid = "COSMORS"

        # The route_string must NOT produce "COSMORS COSMORS(water)".
        route = settings.route_string
        # Exactly one COSMORS token (the bare one from defgrid)
        assert route.upper().count("COSMORS") == 1
        # No parenthesized COSMORS(water) appended on top
        assert (
            "COSMORS(water)" not in route
            or route.upper().count("COSMORS") == 1
        )

    def test_solvent_model_detected_from_route_line_cosmors(self):
        """ORCAFileMixin.solvent_model detects cosmors from route line."""
        import os
        import tempfile

        from chemsmart.io.orca.input import ORCAInput

        inp_content = "!  m062x def2-tzvp COSMORS(water) defgrid2\n* xyz 0 1\nO 0 0 0\n*\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".inp", delete=False
        ) as f:
            f.write(inp_content)
            fname = f.name
        try:
            orca_inp = ORCAInput(fname)
            assert orca_inp.solvent_model == "cosmors"
            assert orca_inp.solvent_id == "water"
        finally:
            os.unlink(fname)

    def test_solvent_model_detected_from_route_line_smd(self):
        """ORCAFileMixin.solvent_model detects smd from route line."""
        import os
        import tempfile

        from chemsmart.io.orca.input import ORCAInput

        inp_content = (
            "!  m062x def2-tzvp SMD(cyclohexane)\n* xyz 0 1\nO 0 0 0\n*\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".inp", delete=False
        ) as f:
            f.write(inp_content)
            fname = f.name
        try:
            orca_inp = ORCAInput(fname)
            assert orca_inp.solvent_model == "smd"
            assert orca_inp.solvent_id == "cyclohexane"
        finally:
            os.unlink(fname)

    def test_solvent_id_not_from_solventfilename_line(self):
        """solvent_id must not be extracted from a solventfilename line."""
        import os
        import tempfile

        from chemsmart.io.orca.input import ORCAInput

        inp_content = (
            "!  m062x def2-tzvp COSMORS(water)\n"
            "%cosmors\n"
            '  solventfilename "water"\n'
            "  temp 298.15\n"
            "end\n"
            "* xyz 0 1\nO 0 0 0\n*\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".inp", delete=False
        ) as f:
            f.write(inp_content)
            fname = f.name
        try:
            orca_inp = ORCAInput(fname)
            # solventfilename line must NOT be mistaken for solvent_id
            assert (
                orca_inp.solvent_id == "water"
            )  # from route line, not from solventfilename
        finally:
            os.unlink(fname)


class TestScfBlock:
    def test_neither_convergence_nor_maxiter_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_scf_block(f)
        assert f.getvalue() == ""

    def test_maxiter_only_defaults_convergence_absent(self):
        writer = _writer(ORCAJobSettings(scf_maxiter=500))
        f = io.StringIO()
        writer._write_scf_block(f)
        assert f.getvalue() == "%scf\n  maxiter 500\nend\n"

    def test_convergence_only_uses_default_maxiter(self):
        writer = _writer(ORCAJobSettings(scf_convergence="tight"))
        f = io.StringIO()
        writer._write_scf_block(f)
        assert (
            f.getvalue() == "%scf\n  maxiter 200\n  convergence tight\nend\n"
        )

    def test_convergence_with_scf_suffix_is_stripped(self):
        writer = _writer(ORCAJobSettings(scf_convergence="TightSCF"))
        f = io.StringIO()
        writer._write_scf_convergence(f)
        assert f.getvalue() == "  convergence tight\n"

    def test_invalid_convergence_raises(self):
        """See BUGS_FOUND.md: the f-string building the ValueError
        message references self.scf_convergence (an attribute of the
        writer, which doesn't exist) instead of
        self.settings.scf_convergence, so this crashes with
        AttributeError before the intended ValueError is ever raised."""
        writer = _writer(ORCAJobSettings(scf_convergence="bogusSCF"))
        with pytest.raises(AttributeError, match="scf_convergence"):
            writer._write_scf_convergence(io.StringIO())


class TestWriteSelfAndRouteQmmm:
    def test_write_self_writes_input_string_verbatim(self):
        writer = _writer(ORCAJobSettings(input_string="! B3LYP def2-SVP\n"))
        f = io.StringIO()
        writer._write_self(f)
        assert f.getvalue() == "! B3LYP def2-SVP\n"

    def test_write_without_target_directory_uses_job_folder(self, tmp_path):
        writer = _writer(
            ORCAJobSettings(
                functional="b3lyp", basis="def2-svp", charge=0, multiplicity=1
            )
        )
        writer.job.folder = str(tmp_path)
        writer._write()
        assert (tmp_path / "test_job.inp").is_file()

    def test_write_dispatches_to_write_self_for_input_string(self, tmp_path):
        writer = _writer(ORCAJobSettings(input_string="! raw input\n"))
        writer.job.folder = str(tmp_path)
        writer._write()
        assert (tmp_path / "test_job.inp").read_text() == "! raw input\n"

    def test_route_section_uses_qmmm_route_string_for_qmmm_settings(self):
        settings = ORCAQMMMJobSettings(jobtype="QMMM", parent_jobtype="opt")
        writer = _writer(settings, molecule=MagicMock(is_monoatomic=False))
        f = io.StringIO()
        writer._write_route_section(f)
        assert f.getvalue() == settings.qmmm_route_string + "\n"


class TestSolventBlockFiltering:
    def test_blank_line_preserved_when_filtering_solvent_name(self):
        """When solvent-name filtering is active (cosmors + solvent_id),
        blank lines in custom_solvent must be preserved, not treated as
        a 'solvent "name"' line to filter."""
        settings = ORCAJobSettings(
            solvent_model="cosmors",
            solvent_id="water",
            custom_solvent=(
                'solvent           "water"\n' "\n" "temp              298.15\n"
            ),
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_solvent_block(f)
        out = f.getvalue()
        assert 'solvent           "water"' not in out
        assert "temp              298.15" in out
        # The blank line is preserved (not filtered) as its own line.
        assert out == "%cosmors\n  \n  temp              298.15\nend\n"


class TestMdciBlock:
    def test_no_cutoff_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_mdci_block(f)
        assert f.getvalue() == ""

    def test_invalid_cutoff_raises(self):
        writer = _writer(ORCAJobSettings(mdci_cutoff="ultra"))
        with pytest.raises(AssertionError):
            writer._write_mdci_block(io.StringIO())

    @pytest.mark.parametrize(
        "cutoff,expected_pairs",
        [
            ("loose", "TCutPairs 1e-3"),
            ("normal", "TCutPairs 1e-4"),
            ("tight", "TCutPairs 1e-5"),
        ],
    )
    def test_cutoff_levels(self, cutoff, expected_pairs):
        writer = _writer(ORCAJobSettings(mdci_cutoff=cutoff))
        f = io.StringIO()
        writer._write_mdci_block(f)
        out = f.getvalue()
        assert out.startswith("%mdci\n")
        assert expected_pairs in out
        assert out.endswith("end\n")

    def test_invalid_density_raises(self):
        writer = _writer(
            ORCAJobSettings(mdci_cutoff="normal", mdci_density="bogus")
        )
        with pytest.raises(AssertionError):
            writer._write_mdci_block(io.StringIO())

    @pytest.mark.parametrize(
        "density,expected",
        [
            ("none", "Density None"),
            ("unrelaxed", "Density Unrelaxed"),
            ("relaxed", "Density Relaxed"),
        ],
    )
    def test_density_levels(self, density, expected):
        writer = _writer(
            ORCAJobSettings(mdci_cutoff="normal", mdci_density=density)
        )
        f = io.StringIO()
        writer._write_mdci_block(f)
        assert expected in f.getvalue()


class TestElpropBlock:
    def test_neither_dipole_nor_quadrupole_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_elprop_block(f)
        assert f.getvalue() == ""

    def test_dipole_true_quadrupole_false(self):
        writer = _writer(ORCAJobSettings(dipole=True))
        f = io.StringIO()
        writer._write_elprop_block(f)
        out = f.getvalue()
        assert "Dipole True" in out
        assert "Quadrupole False" in out

    def test_quadrupole_true_dipole_false(self):
        writer = _writer(ORCAJobSettings(quadrupole=True))
        f = io.StringIO()
        writer._write_elprop_block(f)
        out = f.getvalue()
        assert "Dipole False" in out
        assert "Quadrupole True" in out


class TestQmmmBlockDispatch:
    def test_non_qmmm_settings_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_qmmm_block(f)
        assert f.getvalue() == ""

    def test_qmmm_settings_delegates_to_settings_qmmm_block(self):
        settings = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_qmmm_block(f)
        assert f.getvalue() == f"{settings.qmmm_block}\n"


class TestModredBlock:
    def test_list_format_writes_constraints(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_modred_if_list(f, [[1, 2], [3, 4, 5]])
        out = f.getvalue()
        assert out.startswith("  Constraints\n")
        assert "{B 0 1 C}" in out
        assert "{A 2 3 4 C}" in out
        assert out.rstrip().endswith("end")

    def test_dict_format_scan_only(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        modred = {
            "coords": [[1, 2]],
            "num_steps": [10],
            "dist_start": [1.5],
            "dist_end": [3.5],
        }
        writer._write_modred_if_dict(f, modred)
        out = f.getvalue()
        assert "Constraints" not in out
        assert "  Scan\n" in out
        assert "B 0 1 = 1.5, 3.5, 10" in out
        assert "bond distance" in out

    def test_dict_format_with_constrained_coordinates(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        modred = {
            "constrained_coordinates": [[5, 6]],
            "coords": [[1, 2, 3]],
            "num_steps": [5],
            "dist_start": [90.0],
            "dist_end": [120.0],
        }
        writer._write_modred_if_dict(f, modred)
        out = f.getvalue()
        assert "Constraints" in out
        assert "{B 4 5 C}" in out
        assert "A 0 1 2 = 90.0, 120.0, 5" in out
        assert "angle" in out

    def test_dihedral_scan_variable_name(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        modred = {
            "coords": [[1, 2, 3, 4]],
            "num_steps": [8],
            "dist_start": [0.0],
            "dist_end": [180.0],
        }
        writer._write_modred_if_dict(f, modred)
        out = f.getvalue()
        assert "dihedral" in out

    def test_write_modred_dispatches_by_type(self):
        writer = _writer(ORCAJobSettings())
        f_list = io.StringIO()
        writer._write_modred(f_list, [1, 2])
        assert "Constraints" in f_list.getvalue()

        f_dict = io.StringIO()
        writer._write_modred(
            f_dict,
            {
                "coords": [[1, 2]],
                "num_steps": [4],
                "dist_start": [1.0],
                "dist_end": [2.0],
            },
        )
        assert "Scan" in f_dict.getvalue()

    def test_modred_block_wraps_in_geom_section(self):
        writer = _writer(ORCAJobSettings(modred=[1, 2]))
        f = io.StringIO()
        writer._write_modred_block(f)
        out = f.getvalue()
        assert out.startswith("%geom\n")
        assert out.rstrip().endswith("end")

    def test_no_modred_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_modred_block(f)
        assert f.getvalue() == ""

    def test_modred_of_unsupported_type_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_modred(f, modred="not a list or dict")
        assert f.getvalue() == ""


class TestHessianBlockForTS:
    def test_minimal_ts_settings(self):
        settings = ORCATSJobSettings()
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        out = f.getvalue()
        assert out.startswith("%geom\n")
        assert "Calc_Hess True" in out
        assert f"NumHess {settings.numhess}" in out
        assert f"Recalc_Hess {settings.recalc_hess}" in out
        assert out.rstrip().endswith("end")

    def test_inhess_requires_existing_file(self, tmp_path):
        hess_file = tmp_path / "guess.hess"
        hess_file.write_text("fake hessian")
        settings = ORCATSJobSettings(
            inhess=True, inhess_filename=str(hess_file)
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        out = f.getvalue()
        assert "InHess Read" in out
        assert f'InHessName "{hess_file}"' in out

    def test_inhess_without_filename_raises(self):
        settings = ORCATSJobSettings(inhess=True)
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="No Hessian file"):
            writer._write_hessian_block_for_ts(io.StringIO())

    def test_inhess_missing_file_raises(self):
        settings = ORCATSJobSettings(
            inhess=True, inhess_filename="/no/such/file.hess"
        )
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="is not found"):
            writer._write_hessian_block_for_ts(io.StringIO())

    def test_hybrid_hess_converts_to_0_indexed(self):
        settings = ORCATSJobSettings(
            hybrid_hess=True, hybrid_hess_atoms=[1, 2, 5]
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        assert "Hybrid_Hess {0 1 4} end" in f.getvalue()

    def test_hybrid_hess_without_atoms_raises(self):
        settings = ORCATSJobSettings(hybrid_hess=True)
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="No atoms"):
            writer._write_hessian_block_for_ts(io.StringIO())

    def test_trust_radius_negative_is_fixed(self):
        settings = ORCATSJobSettings(trust_radius=-0.3)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        out = f.getvalue()
        assert "Trust -0.3" in out
        assert "fixed trust radius" in out

    def test_trust_radius_positive_is_update(self):
        settings = ORCATSJobSettings(trust_radius=0.3)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        out = f.getvalue()
        assert "Trust 0.3" in out
        assert "trust radius update" in out

    def test_trust_radius_zero_writes_neither_comment(self):
        settings = ORCATSJobSettings(trust_radius=0)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        out = f.getvalue()
        assert "Trust 0" in out
        assert "fixed trust radius" not in out
        assert "trust radius update" not in out

    def test_scants_requires_modred(self):
        settings = ORCATSJobSettings(tssearch_type="scants")
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="No modred"):
            writer._write_hessian_block_for_ts(io.StringIO())

    def test_scants_writes_scan_block(self):
        settings = ORCATSJobSettings(
            tssearch_type="scants",
            scants_modred={
                "coords": [[1, 2]],
                "num_steps": [4],
                "dist_start": [1.0],
                "dist_end": [2.0],
            },
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        assert "Scan" in f.getvalue()

    def test_full_scan_writes_flag(self):
        settings = ORCATSJobSettings(full_scan=True)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_hessian_block_for_ts(f)
        assert "fullScan True" in f.getvalue()


class TestHessianBlockDispatch:
    def test_non_ts_settings_no_hessian_block(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_hessian_block(f)
        assert f.getvalue() == ""

    def test_ts_settings_dispatches(self):
        writer = _writer(ORCATSJobSettings())
        f = io.StringIO()
        writer._write_hessian_block(f)
        assert f.getvalue().startswith("%geom\n")


class TestIrcBlockForIrc:
    def test_no_irc_specific_options_writes_nothing(self):
        settings = ORCAIRCJobSettings()
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        assert f.getvalue() == ""

    def test_maxiter_option_written(self):
        settings = ORCAIRCJobSettings(maxiter=20)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        out = f.getvalue()
        assert out.startswith("%irc\n")
        assert "maxiter 20" in out
        assert out.rstrip().endswith("end")

    def test_inithess_read_requires_hess_filename(self, tmp_path):
        hess_file = tmp_path / "start.hess"
        hess_file.write_text("fake hessian")
        settings = ORCAIRCJobSettings(
            inithess="read", hess_filename=str(hess_file)
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        out = f.getvalue()
        assert "inithess read" in out
        assert f'Hess_Filename "{hess_file}"' in out

    def test_inithess_read_without_hess_filename_raises(self):
        settings = ORCAIRCJobSettings(inithess="read")
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="No Hessian file"):
            writer._write_irc_block_for_irc(io.StringIO())

    def test_inithess_calc_anfreq_no_filename_needed(self):
        settings = ORCAIRCJobSettings(inithess="calc_anfreq")
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        assert "inithess calc_anfreq" in f.getvalue()

    def test_monitor_internals_requires_internal_modred(self):
        settings = ORCAIRCJobSettings(monitor_internals=True)
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="internal modred"):
            writer._write_irc_block_for_irc(io.StringIO())

    def test_monitor_internals_writes_coordinates(self):
        settings = ORCAIRCJobSettings(
            monitor_internals=True, internal_modred=[[1, 2], [5, 6]]
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        out = f.getvalue()
        assert "True" in out
        assert "{ B 0 1 }" in out
        assert "{ B 4 5 }" in out

    def test_monitor_internals_false_writes_nothing_for_it(self):
        settings = ORCAIRCJobSettings(monitor_internals=False)
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        assert f.getvalue() == "%irc\nend\n"

    @pytest.mark.parametrize(
        "key,expected_text",
        [
            ("adapt_scale_displ", "Adapt_Scale_Displ True"),
            ("sd_parabolicfit", "SD_ParabolicFit True"),
            ("interpolate_only", "Interpolate_only True"),
            ("do_sd_corr", "Do_SD_Corr True"),
            ("sd_corr_parabolicfit", "SD_Corr_ParabolicFit True"),
        ],
    )
    def test_boolean_flag_options(self, key, expected_text):
        settings = ORCAIRCJobSettings(**{key: True})
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        assert expected_text in f.getvalue()

    @pytest.mark.parametrize(
        "key",
        [
            "adapt_scale_displ",
            "sd_parabolicfit",
            "interpolate_only",
            "do_sd_corr",
            "sd_corr_parabolicfit",
        ],
    )
    def test_boolean_flag_options_false_omitted(self, key):
        settings = ORCAIRCJobSettings(**{key: False})
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_irc_block_for_irc(f)
        # False is not None, so the key is "set", but the specific
        # branches only write when the value is True.
        assert f.getvalue() == "%irc\nend\n"


class TestIrcBlockDispatch:
    def test_non_irc_settings_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_irc_block(f)
        assert f.getvalue() == ""

    def test_irc_settings_dispatches(self):
        writer = _writer(ORCAIRCJobSettings(maxiter=5))
        f = io.StringIO()
        writer._write_irc_block(f)
        assert "maxiter 5" in f.getvalue()


class TestNebBlock:
    def test_missing_nimages_raises(self):
        settings = ORCANEBJobSettings(ending_xyzfile="product.xyz")
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="number of images"):
            _ = writer.neb_block

    def test_missing_geometry_raises(self):
        settings = ORCANEBJobSettings(nimages=8)
        writer = _writer(settings)
        with pytest.raises(AssertionError, match="valid input geometry"):
            _ = writer.neb_block

    def test_ending_xyzfile_basic_block(self):
        settings = ORCANEBJobSettings(
            nimages=8, ending_xyzfile="/some/path/product.xyz"
        )
        writer = _writer(settings)
        block = writer.neb_block
        assert block.startswith("%NEB")
        assert 'NEB_END_XYZFILE "product.xyz"' in block
        assert "NImages 8" in block
        assert "PREOPT_ENDS False" in block
        assert block.endswith("end")

    def test_preopt_ends_true(self):
        settings = ORCANEBJobSettings(
            nimages=8,
            ending_xyzfile="product.xyz",
            preopt_ends=True,
        )
        writer = _writer(settings)
        assert "PREOPT_ENDS True" in writer.neb_block

    def test_intermediate_xyzfile_included(self):
        settings = ORCANEBJobSettings(
            nimages=8,
            ending_xyzfile="product.xyz",
            intermediate_xyzfile="/some/ts_guess.xyz",
        )
        writer = _writer(settings)
        assert 'NEB_TS_XYZFILE "ts_guess.xyz"' in writer.neb_block

    def test_restart_file_takes_precedence_and_skips_preopt(self):
        settings = ORCANEBJobSettings(
            nimages=8,
            ending_xyzfile="product.xyz",
            restarting_xyzfile="/some/restart.allxyz",
        )
        writer = _writer(settings)
        block = writer.neb_block
        assert 'Restart_ALLXYZFile "restart.allxyz"' in block
        assert "NEB_END_XYZFILE" not in block
        assert "PREOPT_ENDS" not in block

    def test_write_neb_block_for_neb_matches_property(self):
        settings = ORCANEBJobSettings(nimages=8, ending_xyzfile="product.xyz")
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_neb_block_for_neb(f)
        written = f.getvalue().rstrip("\n")
        assert written == writer.neb_block

    def test_write_neb_block_for_neb_with_restart_file(self):
        """The file-writing method (unlike the neb_block property test
        above, which uses ending_xyzfile) also needs direct coverage
        of its own restarting_xyzfile branch."""
        settings = ORCANEBJobSettings(
            nimages=8,
            ending_xyzfile="product.xyz",
            restarting_xyzfile="/some/restart.allxyz",
        )
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_neb_block_for_neb(f)
        written = f.getvalue()
        assert 'Restart_ALLXYZFile "restart.allxyz"' in written
        assert "NEB_END_XYZFILE" not in written
        assert "PREOPT_ENDS" not in written


class TestNebBlockDispatch:
    def test_non_neb_settings_writes_nothing(self):
        writer = _writer(ORCAJobSettings())
        f = io.StringIO()
        writer._write_neb_block(f)
        assert f.getvalue() == ""

    def test_neb_settings_dispatches(self):
        settings = ORCANEBJobSettings(nimages=4, ending_xyzfile="product.xyz")
        writer = _writer(settings)
        f = io.StringIO()
        writer._write_neb_block(f)
        assert f.getvalue().startswith("%NEB\n")


class TestConstrainedAtoms:
    def test_no_frozen_atoms_writes_nothing(self):
        molecule = MagicMock(frozen_atoms=None)
        writer = _writer(ORCAJobSettings(), molecule=molecule)
        f = io.StringIO()
        writer._write_constrained_atoms(f)
        assert f.getvalue() == ""

    def test_frozen_atoms_written_0_indexed(self):
        molecule = MagicMock(frozen_atoms=[-1, 0, -1])
        writer = _writer(ORCAJobSettings(), molecule=molecule)
        f = io.StringIO()
        writer._write_constrained_atoms(f)
        out = f.getvalue()
        assert out.startswith("%geom\n")
        assert "{ C 0 C }" in out
        assert "{ C 2 C }" in out
        assert "{ C 1 C }" not in out
        assert out.rstrip().endswith("end")

    def test_invert_constraints_flag(self):
        molecule = MagicMock(frozen_atoms=[-1])
        writer = _writer(
            ORCAJobSettings(invert_constraints=True), molecule=molecule
        )
        f = io.StringIO()
        writer._write_constrained_atoms(f)
        assert "InvertConstraints True" in f.getvalue()


class TestChargeAndMultiplicity:
    def test_direct_charge_and_multiplicity(self):
        writer = _writer(ORCAJobSettings(charge=1, multiplicity=2))
        f = io.StringIO()
        writer._write_charge_and_multiplicity(f)
        assert f.getvalue() == "* xyz 1 2\n"

    def test_missing_both_raises(self):
        writer = _writer(ORCAJobSettings())
        with pytest.raises(AssertionError, match="Charge and multiplicity"):
            writer._write_charge_and_multiplicity(io.StringIO())

    def test_qmmm_without_charge_high_raises(self):
        # The QMMM coordinate line ("* xyz ...") must reflect the
        # high-level (QM) region's own charge/multiplicity, not the
        # whole-system total -- using total here could silently write
        # a chemically wrong charge for the QM subsystem. So when only
        # charge_total/mult_total are supplied (no charge_high/mult_high),
        # the writer must raise rather than fall back to total, even
        # though ORCAQMMMJobSettings.__init__ separately resolves
        # settings.charge/settings.multiplicity from total for other
        # purposes. See tests/test_ORCAIO.py::TestORCAQMMMJobSettings::
        # test_qmmm_xyz_line_does_not_fall_back_to_total_charge for the
        # sibling assertion of this same contract.
        settings = ORCAQMMMJobSettings(
            jobtype="QMMM",
            low_level_method="ff.prms",
            charge_total=0,
            mult_total=1,
        )
        writer = _writer(settings)
        f = io.StringIO()
        with pytest.raises(AssertionError, match="Charge and multiplicity"):
            writer._write_charge_and_multiplicity(f)
