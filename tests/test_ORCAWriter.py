import os
import re
import shutil
from filecmp import cmp

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.orca import (
    ORCAModredJob,
    ORCAOptJob,
    ORCAScanJob,
    ORCASinglePointJob,
    ORCATSJob,
)
from chemsmart.jobs.orca.neb import ORCANEBJob
from chemsmart.jobs.orca.writer import ORCAInputWriter
from chemsmart.settings.orca import ORCAJobSettings, ORCAProjectSettings


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
