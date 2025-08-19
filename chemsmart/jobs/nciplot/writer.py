import logging
import os.path

from chemsmart.jobs.writer import InputWriter

logger = logging.getLogger(__name__)


class NCIPLOTInputWriter(InputWriter):
    """Class that writes NCIPLOT input files for a job."""

    def write(self, **kwargs):
        self._write(**kwargs)

    def _write(self, target_directory=None):
        if target_directory is not None:
            if not os.path.exists(target_directory):
                os.makedirs(target_directory)
            folder = target_directory
        else:
            folder = self.job.folder
        job_inputfile = os.path.join(folder, f"{self.job.label}.nci")
        logger.debug(f"Writing NCIPLOT input file: {job_inputfile}")
        f = open(job_inputfile, "w")
        self._write_all(f)
        logger.info(f"Finished writing NCIPLOT input file: {job_inputfile}")
        f.close()

    def _write_all(self, f):
        self._write_filenames(f)
        self._write_rthres(f)
        self._write_ligand(f)
        self._write_radius(f)
        self._write_intermolecular(f)
        self._write_increments(f)
        self._write_fragments(f)
        self._write_cutoffs(f)
        self._write_cutplot(f)
        self._write_dgrid(f)
        self._write_integrate(f)
        self._write_ranges(f)
        self._write_grid_quality(f)

    def _write_filenames(self, f):
        logger.debug("Writing NCIPLOT input files.")
        if self.job.filenames is None:
            # this is the case when the job is created from PubChem
            logger.debug("No filenames provided for NCIPLOT job.")
            f.write("1\n")
            f.write(f"{self.job.label}.xyz\n")
        else:
            number_of_files = len(self.job.filenames)
            if number_of_files == 0:
                raise ValueError(
                    "No filenames provided for NCIPLOT job. Please provide at least one file."
                )
            else:
                logger.debug(f"Number of files: {number_of_files}")
                f.write(f"{number_of_files}\n")
                logger.debug(f"Filenames: {self.job.filenames}")
                for file in self.job.filenames:
                    logger.debug(f"Writing filename: {file}")
                    if not os.path.exists(file):
                        raise FileNotFoundError(
                            f"File {os.path.abspath(file)} does not exist. Please check the file path."
                        )
                    f.write(f"{file}\n")

    def _write_rthres(self, f):
        """Write the rthres section for the input file."""
        rthres = self.settings.rthres
        if rthres is not None:
            logger.debug("Writing rthres section.")
            f.write(f"RTHRES {rthres}\n")

    def _write_ligand(self, f):
        """Write the ligand section for the input file."""
        ligand_file_number = self.settings.ligand_file_number
        ligand_radius = self.settings.ligand_radius
        # check that if one is not None, the other is not None
        if ligand_file_number is not None and ligand_radius is not None:
            logger.debug("Writing ligand section.")
            f.write(f"LIGAND {ligand_file_number} {ligand_radius}\n")
        elif ligand_file_number is not None or ligand_radius is not None:
            raise ValueError(
                "Both ligand_file_number and ligand_radius must be provided or both must be None."
            )
        else:
            logger.debug("No ligand section written, both values are None.")

    def _write_radius(self, f):
        """Write the radius section for the input file."""
        radius_positions = self.settings.radius_positions
        radius_r = self.settings.radius_r
        # check that if one is not None, the other is not None
        if radius_positions is not None and radius_r is not None:
            logger.debug("Writing radius section.")
            radius_line = "RADIUS "
            radius_positions = radius_positions.replace("(", "")
            radius_positions = radius_positions.replace(")", "")
            logger.debug(f"radius_positions: {radius_positions}")
            coords = radius_positions.split(",")
            if len(coords) != 3:
                raise ValueError(
                    "Expected exactly 3 coordinates in 'x,y,z' format"
                )
            # Convert to floats to validate numeric values
            for c in coords:
                try:
                    float(c)
                    radius_line += f"{c} "
                except ValueError:
                    raise ValueError(f"Invalid coordinate value: {c}")
            radius_line += f"{radius_r}\n"
            f.write(radius_line)
        elif radius_positions is not None or radius_r is not None:
            raise ValueError(
                "Both radius_positions and radius_r must be provided or both must be None."
            )
        else:
            logger.debug("No radius section written, both values are None.")

    def _write_intermolecular(self, f):
        """Write the intermolecular section for the input file.

        Args:
            f: File object to write to.
        """
        intercut1 = self.settings.intercut1
        intercut2 = self.settings.intercut2

        # check if either intercut1 or intercut2 is provided
        if intercut1 is not None or intercut2 is not None:
            logger.debug("Writing intermolecular section.")
            f.write("INTERMOLECULAR\n")

            # Validate and convert intercut values to floats
            intercut1_val = float(intercut1) if intercut1 is not None else 0.95
            intercut2_val = float(intercut2) if intercut2 is not None else 0.75

            # ensure values are positive
            if intercut1_val <= 0 or intercut2_val <= 0:
                raise ValueError("INTERCUT values must be positive")

            f.write(f"INTERCUT {intercut1_val} {intercut2_val}\n")
            logger.debug(f"Wrote INTERCUT {intercut1_val} {intercut2_val}")
        else:
            logger.debug(
                "No intermolecular section written, both intercut values are None."
            )

    def _write_increments(self, f):
        """Write the increments section for the input file."""
        increments = self.settings.increments
        if increments is not None:
            logger.debug("Writing increments section.")
            increments_line = "INCREMENTS "
            increments = increments.replace("(", "")
            increments = increments.replace(")", "")
            incr = increments.split(",")

            for i in incr:
                try:
                    float(i)
                    increments_line += f"{i} "
                except ValueError:
                    raise ValueError(f"Invalid increment value: {i}")
            increments_line = increments_line.strip() + "\n"
            f.write(increments_line)
        else:
            logger.debug("No increments section written, increments is None.")

    def _write_fragments(self, f):
        """Write the fragments section for the input file.
        # FROM NCIPLOT README: https://github.com/aoterodelaroza/nciplot/blob/master/README
        # By default, the fragments to which INTERMOLECULAR applies are the
        # molecules defined in input. This behavior can be modified using the
        # **FRAGMENT** keyword:
        #
        #     FRAGMENT
        #       file1.i at1.i at2.i at3.i ...
        #       file2.i at4.i ...
        #     END | ENDFRAGMENT
        #
        # Several FRAGMENT environments can appear in the same input, each
        # defining a different fragment. The same atom should appear only in one
        # of the defined fragments. If the combined list of all the atoms inside
        # a FRAGMENT environment is not exhaustive, then the remaining atoms are
        # not considered. Inside a FRAGMENT environment, atoms can be defined by
        # choosing the molecule (ifile1.i) and then the list of atoms within
        # that molecule (at1.i at2.i etc.). If the FRAGMENT keyword is used,
        # then the INTERMOLECULAR keyword is automatically activated. Note that
        # if the EXC keyword is used, then the criterion applies also to
        # that cube. That is, if the considered point is intramolecular (more
        # than 95% of the density coming from only one fragment), then the value
        # of the exc is set to a very high value to remove the point from the
        # isosurface representation.
        #
        # In the case when the are atoms not assigned to any fragment, then the
        # second RHOCUT2 parameter (g.r) applies. A point is written to the dat
        # and cube files if:
        #
        # $\rho_i < f.r * (\sum \rho_i)$
        #
        # where i runs over fragments and:
        #
        # $(\sum \rho_i) > g.r * \rho_{total}$
        #
        # This way, only the points close to the selected fragments are
        # considered for plotting.
        ifile atom1, atom2, ..., atomn.
        Defining fragments from the .xyz files. ifile is the file label given
        by the input order and atom n is the atomic label in ifile"""
        fragments = self.settings.fragments
        if fragments is not None:
            logger.debug("Writing fragments section.")
            f.write("FRAGMENTS\n")
            # Ensure fragments is a dictionary
            if isinstance(fragments, dict):
                for key, value in fragments.items():
                    if isinstance(value, (list, tuple)):
                        # Convert each value to a string of space-separated atoms
                        atom_string = " ".join(map(str, value))
                        f.write(f"  {key} {atom_string}\n")
                    else:
                        raise ValueError(
                            f"Fragment {key} must be a list or tuple of atoms."
                        )
                f.write("END\n")
            else:
                raise ValueError("Fragments must be a dictionary.")
        else:
            logger.debug("No fragments section written, fragments is None.")

    def _write_cutoffs(self, f):
        """Write the cutoffs section for the input file."""
        cutoff_density_dat = self.settings.cutoff_density_dat
        cutoff_rdg_dat = self.settings.cutoff_rdg_dat

        # check if either cutoff_density_dat or cutoff_rdg_dat is provided
        if cutoff_density_dat is not None or cutoff_rdg_dat is not None:
            logger.debug("Writing cutoffs section.")

            # Validate and convert intercut values to floats
            cutoff_density_dat = (
                float(cutoff_density_dat)
                if cutoff_density_dat is not None
                else 0.5
            )  # default value
            cutoff_rdg_dat = (
                float(cutoff_rdg_dat) if cutoff_rdg_dat is not None else 1.0
            )  # default value

            # ensure values are positive
            if cutoff_density_dat <= 0 or cutoff_rdg_dat <= 0:
                raise ValueError("CUTOFFS values must be positive")

            f.write(f"CUTOFFS {cutoff_density_dat} {cutoff_rdg_dat}\n")
            logger.debug(
                f"Wrote CUTOFFS {cutoff_density_dat} {cutoff_rdg_dat}"
            )
        else:
            logger.debug(
                "No CUTOFFS section written, both CUTOFFS values are None."
            )

    def _write_cutplot(self, f):
        """Write the CUTPLOT section for the input file.

        Args:
            f: File object to write to.
            (promolecular or SCF).
        """
        cutoff_density_cube = self.settings.cutoff_density_cube
        cutoff_rdg_cube = self.settings.cutoff_rdg_cube

        # Determine default values based on input file type
        density = ""
        if self.job.filenames is None:
            logger.debug(
                "No filenames provided for NCIPLOT job."
                "Job likely generated from PubChem structure.\n"
                "Promolecular NCIPLOT job will be created."
            )
            density = "promolecular"
        else:
            logger.debug(
                f"Filenames provided: {self.job.filenames} for NCIPLOT job."
                "Determining density type based on file extension."
            )
            if not isinstance(self.job.filenames, (list, tuple)):
                raise TypeError(
                    f"Expected self.job.filenames to be a list or tuple, got {type(self.job.filenames).__name__}"
                )
            if len(self.job.filenames) == 0:
                raise ValueError(
                    "No filenames provided for NCIPLOT job. Please provide at least one file."
                )
            else:
                logger.debug(f"Number of files: {len(self.job.filenames)}")
                first_filename = self.job.filenames[0]
                if first_filename.endswith(".xyz"):
                    density = "promolecular"
                elif first_filename.endswith(
                    ".wfn"
                ) or first_filename.endswith(".wfx"):
                    density = "SCF"

        if density == "promolecular":
            default_density = 0.07  # Promolecular default for r1
            default_rdg = 0.3  # Promolecular default for r2
        else:
            default_density = 0.05  # SCF default for r1
            default_rdg = 0.5  # SCF default for r2

        # Check if either cutoff_density_cube or cutoff_rdg_cube is provided
        if cutoff_density_cube is not None or cutoff_rdg_cube is not None:
            logger.debug("Writing CUTPLOT section.")

            # Validate and convert cutoff values to floats
            cutoff_density_cube = (
                float(cutoff_density_cube)
                if cutoff_density_cube is not None
                else default_density
            )
            cutoff_rdg_cube = (
                float(cutoff_rdg_cube)
                if cutoff_rdg_cube is not None
                else default_rdg
            )

            # Ensure values are positive
            if cutoff_density_cube <= 0 or cutoff_rdg_cube <= 0:
                raise ValueError("CUTPLOT values must be positive")

            f.write(f"CUTPLOT {cutoff_density_cube} {cutoff_rdg_cube}\n")
            logger.debug(
                f"Wrote CUTPLOT {cutoff_density_cube} {cutoff_rdg_cube}"
            )
        else:
            logger.debug(
                "No CUTPLOT section written, both CUTPLOT values are None."
            )

    def _write_dgrid(self, f):
        """Write the DGRID section for the input file."""
        dgrid = self.settings.dgrid
        if dgrid:
            logger.debug("Writing DGRID section.")
            f.write("DGRID\n")
        else:
            logger.debug("No DGRID section written, dgrid is None.")

    def _write_integrate(self, f):
        """Write the INTEGRATE section for the input file."""
        integrate = self.settings.integrate
        if integrate:
            logger.debug("Writing INTEGRATE section.")
            f.write("INTEGRATE\n")
        else:
            logger.debug("No INTEGRATE section written, integrate is None.")

    def _write_ranges(self, f):
        """Write the RANGES section for the input file."""
        ranges = self.settings.ranges
        if ranges is not None:
            number_of_ranges = len(ranges)
            logger.debug("Writing RANGES section.")
            f.write(f"RANGE {number_of_ranges}\n")
            for r in ranges:
                if isinstance(r, (list, tuple)) and len(r) == 2:
                    # Ensure both values are floats
                    try:
                        r1, r2 = float(r[0]), float(r[1])
                        f.write(f"{r1} {r2}\n")
                        logger.debug(f"Wrote range: {r1} {r2}")
                    except ValueError as e:
                        raise ValueError(
                            f"Invalid range values: {r}. Error: {e}"
                        )
                else:
                    raise ValueError(
                        f"Each range must be a list or tuple of two numeric values. Invalid range: {r}"
                    )
        else:
            logger.debug("No RANGES section written, both values are None.")

    def _write_grid_quality(self, f):
        """Write the GRID QUALITY section for the input file."""
        grid_quality = self.settings.grid_quality
        if grid_quality is not None:
            logger.debug("Writing GRID QUALITY section.")
            f.write(f"{grid_quality.upper()}\n")
        else:
            logger.debug("No GRID QUALITY is given, default grids are used.")
