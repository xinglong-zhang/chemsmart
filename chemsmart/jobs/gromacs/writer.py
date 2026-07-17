from __future__ import annotations

"""
GROMACS input writer.

This module writes GROMACS .mdp files from ChemSmart GROMACS job settings.

The current implementation supports automatic MDP generation for:
- EM  / gmxem
- NVT / gmxnvt
- NPT / gmxnpt

Topology preparation and structure conversion are intentionally not handled
here. Those belong to the GROMACS workflow setup layer.
"""

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class GromacsInputWriter:
    """
    Write GROMACS input files for a GROMACS job.

    If the user already supplies job.mdp_file, this writer does not overwrite
    it. If job.mdp_file is missing, the writer creates a default .mdp file in
    the job folder and assigns it back to job.mdp_file.

    This keeps user-provided MDP files as an advanced override while allowing
    ChemSmart to generate reasonable default inputs automatically.
    """

    def __init__(self, job, jobrunner=None):
        self.job = job
        self.jobrunner = jobrunner

    def write(self, target_directory=None):
        """
        Write missing GROMACS input files.

        Returns:
            Path | None: The written .mdp file path, or None if no file was
            written because job.mdp_file was already provided.
        """
        return self.write_mdp(target_directory=target_directory)

    def write_mdp(self, target_directory=None):
        """
        Write an .mdp file if job.mdp_file is not already set.

        Existing user-provided MDP files are respected and never overwritten.
        """
        if self.job.mdp_file is not None:
            logger.info(
                "Using user-provided GROMACS MDP file: %s",
                self.job.mdp_file,
            )
            return None

        folder = Path(target_directory or self.job.folder)
        folder.mkdir(parents=True, exist_ok=True)

        job_type = self.job.TYPE.lower()

        if job_type == "gmxem":
            mdp_path = folder / "em.mdp"
            content = self._build_em_mdp()
        elif job_type == "gmxnvt":
            mdp_path = folder / "nvt.mdp"
            content = self._build_nvt_mdp()
        elif job_type == "gmxnpt":
            mdp_path = folder / "npt.mdp"
            content = self._build_npt_mdp()
        else:
            raise ValueError(
                f"Unsupported GROMACS job type for MDP writing: {job_type}"
            )

        logger.info("Writing GROMACS MDP file: %s", mdp_path)
        mdp_path.write_text(content, encoding="utf-8")

        self.job.mdp_file = mdp_path
        return mdp_path

    def _build_em_mdp(self):
        """
        Build a default EM .mdp file.
        """
        emtol = self._value_or_default(self.job.emtol, 1000.0)
        emstep = self._value_or_default(self.job.emstep, 0.01)
        nsteps = self._value_or_default(self.job.nsteps, 50000)

        return self._format_mdp(
            {
                "integrator": "steep",
                "emtol": emtol,
                "emstep": emstep,
                "nsteps": nsteps,
                "cutoff-scheme": "Verlet",
                "nstlist": 10,
                "rcoulomb": 1.0,
                "rvdw": 1.0,
                "coulombtype": "PME",
                "pbc": "xyz",
            }
        )

    def _build_nvt_mdp(self):
        """
        Build a default NVT .mdp file.
        """
        temperature = self._value_or_default(self.job.temperature, 300)
        timestep = self._value_or_default(self.job.timestep, 0.002)
        nsteps = self._value_or_default(self.job.nsteps, 50000)
        constraints = self._value_or_default(self.job.constraints, "h-bonds")
        constraint_algorithm = self._value_or_default(
            self.job.constraint_algorithm,
            "lincs",
        )
        thermostat = self._value_or_default(self.job.thermostat, "V-rescale")
        tau_t = self._value_or_default(self.job.tau_t, 0.1)
        tc_grps = self._value_or_default(self.job.tc_grps, "System")

        return self._format_mdp(
            {
                "integrator": "md",
                "nsteps": nsteps,
                "dt": timestep,
                "continuation": "no",
                "constraint_algorithm": constraint_algorithm,
                "constraints": constraints,
                "cutoff-scheme": "Verlet",
                "nstlist": 10,
                "rcoulomb": 1.0,
                "rvdw": 1.0,
                "coulombtype": "PME",
                "pbc": "xyz",
                "tcoupl": thermostat,
                "tc-grps": tc_grps,
                "tau_t": tau_t,
                "ref_t": temperature,
                "pcoupl": "no",
                "gen_vel": "yes",
                "gen_temp": temperature,
                "gen_seed": -1,
            }
        )

    def _build_npt_mdp(self):
        """
        Build a default NPT .mdp file.
        """
        temperature = self._value_or_default(self.job.temperature, 300)
        pressure = self._value_or_default(self.job.pressure, 1.0)
        timestep = self._value_or_default(self.job.timestep, 0.002)
        nsteps = self._value_or_default(self.job.nsteps, 50000)
        constraints = self._value_or_default(self.job.constraints, "h-bonds")
        constraint_algorithm = self._value_or_default(
            self.job.constraint_algorithm,
            "lincs",
        )
        thermostat = self._value_or_default(self.job.thermostat, "V-rescale")
        barostat = self._value_or_default(
            self.job.barostat,
            "Parrinello-Rahman",
        )
        tau_t = self._value_or_default(self.job.tau_t, 0.1)
        tc_grps = self._value_or_default(self.job.tc_grps, "System")
        tau_p = self._value_or_default(self.job.tau_p, 2.0)
        compressibility = self._value_or_default(
            self.job.compressibility,
            "4.5e-5",
        )

        return self._format_mdp(
            {
                "integrator": "md",
                "nsteps": nsteps,
                "dt": timestep,
                "continuation": "yes",
                "constraint_algorithm": constraint_algorithm,
                "constraints": constraints,
                "cutoff-scheme": "Verlet",
                "nstlist": 10,
                "rcoulomb": 1.0,
                "rvdw": 1.0,
                "coulombtype": "PME",
                "pbc": "xyz",
                "tcoupl": thermostat,
                "tc-grps": tc_grps,
                "tau_t": tau_t,
                "ref_t": temperature,
                "pcoupl": barostat,
                "pcoupltype": "isotropic",
                "tau_p": tau_p,
                "ref_p": pressure,
                "compressibility": compressibility,
                "gen_vel": "no",
            }
        )

    @staticmethod
    def _value_or_default(value, default=None):
        """
        Return default when the given value is None.
        """
        return default if value is None else value

    @staticmethod
    def _format_mdp(options):
        """
        Format a dictionary of MDP options into GROMACS .mdp text.
        """
        lines = []

        for key, value in options.items():
            if value is None:
                continue
            lines.append(f"{key:<24} = {value}")

        return "\n".join(lines) + "\n"