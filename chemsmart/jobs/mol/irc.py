import os

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.mol.movie import PyMOLMovieJob
from chemsmart.jobs.runner import JobRunner


class PyMOLIRCMovieJob(PyMOLMovieJob):
    """
    PyMOL job for creating Intrinsic Reaction Coordinate (IRC) movies.

    Specialized movie job for visualizing reaction pathways by
    animating molecular structures along IRC trajectories, showing
    the progression from reactants to products through the transition state.

    Attributes:
        TYPE (str): Job type identifier ('pymol_ircmovie').
        molecules (list[Molecule]): Ordered IRC trajectory structures used
            to generate the animation.
        label (str): Job identifier used for file naming and outputs.
        jobrunner (JobRunner): Execution backend for running the job.
    """

    TYPE = "pymol_ircmovie"

    def __init__(
        self,
        molecules,
        label,
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialize a PyMOL IRC movie job.

        Sets up the job for creating animated visualizations of
        intrinsic reaction coordinate trajectories showing molecular
        changes along the reaction pathway.

        Args:
            molecules: List of Molecule objects representing IRC trajectory.
            label: Job identifier string.
            jobrunner: Job execution runner (default: None).
            **kwargs: Additional arguments passed to parent PyMOLMovieJob.
        """
        super().__init__(
            molecule=molecules, label=label, jobrunner=jobrunner, **kwargs
        )

    @classmethod
    def from_files(
        cls,
        reactant_file,
        product_file,
        all_file,
        label,
        jobrunner=None,
        **kwargs,
    ):
        """
        Create PyMOL IRC movie job from trajectory files.

        Creates an IRC animation job from either separate reactant/product
        files or a complete IRC trajectory file, handling molecular
        structure loading and sequence ordering.

        Args:
            reactant_file: Path to reactant trajectory file.
            product_file: Path to product trajectory file.
            all_file: Path to complete IRC trajectory file.
            label: Job identifier string. If `reactant_file`/`product_file` or
                `all_file` is provided, the label is derived from filenames
                and the passed label is ignored.
            jobrunner: Job execution runner (default: auto-created).
            **kwargs: Additional arguments for job configuration.

        Returns:
            PyMOLIRCMovieJob: Configured IRC movie job instance.

        Raises:
            ValueError: If file inputs are invalid or conflicting.
        """
        if not reactant_file and not product_file and not all_file:
            raise ValueError(
                "No reactant or product file or full irc file provided."
            )
        elif reactant_file and product_file and all_file:
            raise ValueError(
                "You can only provide reactant and product files or "
                "full irc file."
            )
        elif reactant_file and product_file:
            # find the overlap of the two file names
            reactant_file = os.path.basename(reactant_file)
            product_file = os.path.basename(product_file)
            label = os.path.commonprefix([reactant_file, product_file])

            reactant_molecules = Molecule.from_filepath(
                reactant_file, index="::-1"
            )
            product_molecules = Molecule.from_filepath(product_file, index=":")
            molecules = reactant_molecules + product_molecules
        elif all_file:
            label = os.path.splitext(os.path.basename(all_file))[0]
            all_molecules = Molecule.from_filepath(all_file, index=":")
            molecules = all_molecules
        else:
            raise ValueError(
                "You can only provide reactant and product files or "
                "full irc file."
            )

        # Create jobrunner if not provided
        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecules=molecules,
                    label=label,
                    jobrunner=None,
                    **kwargs,
                ),
                server=kwargs.get("server"),
                scratch=kwargs.get("scratch"),
                fake=kwargs.get("fake", False),
                **kwargs,
            )

        return cls(
            molecules=molecules, label=label, jobrunner=jobrunner, **kwargs
        )
