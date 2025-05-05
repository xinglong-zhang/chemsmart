import os

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.mol.movie import PyMOLMovieJob
from chemsmart.jobs.runner import JobRunner


class PyMOLIRCMovieJob(PyMOLMovieJob):
    TYPE = "pymol_ircmovie"

    def __init__(
        self,
        molecules,
        label,
        jobrunner=None,
        **kwargs,
    ):
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
        Create a PyMOLIRCMovieJob instance from reactant and product files.
        """
        if not reactant_file and not product_file and not all_file:
            raise ValueError(
                "No reactant or product file or full irc file provided."
            )
        elif reactant_file and product_file and all_file:
            raise ValueError(
                "You can only provide reactant and product files or full irc file."
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
                "You can only provide reactant and product files or full irc file."
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
