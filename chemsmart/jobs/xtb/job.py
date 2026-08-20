"""ChemSmart job objects for the supported xTB execution surface."""

import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.job import Job
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class XTBJob(Job):
    """Base class for xTB jobs compiled from validated ChemSmart settings."""

    PROGRAM = "xtb"
    _TYPE_TO_JOBTYPE = {
        "xtbsp": "sp",
        "xtbopt": "opt",
        "xtbhess": "hess",
    }

    @staticmethod
    def _has_explicit_gpu_request(jobrunner):
        """Distinguish a GPU request from a server's GPU inventory."""

        if jobrunner is None:
            return False
        requested = getattr(jobrunner, "num_gpus", None)
        marker = getattr(jobrunner, "xtb_gpu_request_explicit", None)
        if marker is not None:
            return bool(marker) and requested != 0
        server = getattr(jobrunner, "server", None)
        inherited = getattr(server, "num_gpus", None)
        if inherited is not None and requested == inherited:
            return False
        return (
            isinstance(requested, (int, float))
            and not isinstance(requested, bool)
            and requested > 0
        )

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        source_filename=None,
        project_reference=None,
        project_source_file=None,
        source_index=None,
        **kwargs,
    ):
        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of {Molecule}, got "
                f"{type(molecule).__name__}."
            )
        if not isinstance(settings, XTBJobSettings):
            raise ValueError(
                f"Settings must be instance of {XTBJobSettings}, got "
                f"{type(settings).__name__}."
            )
        if self._has_explicit_gpu_request(jobrunner):
            raise ValueError(
                "xTB is CPU-only in ChemSmart; an explicit GPU resource "
                "request is not allowed."
            )

        if label is None:
            label = molecule.get_chemical_formula(empirical=True)
        expected_jobtype = self._TYPE_TO_JOBTYPE.get(self.TYPE)
        settings.validate(expected_jobtype=expected_jobtype)
        settings.validate_for_molecule(molecule)

        self.label = label
        self._submission_execution_cwd = os.path.abspath(os.getcwd())
        super().__init__(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
        self._execution_root = os.path.abspath(self.folder)
        # Settings are the canonical execution state. Keep the job-owned
        # geometry copy aligned so downstream evidence cannot expose stale
        # charge/multiplicity metadata from the source molecule.
        self.molecule = molecule.copy()
        self.molecule.charge = settings.charge
        self.molecule.multiplicity = settings.multiplicity
        self.settings = settings.copy()
        self.source_filename = (
            None
            if source_filename is None
            else os.path.abspath(os.fspath(source_filename))
        )
        self.project_reference = project_reference
        self.project_source_file = (
            None
            if project_source_file is None
            else os.path.abspath(os.fspath(project_source_file))
        )
        self.source_index = source_index
        from chemsmart.jobs.xtb.validation import (
            bind_xtb_declared_artifacts,
        )

        self.declared_provenance_binding = bind_xtb_declared_artifacts(self)

    @classmethod
    def settings_class(cls) -> Type[XTBJobSettings]:
        return XTBJobSettings

    @property
    def xyzfile(self):
        return os.path.join(self.folder, f"{self.label}.xyz")

    @property
    def inputfile(self):
        return self.xyzfile

    @property
    def outputfile(self):
        return os.path.join(self.folder, f"{self.label}.out")

    @property
    def errfile(self):
        return os.path.join(self.folder, f"{self.label}.err")

    @property
    def environment_receiptfile(self):
        return os.path.join(
            self.folder, f"{self.label}.xtb-environment-receipt.json"
        )

    @property
    def result_receiptfile(self):
        return os.path.join(
            self.folder, f"{self.label}.xtb-result-receipt.json"
        )

    @property
    def preview_receiptfile(self):
        return os.path.join(self.folder, "xtb-preview-receipt-v1.json")

    @property
    def execution_root(self):
        """Stable parent under which private per-run workspaces are created."""

        return self._execution_root

    def bind_workspace(self, workspace):
        """Bind all path properties to one freshly created private workspace."""

        self.folder = os.path.abspath(os.fspath(workspace))

    @property
    def submission_execution_cwd(self):
        """Directory from which a reconstructed ``run`` must be invoked."""

        return self._submission_execution_cwd

    @staticmethod
    def _replace_cli_option(args, names, value, *, insert_before):
        """Replace one Click option pair in an already reconstructed argv."""

        rewritten = list(args)
        index = 0
        while index < len(rewritten):
            if rewritten[index] in names:
                del rewritten[index : index + 2]
                continue
            index += 1
        if value is None:
            return rewritten
        try:
            insert_at = (
                len(rewritten) - 1 - rewritten[::-1].index(insert_before)
            )
        except ValueError as exc:
            raise ValueError(
                f"Cannot bind xTB submission option before {insert_before!r}."
            ) from exc
        rewritten[insert_at:insert_at] = [names[0], str(value)]
        return rewritten

    def rewrite_submission_cli_args(self, cli_args):
        """Bind a scheduler-side reparse to this exact xTB job and source."""

        args = list(cli_args)
        leaf = self._TYPE_TO_JOBTYPE[self.TYPE]
        args = self._replace_cli_option(
            args,
            ("--append-label", "-a"),
            None,
            insert_before=leaf,
        )
        args = self._replace_cli_option(
            args,
            ("--filename", "-f"),
            self.source_filename,
            insert_before=leaf,
        )
        if self.project_reference is not None:
            args = self._replace_cli_option(
                args,
                ("--project", "-p"),
                self.project_reference,
                insert_before=leaf,
            )
        args = self._replace_cli_option(
            args,
            ("--index", "-i"),
            self.source_index,
            insert_before=leaf,
        )
        return self._replace_cli_option(
            args,
            ("--label", "-l"),
            self.label,
            insert_before=leaf,
        )

    def _determine_folder(self):
        folder = os.path.abspath(os.getcwd())
        if self.label:
            folder = os.path.join(folder, self.label)
            os.makedirs(folder, exist_ok=True)
        return folder

    def _backup_files(self, **kwargs):
        folder = self._create_backup_folder_name()
        self.backup_file(self.xyzfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)
        self.backup_file(self.errfile, folder=folder, **kwargs)
        self.backup_file(self.environment_receiptfile, folder=folder, **kwargs)
        self.backup_file(self.result_receiptfile, folder=folder, **kwargs)

    def _job_is_complete(self):
        """Require a hash-bound deterministic result receipt."""

        from chemsmart.jobs.xtb.validation import (
            load_validated_result_receipt,
        )

        self.validation_receipt = load_validated_result_receipt(self)
        return self.validation_receipt is not None

    def _output(self):
        if not os.path.exists(self.outputfile):
            logger.debug(f"xTB output file not found: {self.outputfile}")
            return None
        try:
            from chemsmart.io.xtb.output import XTBOutput

            return XTBOutput(folder=self.folder)
        except Exception as exc:
            logger.error(
                f"Error reading xTB output folder {self.folder}: {exc}"
            )
            return None

    def _run(self, **kwargs):
        logger.info(f"Running xTB job {self} with runner {self.jobrunner}")
        self.jobrunner.run(self, **kwargs)

    @classmethod
    def from_filename(
        cls,
        filename,
        settings=None,
        index="-1",
        label=None,
        jobrunner=None,
        **kwargs,
    ):
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        molecule = molecules[string2index_1based(index)]
        return cls(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
