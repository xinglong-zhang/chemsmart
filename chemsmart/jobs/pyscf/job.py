"""PySCF job classes.

Unlike Gaussian and ORCA, PySCF is a Python library rather than an external
binary. ChemSmart therefore *generates* the calculation script and controls
both ends of the interface, which changes the file contract:

``label.py``   the generated, standalone-runnable driver script
``label.out``  PySCF's own log (human debugging; nothing parses it)
``label.h5``   the structured results -- the sole programmatic contract
``label.err``  stderr from the subprocess
"""

import json
import logging
import os
from typing import Type

from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.pyscf.output import (
    PySCFArtifactBindingError,
    pyscf_source_artifact_binding,
)
from chemsmart.jobs.job import Job
from chemsmart.jobs.pyscf.environment import canonical_sha256, sha256_file
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.jobs.pyscf.validation import verify_provenance
from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class PySCFJob(Job):
    """Base PySCF job.

    Attributes:
        PROGRAM (str): Program identifier ('PySCF').
        molecule (Molecule): Molecular structure used for the calculation.
        settings (PySCFJobSettings): Configuration options for the job.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    PROGRAM = "PySCF"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        super().__init__(
            molecule=molecule, label=label, jobrunner=jobrunner, **kwargs
        )

        if not isinstance(settings, PySCFJobSettings):
            raise ValueError(
                f"Settings must be instance of {PySCFJobSettings} for "
                f"{self}, but is {settings} instead!"
            )

        if not isinstance(molecule, Molecule):
            raise ValueError(
                f"Molecule must be instance of Molecule for {self}, but is "
                f"{molecule} instead!"
            )

        self.molecule = molecule.copy() if molecule is not None else None
        self.settings = settings.copy()
        # Preserve an electronic state carried by the source artifact unless
        # the project or CLI already supplied an explicit replacement.  This
        # also keeps direct programmatic job construction consistent with the
        # CLI's per-molecule resolution path.
        if self.settings.charge is None and self.molecule.charge is not None:
            self.settings.charge = self.molecule.charge
        if (
            self.settings.multiplicity is None
            and self.molecule.multiplicity is not None
        ):
            self.settings.multiplicity = self.molecule.multiplicity

        if label is None:
            label = molecule.get_chemical_formula(empirical=True)
        self.label = label

    @classmethod
    def settings_class(cls) -> Type[PySCFJobSettings]:
        return PySCFJobSettings

    @property
    def inputfile(self):
        """Path to the generated driver script."""
        return os.path.join(self.folder, self.label + ".py")

    @property
    def outputfile(self):
        """Path to PySCF's own log file."""
        return os.path.join(self.folder, self.label + ".out")

    @property
    def resultsfile(self):
        """Path to the structured results file.

        This is the only file read programmatically. ``label.out`` exists for
        humans, and is deliberately never parsed.
        """
        return os.path.join(self.folder, self.label + ".h5")

    @property
    def errfile(self):
        return os.path.join(self.folder, self.label + ".err")

    @property
    def environment_receiptfile(self):
        return os.path.join(self.folder, self.label + ".environment.json")

    @property
    def input_receiptfile(self):
        return os.path.join(self.folder, self.label + ".input.json")

    @property
    def run_receiptfile(self):
        return os.path.join(self.folder, self.label + ".receipt.json")

    def _backup_files(self, **kwargs):
        """Back up the script, log and results together.

        All three are needed to reproduce or audit a run, so backing up a
        subset would leave an unusable record.
        """
        folder = self._create_backup_folder_name()
        self.backup_file(self.inputfile, folder=folder, **kwargs)
        self.backup_file(self.outputfile, folder=folder, **kwargs)
        self.backup_file(self.resultsfile, folder=folder, **kwargs)
        self.backup_file(self.errfile, folder=folder, **kwargs)
        self.backup_file(
            self.environment_receiptfile, folder=folder, **kwargs
        )
        self.backup_file(self.input_receiptfile, folder=folder, **kwargs)
        self.backup_file(self.run_receiptfile, folder=folder, **kwargs)

    @staticmethod
    def _read_bound_receipt(path):
        with open(path, encoding="utf-8") as handle:
            receipt = json.load(handle)
        if not isinstance(receipt, dict):
            return None
        observed = receipt.get("receipt_sha256")
        body = dict(receipt)
        body.pop("receipt_sha256", None)
        if not observed or observed != canonical_sha256(body):
            return None
        return receipt

    def _current_geometry_sha256(self):
        """Digest the exact geometry and effective electronic state now held."""
        charge = self.settings.charge
        if charge is None:
            charge = self.molecule.charge
        multiplicity = self.settings.multiplicity
        if multiplicity is None:
            multiplicity = self.molecule.multiplicity
        payload = {
            "symbols": list(self.molecule.chemical_symbols),
            "positions": [
                [float(x), float(y), float(z)]
                for x, y, z in self.molecule.positions
            ],
            "unit": "Angstrom",
            "charge": int(charge),
            "multiplicity": int(multiplicity),
        }
        return canonical_sha256(payload)

    def _job_is_complete(self):
        """Return whether the engine artifact reports normal termination.

        This is deliberately independent of modern ChemSmart provenance.
        ``Job.run`` uses completion only to avoid rerunning an already
        completed calculation.  Missing or stale evidence lowers the state to
        unvalidated; it must not silently launch the chemistry again.
        """

        output = self._output()
        return output is not None and output.normal_termination is True

    def is_validated(self):
        """Return true only for a hash-bound, deterministically validated run."""
        try:
            run_receipt = self._read_bound_receipt(self.run_receiptfile)
            input_receipt = self._read_bound_receipt(self.input_receiptfile)
            environment_receipt = self._read_bound_receipt(
                self.environment_receiptfile
            )
            if not all((run_receipt, input_receipt, environment_receipt)):
                return False
            if (
                run_receipt.get("state") != "validated"
                or run_receipt.get("fake") is not False
                or run_receipt.get("child_returncode") != 0
                or run_receipt.get("findings")
            ):
                return False
            run_id = run_receipt.get("run_id")
            run_nonce = run_receipt.get("run_nonce")
            if not run_id or not run_nonce:
                return False
            if input_receipt.get("run_id") != run_id or input_receipt.get(
                "run_nonce"
            ) != run_nonce:
                return False
            if environment_receipt.get(
                "run_id"
            ) != run_id or environment_receipt.get("run_nonce") != run_nonce:
                return False
            if not os.path.isfile(self.inputfile) or not os.path.isfile(
                self.resultsfile
            ):
                return False
            if run_receipt.get("script_sha256") != sha256_file(
                self.inputfile
            ):
                return False
            if run_receipt.get("result_sha256") != sha256_file(
                self.resultsfile
            ):
                return False
            if run_receipt.get("input_receipt_sha256") != input_receipt.get(
                "receipt_sha256"
            ):
                return False
            if run_receipt.get(
                "environment_receipt_sha256"
            ) != environment_receipt.get("receipt_sha256"):
                return False
            if run_receipt.get(
                "input_geometry_sha256"
            ) != self._current_geometry_sha256():
                return False
            try:
                input_artifact = pyscf_source_artifact_binding(self.molecule)
            except PySCFArtifactBindingError:
                return False
            observed_artifact_kind = (
                input_artifact["kind"] if input_artifact else None
            )
            observed_artifact_path = (
                input_artifact["path"] if input_artifact else None
            )
            observed_artifact_sha256 = (
                input_artifact["sha256"] if input_artifact else None
            )
            if run_receipt.get(
                "input_artifact_kind"
            ) != observed_artifact_kind:
                return False
            if run_receipt.get(
                "input_artifact_path"
            ) != observed_artifact_path:
                return False
            if run_receipt.get(
                "input_artifact_sha256"
            ) != observed_artifact_sha256:
                return False

            expected = {
                "run_id": run_receipt.get("run_id"),
                "run_nonce": run_receipt.get("run_nonce"),
                "script_sha256": run_receipt.get("script_sha256"),
                "input_receipt_sha256": run_receipt.get(
                    "input_receipt_sha256"
                ),
                "environment_receipt_sha256": run_receipt.get(
                    "environment_receipt_sha256"
                ),
                "input_geometry_sha256": run_receipt.get(
                    "input_geometry_sha256"
                ),
                "requested_settings_sha256": run_receipt.get(
                    "requested_settings_sha256"
                ),
                "project_yaml_digest": run_receipt.get(
                    "project_yaml_sha256"
                ),
                "require_applied_settings_sha256": True,
            }
            if "input_artifact_kind" in run_receipt:
                expected["input_artifact_kind"] = run_receipt.get(
                    "input_artifact_kind"
                )
            if "input_artifact_sha256" in run_receipt:
                expected["input_artifact_sha256"] = run_receipt.get(
                    "input_artifact_sha256"
                )
            return not verify_provenance(
                self.settings,
                self.resultsfile,
                expected_receipt=expected,
            )
        except (OSError, TypeError, ValueError, json.JSONDecodeError):
            return False

    def _output(self):
        """Return a ``PySCFOutput``, or None when results are absent.

        Keyed on the results file rather than the log: a run that produced a
        log but no ``.h5`` did not finish, and must not read as complete.
        """
        if not os.path.exists(self.resultsfile):
            logger.debug(f"Results file not found: {self.resultsfile}")
            return None

        try:
            from chemsmart.io.pyscf.output import PySCFOutput

            return PySCFOutput(self.resultsfile)
        except (OSError, KeyError, ValueError) as e:
            logger.error(f"Error creating PySCFOutput object: {e}")
            return None

    def _run(self, **kwargs):
        logger.info(f"Running PySCFJob {self} with jobrunner {self.jobrunner}")
        self.jobrunner.run(self, **kwargs)

    @classmethod
    def from_filename(
        cls,
        filename,
        settings=None,
        index="-1",
        label=None,
        jobrunner=None,
        keywords=("charge", "multiplicity"),
        **kwargs,
    ):
        """Create a PySCF job from a molecular structure file."""
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        logger.info(f"Number of molecules read: {len(molecules)}")
        molecules = molecules[string2index_1based(index)]

        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecules,
                    settings=settings,
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
            molecule=molecules,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

    @classmethod
    def from_pubchem(
        cls, identifier, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """Create a PySCF job from a PubChem identifier."""
        molecule = Molecule.from_pubchem(identifier=identifier)

        if jobrunner is None:
            jobrunner = JobRunner.from_job(
                cls(
                    molecule=molecule,
                    settings=settings,
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
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )


class PySCFGeneralJob(PySCFJob):
    """General PySCF job.

    Exists so that specialised subclasses can create and run a plain PySCF
    job internally without recursing into their own type.
    """

    TYPE = "pyscfjob"

    def __init__(self, molecule, settings=None, label=None, **kwargs):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
