"""PySCF job runners.

PySCF runs as a **subprocess** even when it is importable in ChemSmart's own
environment. Four reasons, none of which is isolation for its own sake:

1. ``OMP_NUM_THREADS`` and friends must be set *before* an interpreter starts
   to bind libcint and BLAS thread pools. In-process, ChemSmart's numpy has
   already initialised OpenMP and ``lib.num_threads()`` cannot fix BLAS.
2. Peak RSS (measured: 620 MB for an H2O Hessian, GBs for real systems) is
   returned to the OS on process exit instead of fragmenting ChemSmart's heap
   for the rest of a batch.
3. A libcint/libxc segfault or an OOM kill takes down only the child.
4. GPU4PySCF needs cupy + CUDA and can never share ChemSmart's environment.
   Selecting a different interpreter per program is the only mechanism that
   permits it, and it must be the same code path as CPU.
"""

import json
import logging
import math
import os
import secrets
import shlex
import subprocess
import uuid
from contextlib import suppress
from dataclasses import asdict, is_dataclass
from functools import lru_cache
from glob import glob
from numbers import Real
from shutil import copy, move

from chemsmart.io.pyscf.output import (
    PYSCF_SOURCE_ARTIFACT_INFO_KEY,
    PySCFArtifactBindingError,
    pyscf_source_artifact_binding,
    read_pyscf_h5,
)
from chemsmart.jobs.pyscf.environment import (
    canonical_sha256,
    environment_blockers,
    probe_compute_environment,
    sha256_file,
    write_json_receipt,
)
from chemsmart.jobs.pyscf.validation import (
    FREQUENCY_VALIDATION_SCHEMA_VERSION,
    RESULT_VALIDATION_SCHEMA_VERSION,
    preflight,
    validate_pyscf_result,
    verify_provenance,
)
from chemsmart.jobs.pyscf.writer import (
    PySCFScriptWriter,
    applied_pyscf_spec,
    write_pyscf_h5,
)
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import PySCFExecutable
from chemsmart.utils.process_observation import (
    launch_failure_observation,
    observe_process,
)

logger = logging.getLogger(__name__)

INPUT_RECEIPT_SCHEMA_VERSION = "chemsmart.pyscf-input.v1"
RUN_RECEIPT_SCHEMA_VERSION = "chemsmart.pyscf-run.v1"

# A fake run is a preview, so only facts that require importing the exact
# compute environment may be deferred. Scientific contradictions remain hard
# failures even when no engine is launched.
_PREVIEW_DEFERRED_RULES = frozenset(
    {
        "pyscf.solvent.database_unavailable",
        "pyscf.solver.uncallable",
        "pyscf.dispersion.dependency_missing",
        "pyscf.gpu.dependency_missing",
        "pyscf.gpu.cuda_unavailable",
        "pyscf.gpu.cutensor_incompatible",
        "pyscf.gpu.basis_angular_momentum",
        "pyscf.gpu.aux_basis_angular_momentum",
        "pyscf.gpu.functional_unverified",
        "pyscf.td.preview_only_capability",
    }
)


class PySCFPreflightError(RuntimeError):
    """Raised before launch when deterministic PySCF preflight is red."""

    def __init__(self, findings, *, receipt_path=None):
        self.findings = tuple(findings)
        self.receipt_path = receipt_path
        super().__init__(
            "PySCF preflight failed: "
            + ", ".join(
                str(finding.get("rule_id", "unknown"))
                for finding in self.findings
            )
        )


class PySCFResultValidationError(RuntimeError):
    """Raised after launch when exit or immutable result evidence is red."""

    def __init__(self, findings, *, receipt_path=None):
        self.findings = tuple(findings)
        self.receipt_path = receipt_path
        super().__init__(
            "PySCF result validation failed: "
            + ", ".join(
                str(finding.get("rule_id", "unknown"))
                for finding in self.findings
            )
        )


def _run_receipt_state(
    *, fake, findings, engine_complete, result_validation_state
):
    """Separate engine completion from scientific result validation."""

    if findings:
        return "failed"
    if fake:
        return "previewed"
    if result_validation_state == "validated":
        return "validated"
    if engine_complete and result_validation_state == "unclassified":
        # A direct Hessian command can prove an intact, completed engine
        # artifact without proving whether that geometry is a minimum or a
        # transition state.  Preserve the useful execution outcome without
        # fabricating an expected imaginary-mode count.
        return "engine_complete"
    return "failed"


def _json_safe(value):
    """Return a JSON-safe public representation without arbitrary reprs."""
    if is_dataclass(value):
        return _json_safe(asdict(value))
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        items = [_json_safe(item) for item in value]
        return sorted(
            items,
            key=lambda item: json.dumps(
                item, sort_keys=True, separators=(",", ":")
            ),
        )
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if hasattr(value, "tolist"):
        return _json_safe(value.tolist())
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return type(value).__name__


def _finalize_receipt(path, payload):
    body = _json_safe(dict(payload))
    body.pop("receipt_sha256", None)
    body["receipt_sha256"] = canonical_sha256(body)
    write_json_receipt(path, body)
    return body


def _receipt_file_sha256(path):
    """Return a receipt's verified embedded digest, or None if substituted."""
    payload = _load_bound_receipt(path)
    return payload.get("receipt_sha256") if payload else None


def _load_bound_receipt(path):
    """Load a receipt only when its embedded digest matches exact bytes."""
    try:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
        observed = payload.get("receipt_sha256")
        body = dict(payload)
        body.pop("receipt_sha256", None)
        if observed and observed == canonical_sha256(body):
            return payload
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        pass
    return None


class PySCFJobRunner(JobRunner):
    """Run a PySCF job by executing its generated driver script.

    Attributes:
        PROGRAM (str): Program identifier ('pyscf').
        JOBTYPES (list): Job types handled by this runner.
        FAKE (bool): False; the fake variant is a separate subclass.
        SCRATCH (bool): False by default -- a PySCF job writes three small
            files and gains nothing from a scratch round trip, unlike ORCA's
            large temporary integral files.
    """

    PROGRAM = "pyscf"
    JOBTYPES = [
        "pyscf_sp",
        "pyscf_opt",
        "pyscf_hess",
        "pyscf_td",
        "pyscfjob",
    ]
    FAKE = False
    SCRATCH = False
    NODE_TIMEOUT_SECONDS = 600
    PROCESS_SAMPLE_INTERVAL_SECONDS = 0.1

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )
        logger.debug(f"PySCF jobrunner num cores: {self.num_cores}")
        logger.debug(f"PySCF jobrunner num gpus: {self.num_gpus}")
        logger.debug(f"PySCF jobrunner mem gb: {self.mem_gb}")
        logger.debug(f"PySCF jobrunner scratch: {self.scratch}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """Return the PySCF executable configuration for this server."""
        try:
            return PySCFExecutable.from_servername(servername=self.server.name)
        except FileNotFoundError as e:
            logger.error(
                f"No server file {self.server} is found: {e}\n"
                f"Available servers are: {PySCFExecutable.available_servers}"
            )
            raise

    def _prerun(self, job):
        if not getattr(self, "_run_id", None):
            self._begin_run_identity()
        self._assign_variables(job)
        self._set_receipt_paths(job)
        job.settings.validate()
        self._validate_resources(job)
        try:
            self._input_artifact_binding = (
                self._current_input_artifact_binding(job)
            )
        except PySCFArtifactBindingError as exc:
            finding = _json_safe(exc.as_finding())
            self._quarantine_stale_targets(
                job,
                protected_paths=self._declared_input_artifact_paths(job),
            )
            receipt = _finalize_receipt(
                self.environment_receiptfile,
                {
                    "schema_version": "chemsmart.pyscf-environment.v1",
                    "run_id": self._run_id,
                    "run_nonce": self._run_nonce,
                    "status": "blocked_input_artifact",
                    "execution_ready": False,
                    "preflight_findings": [finding],
                    "quarantined_targets": self._quarantined_targets,
                },
            )
            self._environment_receipt = receipt
            raise PySCFPreflightError(
                [finding], receipt_path=self.environment_receiptfile
            ) from exc
        self._quarantine_stale_targets(
            job,
            protected_paths=(
                (self._input_artifact_binding["path"],)
                if self._input_artifact_binding
                else ()
            ),
        )

        if self.FAKE:
            findings = [
                _json_safe(item)
                for item in preflight(job.settings, job.molecule, {})
            ]
            hard_findings = [
                finding
                for finding in findings
                if finding.get("rule_id") not in _PREVIEW_DEFERRED_RULES
            ]
            receipt = _finalize_receipt(
                self.environment_receiptfile,
                {
                    "schema_version": "chemsmart.pyscf-environment.v1",
                    "run_id": self._run_id,
                    "run_nonce": self._run_nonce,
                    "status": "deferred_preview",
                    "execution_ready": False,
                    "input_artifact_kind": (
                        self._input_artifact_binding["kind"]
                        if self._input_artifact_binding
                        else None
                    ),
                    "input_artifact_path": (
                        self._input_artifact_binding["path"]
                        if self._input_artifact_binding
                        else None
                    ),
                    "input_artifact_sha256": (
                        self._input_artifact_binding["sha256"]
                        if self._input_artifact_binding
                        else None
                    ),
                    "reason": (
                        "fake preview does not import the compute interpreter"
                    ),
                    "preflight_findings": findings,
                    "deferred_rule_ids": sorted(
                        {
                            finding["rule_id"]
                            for finding in findings
                            if finding["rule_id"] in _PREVIEW_DEFERRED_RULES
                        }
                    ),
                },
            )
            self._environment_receipt = receipt
            if hard_findings:
                raise PySCFPreflightError(
                    hard_findings, receipt_path=self.environment_receiptfile
                )
            return

        child_env = self._update_os_environ(job)
        receipt = probe_compute_environment(
            interpreter=self.executable.get_executable(),
            job=job,
            env=child_env,
            receipt_path=self.environment_receiptfile,
        )
        findings = environment_blockers(receipt, engine=job.settings.engine)
        findings.extend(
            _json_safe(item)
            for item in preflight(job.settings, job.molecule, receipt)
        )
        receipt["preflight_findings"] = findings
        receipt["execution_ready"] = not findings
        receipt["run_id"] = self._run_id
        receipt["run_nonce"] = self._run_nonce
        receipt["input_artifact_kind"] = (
            self._input_artifact_binding["kind"]
            if self._input_artifact_binding
            else None
        )
        receipt["input_artifact_path"] = (
            self._input_artifact_binding["path"]
            if self._input_artifact_binding
            else None
        )
        receipt["input_artifact_sha256"] = (
            self._input_artifact_binding["sha256"]
            if self._input_artifact_binding
            else None
        )
        receipt = _finalize_receipt(self.environment_receiptfile, receipt)
        self._environment_receipt = receipt
        if findings:
            raise PySCFPreflightError(
                findings, receipt_path=self.environment_receiptfile
            )

    def _assign_variables(self, job):
        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

        if self.executable and self.executable.local_run is not None:
            job.local = self.executable.local_run

    def _set_up_variables_in_scratch(self, job):
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        self._set_file_paths(job, scratch_job_dir)

    def _set_up_variables_in_job_directory(self, job):
        self.running_directory = job.folder
        self._set_file_paths(job, job.folder)

    def _set_file_paths(self, job, directory):
        self.job_inputfile = os.path.abspath(
            os.path.join(directory, job.label + ".py")
        )
        self.job_outputfile = os.path.abspath(
            os.path.join(directory, job.label + ".out")
        )
        self.job_resultsfile = os.path.abspath(
            os.path.join(directory, job.label + ".h5")
        )
        self.job_errfile = os.path.abspath(
            os.path.join(directory, job.label + ".err")
        )
        logger.debug(f"Running directory: {self.running_directory}")

    def _set_receipt_paths(self, job):
        # Controller receipts live in the durable job folder even when the
        # calculation itself uses scratch. A preflight failure must not strand
        # its only evidence in an ephemeral directory.
        self.environment_receiptfile = os.path.abspath(
            os.path.join(job.folder, job.label + ".environment.json")
        )
        self.input_receiptfile = os.path.abspath(
            os.path.join(job.folder, job.label + ".input.json")
        )
        self.run_receiptfile = os.path.abspath(
            os.path.join(job.folder, job.label + ".receipt.json")
        )

    def _quarantine_stale_targets(self, job, protected_paths=()):
        """Move only exact prior run targets into a recoverable quarantine."""
        if getattr(self, "_targets_quarantined", False):
            return
        protected = {
            os.path.realpath(os.path.abspath(path))
            for path in protected_paths
            if isinstance(path, str) and path
        }
        candidates = (
            job.inputfile,
            job.outputfile,
            job.resultsfile,
            job.errfile,
            self.environment_receiptfile,
            self.input_receiptfile,
            self.run_receiptfile,
            self.job_inputfile,
            self.job_outputfile,
            self.job_resultsfile,
            self.job_errfile,
        )
        unique = []
        seen = set()
        for value in candidates:
            path = os.path.abspath(value)
            if path not in seen:
                unique.append(path)
                seen.add(path)
        quarantine_root = os.path.join(
            job.folder,
            ".chemsmart-stale",
            job.label,
            self._run_id,
        )
        quarantined = []
        for index, source in enumerate(unique):
            if os.path.realpath(source) in protected:
                continue
            if not os.path.isfile(source):
                continue
            os.makedirs(quarantine_root, exist_ok=True)
            destination = os.path.join(
                quarantine_root,
                f"{index:02d}-{os.path.basename(source)}",
            )
            source_sha256 = sha256_file(source)
            move(source, destination)
            quarantined.append(
                {
                    "source": source,
                    "source_sha256": source_sha256,
                    "quarantined_to": os.path.relpath(
                        destination, job.folder
                    ),
                }
            )
        self._quarantined_targets = quarantined
        self._targets_quarantined = True

    def _validate_resources(self, job):
        def positive_integer(name, value):
            if (
                isinstance(value, bool)
                or not isinstance(value, int)
                or value <= 0
            ):
                raise ValueError(
                    f"{name} must be a positive integer, got {value!r}"
                )

        positive_integer("num_cores", self.num_cores)
        if self.mem_gb is not None and (
            isinstance(self.mem_gb, bool)
            or not isinstance(self.mem_gb, Real)
            or not math.isfinite(float(self.mem_gb))
            or float(self.mem_gb) <= 0
        ):
            raise ValueError(
                "mem_gb must be a positive finite scalar, got "
                f"{self.mem_gb!r}"
            )
        if isinstance(self.num_gpus, bool) or not isinstance(
            self.num_gpus, int
        ):
            raise ValueError(
                f"num_gpus must be a non-negative integer, got {self.num_gpus!r}"
            )
        if self.num_gpus < 0:
            raise ValueError(
                f"num_gpus must be a non-negative integer, got {self.num_gpus!r}"
            )
        if (
            job.settings.engine == "gpu"
            and self.num_gpus <= 0
            and not self.FAKE
        ):
            raise ValueError(
                "engine='gpu' requires a positive resolved num_gpus; "
                "ChemSmart will not invent a GPU device"
            )

    def _current_input_artifact_binding(self, job):
        """Verify the source HDF5 and reject source/target aliasing."""
        binding = pyscf_source_artifact_binding(job.molecule)
        if binding is None:
            return None
        source = os.path.realpath(binding["path"])
        targets = {
            os.path.realpath(os.path.abspath(path))
            for path in (
                job.inputfile,
                job.outputfile,
                job.resultsfile,
                job.errfile,
                self.job_inputfile,
                self.job_outputfile,
                self.job_resultsfile,
                self.job_errfile,
            )
        }
        if source in targets:
            raise PySCFArtifactBindingError(
                "pyscf.input_artifact.target_collision",
                "input_artifact.path",
                "source artifact distinct from every output target",
                source,
                evidence_ref=f"file:{source}",
            )
        return binding

    @staticmethod
    def _declared_input_artifact_paths(job):
        """Return raw declared source paths without trusting the binding."""
        info = getattr(job.molecule, "info", None)
        binding = (
            info.get(PYSCF_SOURCE_ARTIFACT_INFO_KEY)
            if isinstance(info, dict)
            else None
        )
        path = binding.get("path") if isinstance(binding, dict) else None
        return (path,) if isinstance(path, str) and path else ()

    def _write_input(self, job):
        """Generate the driver script in the running directory."""
        writer = PySCFScriptWriter(job=job)
        try:
            binding = self._current_input_artifact_binding(job)
            config = writer.build_config()
        except PySCFArtifactBindingError as exc:
            finding = _json_safe(exc.as_finding())
            _finalize_receipt(
                self.input_receiptfile,
                {
                    "schema_version": INPUT_RECEIPT_SCHEMA_VERSION,
                    "run_id": self._run_id,
                    "run_nonce": self._run_nonce,
                    "state": "blocked",
                    "findings": [finding],
                },
            )
            raise PySCFPreflightError(
                [finding], receipt_path=self.input_receiptfile
            ) from exc
        if binding != getattr(self, "_input_artifact_binding", None):
            finding = {
                "rule_id": "pyscf.input_artifact.binding_changed",
                "field": "input_artifact",
                "expected": getattr(self, "_input_artifact_binding", None),
                "observed": binding,
                "evidence_ref": "molecule:source_artifact",
            }
            _finalize_receipt(
                self.input_receiptfile,
                {
                    "schema_version": INPUT_RECEIPT_SCHEMA_VERSION,
                    "run_id": self._run_id,
                    "run_nonce": self._run_nonce,
                    "state": "blocked",
                    "findings": [finding],
                },
            )
            raise PySCFPreflightError(
                [finding], receipt_path=self.input_receiptfile
            )
        writer.write(target_directory=self.running_directory, config=config)
        receipt = _finalize_receipt(
            self.input_receiptfile,
            {
                "schema_version": INPUT_RECEIPT_SCHEMA_VERSION,
                "run_id": self._run_id,
                "run_nonce": self._run_nonce,
                "state": "previewed",
                "label": job.label,
                "script_path": os.path.basename(self.job_inputfile),
                "script_sha256": sha256_file(self.job_inputfile),
                "config_sha256": canonical_sha256(config),
                "input_geometry_sha256": config[
                    "input_geometry_sha256"
                ],
                "input_artifact_kind": config.get("input_artifact_kind"),
                "input_artifact_path": (
                    binding["path"] if binding else None
                ),
                "input_artifact_sha256": config.get(
                    "input_artifact_sha256"
                ),
                "project_yaml_sha256": config.get("project_yaml_digest"),
                "requested_settings_sha256": config[
                    "requested_settings_sha256"
                ],
                "environment_receipt_sha256": self._environment_receipt[
                    "receipt_sha256"
                ],
                "compute_interpreter": (
                    None
                    if self.FAKE
                    else self._environment_receipt.get("interpreter")
                ),
                "stages": config["stages"],
                "engine": config["engine"],
                "materializations": config.get("materializations", {}),
                "quarantined_targets": self._quarantined_targets,
            },
        )
        self._input_config = config
        self._input_receipt = receipt

    def _get_command(self, job):
        """Return typed argv for ``<interpreter> <label>.py``.

        For a library backend the executable is the interpreter that owns
        pyscf, resolved from the server YAML ``PYSCF.EXEFOLDER`` when set and
        falling back to the running interpreter.
        """
        del job
        exe = (
            getattr(self, "_environment_receipt", {}).get("interpreter")
            or self.executable.get_executable()
        )
        logger.info(f"PySCF interpreter: {exe}")
        return (str(exe), str(self.job_inputfile))

    def _update_os_environ(self, job):
        """Return the child environment with thread and device limits set.

        These must be in place *before* the interpreter starts: BLAS and
        OpenMP size their pools at load time, so setting them afterwards has
        no effect. This is the substantive reason PySCF runs out-of-process.
        """
        env = super()._update_os_environ(job)
        threads = str(self.num_cores)
        for var in (
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
        ):
            env[var] = threads

        engine = getattr(job.settings, "engine", "cpu")
        if engine == "gpu":
            num_gpus = self.num_gpus
            env.setdefault(
                "CUDA_VISIBLE_DEVICES",
                ",".join(str(i) for i in range(num_gpus)),
            )
        else:
            # Make a CPU run reproducible even on a GPU node: without this a
            # stray gpu4pyscf import could quietly pick up a device.
            env["CUDA_VISIBLE_DEVICES"] = ""
        environment_receipt = getattr(self, "_environment_receipt", None)
        input_receipt = getattr(self, "_input_receipt", None)
        if environment_receipt:
            env["CHEMSMART_PYSCF_ENVIRONMENT_RECEIPT_SHA256"] = (
                environment_receipt["receipt_sha256"]
            )
        if input_receipt:
            env["CHEMSMART_PYSCF_SCRIPT_SHA256"] = input_receipt[
                "script_sha256"
            ]
            env["CHEMSMART_PYSCF_INPUT_RECEIPT_SHA256"] = input_receipt[
                "receipt_sha256"
            ]
        return env

    def run(self, job, **kwargs):
        """Run the hardened PySCF lifecycle and propagate the child status."""
        self._reset_run_state()
        self._begin_run_identity()
        self._prerun(job)
        self._write_input(job)
        command = self._get_command(job)
        env = self._update_os_environ(job)
        try:
            process = self._create_process(job, command=command, env=env)
            self._child_returncode = self._run(process, **kwargs)
        except Exception as exc:
            self._child_returncode = None
            self._process_observation = launch_failure_observation(
                timeout_seconds=self.NODE_TIMEOUT_SECONDS,
                memory_limit_mb=self._process_memory_limit_mb(),
                error_type=type(exc).__name__,
            ).as_dict()
            self._launch_failure = {
                "type": type(exc).__name__,
                "message": str(exc)[:500],
            }
        self._postrun(job)
        try:
            self._postrun_cleanup(job)
        except PySCFResultValidationError as exc:
            if isinstance(self._child_returncode, int) and (
                self._child_returncode != 0
            ):
                raise subprocess.CalledProcessError(
                    self._child_returncode, command
                ) from exc
            raise
        return self._child_returncode

    def _reset_run_state(self):
        for name in (
            "_child_returncode",
            "_environment_receipt",
            "_input_receipt",
            "_input_config",
            "_input_artifact_binding",
            "_launch_failure",
            "_process_observation",
            "_run_receipt",
            "_quarantined_targets",
            "_targets_quarantined",
            "_run_id",
            "_run_nonce",
        ):
            self.__dict__.pop(name, None)

    def _begin_run_identity(self):
        self._run_id = str(uuid.uuid4())
        self._run_nonce = secrets.token_hex(32)

    def _process_memory_limit_mb(self):
        """Return the exact configured process-memory boundary in MiB."""

        return None if self.mem_gb is None else float(self.mem_gb) * 1024.0

    def _run(self, process, **kwargs):
        """Observe the child tree under fixed time and memory boundaries."""

        del kwargs
        result = observe_process(
            process,
            timeout_seconds=self.NODE_TIMEOUT_SECONDS,
            memory_limit_mb=self._process_memory_limit_mb(),
            sample_interval_seconds=self.PROCESS_SAMPLE_INTERVAL_SECONDS,
        )
        self._process_observation = result.observation.as_dict()
        return result.observation.returncode

    def _postrun_cleanup(self, job):
        """Validate immutable result provenance before declaring success."""
        returncode = getattr(self, "_child_returncode", None)
        expected = {
            "run_id": self._run_id,
            "run_nonce": self._run_nonce,
            "script_sha256": self._input_receipt["script_sha256"],
            "input_receipt_sha256": self._input_receipt["receipt_sha256"],
            "environment_receipt_sha256": self._environment_receipt[
                "receipt_sha256"
            ],
            "input_geometry_sha256": self._input_receipt[
                "input_geometry_sha256"
            ],
            "input_artifact_kind": self._input_receipt[
                "input_artifact_kind"
            ],
            "input_artifact_sha256": self._input_receipt[
                "input_artifact_sha256"
            ],
            "requested_settings_sha256": self._input_receipt[
                "requested_settings_sha256"
            ],
            "project_yaml_digest": self._input_receipt[
                "project_yaml_sha256"
            ],
            "require_applied_settings_sha256": not self.FAKE,
            "require_engine_complete": not self.FAKE,
        }
        findings = []
        if self.FAKE:
            # A fake artifact proves compilation and serialization only.  Its
            # deliberately incomplete engine status is not a preview failure,
            # while every other provenance mismatch remains a hard finding.
            findings.extend(
                _json_safe(item)
                for item in verify_provenance(
                    job.settings,
                    job.resultsfile,
                    expected_receipt=expected,
                )
                if item.rule_id != "pyscf.provenance.incomplete_calculation"
            )
        try:
            current_input_artifact = self._current_input_artifact_binding(job)
        except PySCFArtifactBindingError as exc:
            current_input_artifact = None
            findings.append(_json_safe(exc.as_finding()))
        artifact_checks = (
            (
                "script_sha256",
                expected["script_sha256"],
                sha256_file(job.inputfile)
                if os.path.isfile(job.inputfile)
                else None,
                f"file:{os.path.basename(job.inputfile)}",
            ),
            (
                "input_receipt_sha256",
                expected["input_receipt_sha256"],
                _receipt_file_sha256(self.input_receiptfile),
                f"file:{os.path.basename(self.input_receiptfile)}",
            ),
            (
                "environment_receipt_sha256",
                expected["environment_receipt_sha256"],
                _receipt_file_sha256(self.environment_receiptfile),
                f"file:{os.path.basename(self.environment_receiptfile)}",
            ),
            (
                "input_geometry_sha256",
                expected["input_geometry_sha256"],
                job._current_geometry_sha256(),
                "job:current_geometry",
            ),
            (
                "input_artifact_kind",
                expected["input_artifact_kind"],
                (
                    current_input_artifact["kind"]
                    if current_input_artifact
                    else None
                ),
                "molecule:source_artifact",
            ),
            (
                "input_artifact_path",
                self._input_receipt["input_artifact_path"],
                (
                    current_input_artifact["path"]
                    if current_input_artifact
                    else None
                ),
                "molecule:source_artifact",
            ),
            (
                "input_artifact_sha256",
                expected["input_artifact_sha256"],
                (
                    current_input_artifact["sha256"]
                    if current_input_artifact
                    else None
                ),
                "molecule:source_artifact",
            ),
        )
        for field, expected_value, observed_value, evidence_ref in (
            artifact_checks
        ):
            if expected_value != observed_value:
                findings.append(
                    {
                        "rule_id": "pyscf.artifact.binding_mismatch",
                        "field": field,
                        "expected": expected_value,
                        "observed": observed_value,
                        "evidence_ref": evidence_ref,
                    }
                )
        launch_failure = getattr(self, "_launch_failure", None)
        process_observation = getattr(self, "_process_observation", None)
        if launch_failure:
            findings.append(
                {
                    "rule_id": "pyscf.process.launch_failed",
                    "field": "process.launch",
                    "expected": "child process started and observed",
                    "observed": launch_failure,
                    "evidence_ref": f"file:{os.path.basename(job.errfile)}",
                }
            )
        elif process_observation and process_observation.get("timed_out"):
            findings.append(
                {
                    "rule_id": "pyscf.process.timeout",
                    "field": "process.wall_seconds",
                    "expected": {
                        "maximum_seconds": process_observation.get(
                            "timeout_seconds"
                        )
                    },
                    "observed": process_observation.get("wall_seconds"),
                    "evidence_ref": "run:process_observation",
                }
            )
            if not process_observation.get("termination_confirmed"):
                findings.append(
                    {
                        "rule_id": "pyscf.process.termination_ambiguous",
                        "field": "process.termination_confirmed",
                        "expected": True,
                        "observed": process_observation.get(
                            "termination_confirmed"
                        ),
                        "evidence_ref": "run:process_observation",
                    }
                )
        elif process_observation and process_observation.get(
            "memory_limit_exceeded"
        ):
            findings.append(
                {
                    "rule_id": "pyscf.process.memory_limit_exceeded",
                    "field": "process.peak_rss_mb",
                    "expected": {
                        "maximum_mb": self._process_memory_limit_mb()
                    },
                    "observed": process_observation.get("peak_rss_mb"),
                    "evidence_ref": "run:process_observation",
                }
            )
            if not process_observation.get("termination_confirmed"):
                findings.append(
                    {
                        "rule_id": "pyscf.process.termination_ambiguous",
                        "field": "process.termination_confirmed",
                        "expected": True,
                        "observed": process_observation.get(
                            "termination_confirmed"
                        ),
                        "evidence_ref": "run:process_observation",
                    }
                )
        elif returncode != 0:
            findings.append(
                {
                    "rule_id": "pyscf.process.nonzero_exit",
                    "field": "process.returncode",
                    "expected": 0,
                    "observed": returncode,
                    "evidence_ref": f"file:{os.path.basename(job.errfile)}",
                }
            )

        result_sha256 = (
            sha256_file(job.resultsfile)
            if os.path.isfile(job.resultsfile)
            else None
        )
        applied_settings_sha256 = None
        result_spec = {}
        result_provenance = {}
        result_status = {}
        result_values = {}
        if result_sha256 is not None:
            try:
                (
                    result_spec,
                    result_provenance,
                    result_status,
                    result_values,
                ) = read_pyscf_h5(job.resultsfile)
                applied_settings_sha256 = result_spec.get(
                    "applied_settings_sha256"
                )
            except (OSError, KeyError, TypeError, ValueError):
                pass
        if not self.FAKE:
            findings.extend(
                self._runtime_environment_findings(result_provenance)
            )
        property_findings = self._property_findings(job, result_status)
        findings.extend(
            finding for finding in property_findings if finding["required"]
        )
        result_validation = {
            "schema_version": RESULT_VALIDATION_SCHEMA_VERSION,
            "state": "not_evaluated_preview",
            "jobtype": str(job.settings.jobtype).lower(),
            "findings": [],
            "frequency_validation": {
                "schema_version": FREQUENCY_VALIDATION_SCHEMA_VERSION,
                "state": (
                    "not_evaluated_preview"
                    if "hess" in job.stages
                    else "not_applicable"
                ),
                "findings": [],
            },
        }
        if not self.FAKE:
            result_validation = _json_safe(
                validate_pyscf_result(
                    job.resultsfile,
                    settings=job.settings,
                    expected_jobtype=self._input_config["jobtype"],
                    expected_charge=self._input_config["charge"],
                    expected_multiplicity=self._input_config["multiplicity"],
                    expected_symbols=self._input_config["symbols"],
                    expected_positions=self._input_config["positions"],
                    expected_receipt=expected,
                )
            )
            findings.extend(result_validation["findings"])
        frequency_validation = result_validation["frequency_validation"]

        engine_complete = bool(
            not self.FAKE
            and returncode == 0
            and result_status.get("engine_complete") is True
        )
        state = _run_receipt_state(
            fake=bool(self.FAKE),
            findings=findings,
            engine_complete=engine_complete,
            result_validation_state=result_validation.get("state"),
        )
        receipt = _finalize_receipt(
            self.run_receiptfile,
            {
                "schema_version": RUN_RECEIPT_SCHEMA_VERSION,
                "run_id": self._run_id,
                "run_nonce": self._run_nonce,
                "state": state,
                "fake": bool(self.FAKE),
                "child_returncode": returncode,
                "process_observation": process_observation,
                "engine_complete": engine_complete,
                "scientifically_validated": state == "validated",
                "scientific_validation_state": result_validation.get(
                    "state"
                ),
                "script_sha256": expected["script_sha256"],
                "input_receipt_sha256": expected[
                    "input_receipt_sha256"
                ],
                "environment_receipt_sha256": expected[
                    "environment_receipt_sha256"
                ],
                "input_geometry_sha256": expected[
                    "input_geometry_sha256"
                ],
                "input_artifact_kind": expected["input_artifact_kind"],
                "input_artifact_path": self._input_receipt[
                    "input_artifact_path"
                ],
                "input_artifact_sha256": expected[
                    "input_artifact_sha256"
                ],
                "project_yaml_sha256": expected["project_yaml_digest"],
                "requested_settings_sha256": expected[
                    "requested_settings_sha256"
                ],
                "applied_settings_sha256": applied_settings_sha256,
                "result_sha256": result_sha256,
                "findings": findings,
                "property_findings": property_findings,
                "result_validation": result_validation,
                "frequency_validation": frequency_validation,
                "quarantined_targets": self._quarantined_targets,
            },
        )
        self._run_receipt = receipt
        if findings:
            raise PySCFResultValidationError(
                findings, receipt_path=self.run_receiptfile
            )
        if self.scratch and self.delete_scratch and not self.FAKE:
            self._delete_scratch_directory()
        return receipt

    def _runtime_environment_findings(self, provenance):
        receipt = _load_bound_receipt(self.environment_receiptfile)
        if not receipt:
            return [
                {
                    "rule_id": "pyscf.environment.receipt_unreadable",
                    "field": "environment_receipt",
                    "expected": "digest-valid exact environment receipt",
                    "observed": None,
                    "evidence_ref": "environment:receipt",
                }
            ]
        findings = []

        def compare(field, expected, observed, evidence_ref):
            if expected != observed:
                findings.append(
                    {
                        "rule_id": "pyscf.environment.runtime_mismatch",
                        "field": field,
                        "expected": expected,
                        "observed": observed,
                        "evidence_ref": evidence_ref,
                    }
                )

        interpreter = receipt.get("interpreter")
        observed_executable_sha256 = (
            sha256_file(interpreter)
            if interpreter and os.path.isfile(interpreter)
            else None
        )
        compare(
            "interpreter_sha256",
            receipt.get("interpreter_sha256"),
            observed_executable_sha256,
            "environment:interpreter_sha256",
        )
        compare(
            "provenance.interpreter",
            receipt.get("interpreter_observed"),
            provenance.get("interpreter"),
            "h5:/provenance/interpreter",
        )
        dependencies = receipt.get("dependencies") or {}
        for dependency, provenance_field in (
            ("pyscf", "pyscf_version"),
            ("numpy", "numpy_version"),
            ("h5py", "h5py_version"),
        ):
            compare(
                f"provenance.{provenance_field}",
                (dependencies.get(dependency) or {}).get("version"),
                provenance.get(provenance_field),
                f"h5:/provenance/{provenance_field}",
            )
        compare(
            "provenance.libxc_version",
            receipt.get("libxc_version"),
            provenance.get("libxc_version"),
            "h5:/provenance/libxc_version",
        )
        if provenance.get("engine") == "gpu":
            runtime = provenance.get("runtime") or {}
            compare(
                "provenance.gpu4pyscf_version",
                (receipt.get("gpu4pyscf_distribution") or {}).get(
                    "version"
                ),
                provenance.get("gpu4pyscf_version"),
                "h5:/provenance/gpu4pyscf_version",
            )
            for field in (
                "device_count",
                "device_name",
                "device_uuid",
                "cuda_driver_version",
                "cuda_runtime_version",
                "cutensor_runtime",
            ):
                compare(
                    f"provenance.runtime.{field}",
                    receipt.get(field),
                    runtime.get(field),
                    f"h5:/provenance/runtime/{field}",
                )
            runtime_packages = runtime.get("packages") or {}
            for distribution_field in (
                "gpu4pyscf_distribution",
                "cupy_distribution",
                "cutensor_distribution",
            ):
                distribution = receipt.get(distribution_field) or {}
                name = distribution.get("name")
                compare(
                    f"provenance.runtime.packages.{name}",
                    distribution.get("version"),
                    runtime_packages.get(name) if name else None,
                    "h5:/provenance/runtime/packages",
                )
        return findings

    @staticmethod
    def _property_findings(job, status):
        required = set(
            getattr(job, "required_properties", ())
            or getattr(job, "kwargs", {}).get("required_properties", ())
        )
        findings = []
        properties = status.get("properties", {}) if isinstance(status, dict) else {}
        for name in sorted(properties):
            detail = properties[name]
            if not isinstance(detail, dict) or detail.get("status") == "ok":
                continue
            findings.append(
                {
                    "rule_id": "pyscf.property.unavailable",
                    "field": f"properties.{name}",
                    "expected": "available" if name in required else "optional",
                    "observed": detail,
                    "required": name in required,
                    "evidence_ref": f"h5:/status/properties/{name}",
                }
            )
        for name in sorted(required.difference(properties)):
            findings.append(
                {
                    "rule_id": "pyscf.property.missing",
                    "field": f"properties.{name}",
                    "expected": "required property status",
                    "observed": "<missing>",
                    "required": True,
                    "evidence_ref": f"h5:/status/properties/{name}",
                }
            )
        return findings

    def _create_process(self, job, command, env):
        with open(self.job_errfile, "w") as err:
            logger.info(
                f"Command executed: {shlex.join(command)}\n"
                f"PySCF log: {self.job_outputfile}\n"
                f"Results: {self.job_resultsfile}\n"
                f"Errors: {self.job_errfile}"
            )
            return subprocess.Popen(
                command,
                stdout=err,
                stderr=subprocess.STDOUT,
                env=env,
                cwd=self.running_directory,
                start_new_session=True,
            )

    def _postrun(self, job):
        """Copy artifacts back from scratch, if scratch was used."""
        if not self.scratch:
            return
        for path in glob(f"{self.running_directory}/{job.label}*"):
            try:
                copy(path, job.folder)
            except Exception as e:
                logger.error(f"Failed to copy {path} to {job.folder}: {e}")


class FakePySCFJobRunner(PySCFJobRunner):
    """Fake PySCF runner for tests and dry runs.

    Writes the real driver script and a synthetic results file **without
    importing pyscf**, so ``--fake`` exercises the full
    CLI -> settings -> writer -> output -> parse path in any environment,
    including one where PySCF is not installed at all.

    This is deliberately unlike ``FakeORCA``, which fabricates log text for a
    real regex parser to re-parse. Because the structured results file is the
    only programmatic contract here, the fake can populate it directly.
    """

    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)

    def run(self, job, **kwargs):
        self._reset_run_state()
        self._begin_run_identity()
        self._prerun(job)
        self._write_input(job)
        self._write_fake_results(job)
        self._postrun(job)
        self._child_returncode = 0
        self._postrun_cleanup(job)
        return 0

    def _write_fake_results(self, job):
        """Write a synthetic ``label.h5`` and a marker log."""
        import numpy as np

        config = self._input_config

        spec = applied_pyscf_spec(config)
        spec["fake"] = True

        status = {
            "stages": {
                stage: {
                    "state": "not_evaluated_preview",
                    "converged": None,
                }
                for stage in config["stages"]
            },
            "normal_termination": False,
            "engine_complete": False,
            "evaluation_state": "not_evaluated_preview",
            "failure": None,
            "preview_only": True,
            "properties": {
                name: {
                    "status": "unavailable",
                    "failure": {
                        "type": "not_executed",
                        "message": "fake preview performs no calculation",
                    },
                }
                for name in (
                    "forces",
                    "mulliken_charges",
                    "dipole_moment",
                    "point_group",
                )
            },
        }
        provenance = {
            "run_id": config["run_id"],
            "run_nonce": config["run_nonce"],
            "pyscf_version": "fake",
            "gpu4pyscf_version": None,
            "engine": config["engine"],
            "num_threads": config["num_threads"],
            "settings_digest": config["settings_digest"],
            "requested_settings_sha256": config[
                "requested_settings_sha256"
            ],
            "applied_settings_sha256": None,
            "input_geometry_sha256": config[
                "input_geometry_sha256"
            ],
            "input_artifact_kind": config.get("input_artifact_kind"),
            "input_artifact_sha256": config.get("input_artifact_sha256"),
            "project_yaml_digest": config.get("project_yaml_digest"),
            "script_sha256": self._input_receipt["script_sha256"],
            "input_receipt_sha256": self._input_receipt[
                "receipt_sha256"
            ],
            "environment_receipt_sha256": self._environment_receipt[
                "receipt_sha256"
            ],
            "wall_seconds": 0.0,
            "core_seconds": 0.0,
            "fake": True,
            "runtime": {"mean_field_class": None, "preview_only": True},
        }

        # Preserve the exact preview geometry for inspection without
        # fabricating energies, orbital data, atomic numbers, modes, or any
        # other numerical chemistry result.
        results = {
            "positions": np.asarray(config["positions"], dtype=float),
        }

        with open(self.job_outputfile, "w") as handle:
            handle.write(
                "PySCF version 0.0.0-fake\n"
                "PySCF path  (fake run; no calculation was performed)\n"
            )

        write_pyscf_h5(
            self.job_resultsfile,
            spec=spec,
            provenance=provenance,
            status=status,
            results=results,
        )
