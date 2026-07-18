"""Plain-text run report for CHEMSMART Iterate jobs.

This module renders a Gaussian-log-style ``.out`` report describing an
iterate run: the input structures, the resolved algorithm configuration,
combination generation, per-combination execution results, and the final
delivery statistics. The report is intentionally plain text (never JSON) and
is written atomically so a partial file is never left behind.

The structured :class:`CombinationResult` replaces the old
``(label, Molecule | None)`` tuple so the report can distinguish generation
failures, timeouts, and output-write failures.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from datetime import datetime
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from chemsmart.io.molecules.structure import Molecule

# --- Execution status values ---------------------------------------------
STATUS_SUCCESS = "SUCCESS"
STATUS_FAILED = "FAILED"
STATUS_TIMED_OUT = "TIMED OUT"
STATUS_WRITE_FAILED = "WRITE FAILED"

# --- Failure stages ------------------------------------------------------
STAGE_PREPROCESSING = "PREPROCESSING"
STAGE_ALGORITHM = "ALGORITHM"
STAGE_WRITE = "WRITE"

# --- Run-level error codes (Gaussian-style; several may apply at once) ---
# The run no longer carries a single top-level status string. Each distinct
# problem contributes an error code, and a run terminates normally only when
# no code applies. Combination-level statuses (SUCCESS / FAILED / TIMED OUT /
# WRITE FAILED) above are unaffected.
ERROR_CODE_INPUT = "ITR-INPUT-001"
ERROR_CODE_EXEC = "ITR-EXEC-001"
ERROR_CODE_TIMEOUT = "ITR-TIMEOUT-001"
ERROR_CODE_WRITE = "ITR-WRITE-001"
ERROR_CODE_INTERNAL = "ITR-INTERNAL-001"
ERROR_CODE_INTERRUPTED = "ITR-INTERRUPTED-001"

ERROR_CODE_DESCRIPTIONS = {
    ERROR_CODE_INPUT: "Input file load or preprocessing error",
    ERROR_CODE_EXEC: "Combination execution failure",
    ERROR_CODE_TIMEOUT: "Combination execution timeout",
    ERROR_CODE_WRITE: "Structure output write failure",
    ERROR_CODE_INTERNAL: "Unexpected internal error",
    ERROR_CODE_INTERRUPTED: "Run interrupted by user",
}

_WIDTH = 80


@dataclass
class CombinationResult:
    """Structured outcome for a single combination.

    Attributes
    ----------
    combination_number : int
        Stable 1-based index in the original combination order.
    label : str
        Unique combination label.
    execution_status : str
        One of ``SUCCESS``, ``FAILED``, ``TIMED OUT`` or ``WRITE FAILED``.
    molecule : Molecule or None
        The generated structure (``None`` on failure/timeout).
    duration_seconds : float or None
        Wall-clock time spent on the combination.
    failure_stage : str or None
        ``PREPROCESSING`` / ``ALGORITHM`` / ``WRITE`` (input-file load
        failures are tracked as input errors, not per-combination stages).
    error_type : str or None
        Exception class name (concise; no traceback in the report).
    error_message : str or None
        Short human-readable failure reason.
    output_path : str or None
        Path the structure was written to (merged file or per-structure file).
    structure_index : int or None
        1-based index of the structure within its output file/set.
    """

    combination_number: int
    label: str
    execution_status: str
    molecule: Optional["Molecule"] = None
    duration_seconds: Optional[float] = None
    failure_stage: Optional[str] = None
    error_type: Optional[str] = None
    error_message: Optional[str] = None
    output_path: Optional[str] = None
    structure_index: Optional[int] = None


def summarize_results(results: list, input_error_count: int = 0) -> dict:
    """Compute the combination/output statistics and the run status.

    A combination that generated a structure but could not be written counts
    as *generated* for execution stats yet as *not delivered* for output
    stats, so the two views stay internally consistent:

        total     == generated + failed + timed_out
        generated == write_succeeded + write_failed

    ``input_error_count`` is the number of declared skeleton/substituent
    files that failed to load; any such failure makes the run an
    ``INPUT ERROR`` regardless of how many structures were produced.

    The returned ``error_codes`` list holds one Gaussian-style code per
    distinct problem (input / execution / timeout / write); there is no
    single top-level status. ``exit_code`` is 0 only when the run is
    completely clean, and 1 whenever any error code applies -- including a
    run that produced some structures but also had failures or timeouts.
    """
    generated = sum(
        1
        for r in results
        if r.execution_status in (STATUS_SUCCESS, STATUS_WRITE_FAILED)
    )
    write_succeeded = sum(
        1 for r in results if r.execution_status == STATUS_SUCCESS
    )
    write_failed = sum(
        1 for r in results if r.execution_status == STATUS_WRITE_FAILED
    )
    failed = sum(1 for r in results if r.execution_status == STATUS_FAILED)
    timed_out = sum(
        1 for r in results if r.execution_status == STATUS_TIMED_OUT
    )
    total = len(results)

    # Gaussian-style: collect one error code per distinct problem. There is
    # no single top-level status and no priority ordering -- any code present
    # means an error termination.
    error_codes: list[str] = []
    if input_error_count > 0:
        error_codes.append(ERROR_CODE_INPUT)
    if failed > 0:
        error_codes.append(ERROR_CODE_EXEC)
    if timed_out > 0:
        error_codes.append(ERROR_CODE_TIMEOUT)
    if write_failed > 0:
        error_codes.append(ERROR_CODE_WRITE)
    exit_code = 1 if error_codes else 0
    return {
        "total": total,
        "generated": generated,
        "write_succeeded": write_succeeded,
        "write_failed": write_failed,
        "failed": failed,
        "timed_out": timed_out,
        "input_errors": input_error_count,
        "error_codes": error_codes,
        "exit_code": exit_code,
    }


@dataclass
class IterateReport:
    """All data needed to render an iterate run report."""

    run_id: str
    chemsmart_version: str
    rdkit_version: str
    started_at: datetime
    finished_at: datetime
    duration_seconds: float
    working_directory: str
    command_line: str
    config_file: Optional[str]
    config_sha256: str

    skeleton_entries: list = field(default_factory=list)
    substituent_entries: list = field(default_factory=list)
    input_errors: list = field(default_factory=list)

    algorithm_name: str = ""
    algorithm_options: dict = field(default_factory=dict)
    combination_mode: str = ""
    nprocs: int = 1
    timeout_seconds: float = 0.0
    output_mode: str = "merged"

    loaded_skeletons: int = 0
    loaded_substituents: int = 0
    attachment_site_count: int = 0
    total_combinations: int = 0
    per_skeleton_counts: list = field(default_factory=list)

    results: list = field(default_factory=list)
    output_location: Optional[str] = None
    # Optional overrides for abnormal terminations (e.g. INTERNAL ERROR /
    # INTERRUPTED). ``extra_error_codes`` are merged with the codes derived
    # from ``results``; ``exit_code_override`` forces the exit code (e.g. 130
    # for SIGINT); ``error_message`` adds a short human-readable detail line.
    extra_error_codes: list = field(default_factory=list)
    exit_code_override: Optional[int] = None
    error_message: Optional[str] = None

    def render(self) -> str:
        """Render the full plain-text report."""
        stats = summarize_results(self.results, len(self.input_errors))
        error_codes = list(stats["error_codes"])
        for code in self.extra_error_codes:
            if code not in error_codes:
                error_codes.append(code)
        exit_code = (
            self.exit_code_override
            if self.exit_code_override is not None
            else stats["exit_code"]
        )
        lines: list[str] = []
        lines += self._render_header()
        lines += self._render_general()
        lines += self._render_inputs()
        lines += self._render_input_errors()
        lines += self._render_algorithm()
        lines += self._render_generation()
        lines += self._render_execution()
        lines += self._render_failed()
        lines += self._render_timed_out()
        lines += self._render_write_failures()
        lines += self._render_output_structures(stats)
        lines += self._render_final_summary(stats, error_codes, exit_code)
        return "\n".join(lines) + "\n"

    # --- section renderers ----------------------------------------------

    @staticmethod
    def _bar(char: str) -> str:
        return " " + char * (_WIDTH - 1)

    def _section(self, title: str) -> list:
        return ["", self._bar("-"), " " + title, self._bar("-"), ""]

    def _render_header(self) -> list:
        title = "CHEMSMART ITERATE JOB REPORT"
        return [
            self._bar("="),
            " " + title.center(_WIDTH - 1).rstrip(),
            self._bar("="),
        ]

    def _render_general(self) -> list:
        started = self.started_at.strftime("%Y-%m-%d %H:%M:%S")
        finished = self.finished_at.strftime("%Y-%m-%d %H:%M:%S")
        lines = ["", " GENERAL INFORMATION", ""]
        lines += [
            f" Run ID:                  {self.run_id}",
            f" CHEMSMART version:       {self.chemsmart_version}",
            f" RDKit version:           {self.rdkit_version}",
            "",
            f" Start time:              {started}",
            f" Finish time:             {finished}",
            f" Total elapsed time:      {_format_duration(self.duration_seconds)}",
            f" Working directory:       {self.working_directory}",
            f" Command:                 {self.command_line}",
            "",
            f" Configuration file:      {self.config_file or 'N/A'}",
            f" Configuration SHA256:    {self.config_sha256}",
        ]
        return lines

    def _render_inputs(self) -> list:
        lines = self._section("INPUT STRUCTURES")
        lines.append(f" Skeletons: {len(self.skeleton_entries)}")
        lines.append("")
        for i, entry in enumerate(self.skeleton_entries, start=1):
            label = entry.get("label") or f"skeleton{i}"
            lines += [
                f"   [{i}] Label:             {label}",
                f"       File as written:   {entry.get('file_path_raw')}",
                f"       Resolved file:     {entry.get('file_path')}",
                f"       Attachment sites:  {_format_sites(entry)}",
                f"       Skeleton indices:  "
                f"{_format_list(entry.get('skeleton_indices'))}",
                "",
            ]
        lines.append(f" Substituents: {len(self.substituent_entries)}")
        lines.append("")
        for i, entry in enumerate(self.substituent_entries, start=1):
            label = entry.get("label") or f"substituent{i}"
            link = entry.get("link_index")
            link_str = link[0] if link else "N/A"
            lines += [
                f"   [{i}] Label:             {label}",
                f"       File as written:   {entry.get('file_path_raw')}",
                f"       Resolved file:     {entry.get('file_path')}",
                f"       Link index:        {link_str}",
                f"       Groups:            "
                f"{_format_list(entry.get('groups'))}",
                "",
            ]
        return lines

    def _render_input_errors(self) -> list:
        lines = self._section("INPUT ERRORS")
        if not self.input_errors:
            lines.append(" None.")
            return lines
        for i, err in enumerate(self.input_errors, start=1):
            lines += [
                f" [{i}] Type:            {err.get('type')}",
                f"     Label:           {err.get('label')}",
                f"     File as written: {err.get('raw_path')}",
                f"     Resolved file:   {err.get('resolved_path')}",
                f"     Error type:      {err.get('error_type')}",
                f"     Reason:          {err.get('error_message')}",
                "",
            ]
        return lines

    def _render_algorithm(self) -> list:
        lines = self._section("ALGORITHM")
        lines += [
            f" Algorithm:               {self.algorithm_name}",
            f" Combination mode:        {self.combination_mode}",
            "",
            " Algorithm options:",
            "",
        ]
        for key, value in self.algorithm_options.items():
            lines.append(f"   {key + ':':<28} {value}")
        lines += [
            "",
            f" Number of processes:     {self.nprocs}",
            f" Worker timeout:          {self.timeout_seconds} seconds",
            f" Output mode:             {self.output_mode}",
        ]
        return lines

    def _render_generation(self) -> list:
        lines = self._section("COMBINATION GENERATION")
        lines += [
            f" Skeletons loaded:        {self.loaded_skeletons}",
            f" Substituents loaded:     {self.loaded_substituents}",
            f" Attachment sites:        {self.attachment_site_count}",
            f" Total combinations:      {self.total_combinations}",
            "",
            " Combination distribution:",
            "",
            "   Skeleton                          Combinations",
            "   --------------------------------  ------------",
        ]
        for label, count in self.per_skeleton_counts:
            lines.append(f"   {label:<32}  {count}")
        return lines

    def _render_execution(self) -> list:
        lines = self._section("EXECUTION RESULTS")
        lines += [
            "   No.  Status        Time/s    Combination",
            "   ---  ------------  --------  " + "-" * 48,
        ]
        for r in self.results:
            dur = (
                f"{r.duration_seconds:.2f}"
                if r.duration_seconds is not None
                else "-"
            )
            lines.append(
                f"   {r.combination_number:>3}  "
                f"{r.execution_status:<12}  {dur:>8}  {r.label}"
            )
        return lines

    def _render_failed(self) -> list:
        lines = self._section("FAILED COMBINATIONS")
        failed = [
            r for r in self.results if r.execution_status == STATUS_FAILED
        ]
        if not failed:
            lines.append(" None.")
            return lines
        for r in failed:
            dur = (
                f"{r.duration_seconds:.2f}"
                if r.duration_seconds is not None
                else "-"
            )
            lines += [
                f" [Combination {r.combination_number}]",
                "",
                f"   Label:                 {r.label}",
                f"   Failure stage:         {r.failure_stage}",
                f"   Algorithm:             {self.algorithm_name}",
                f"   Elapsed time:          {dur} seconds",
                f"   Error type:            {r.error_type}",
                f"   Reason:                {r.error_message}",
                "",
            ]
        return lines

    def _render_timed_out(self) -> list:
        lines = self._section("TIMED-OUT COMBINATIONS")
        timed = [
            r for r in self.results if r.execution_status == STATUS_TIMED_OUT
        ]
        if not timed:
            lines.append(" None.")
            return lines
        for r in timed:
            dur = (
                f"{r.duration_seconds:.2f}"
                if r.duration_seconds is not None
                else "-"
            )
            lines += [
                f" [Combination {r.combination_number}]",
                "",
                f"   Label:                 {r.label}",
                f"   Elapsed time:          {dur} seconds",
                f"   Timeout threshold:     {self.timeout_seconds} seconds",
                "",
            ]
        return lines

    def _render_write_failures(self) -> list:
        lines = self._section("OUTPUT WRITE FAILURES")
        write_failed = [
            r
            for r in self.results
            if r.execution_status == STATUS_WRITE_FAILED
        ]
        if not write_failed:
            lines.append(" None.")
            return lines
        for r in write_failed:
            lines += [
                f" [Combination {r.combination_number}]",
                "",
                f"   Label:                 {r.label}",
                f"   Output mode:           {self.output_mode}",
                f"   Target path:           {r.output_path}",
                f"   Error type:            {r.error_type}",
                f"   Reason:                {r.error_message}",
                "",
            ]
        return lines

    def _render_output_structures(self, stats: dict) -> list:
        lines = self._section("OUTPUT STRUCTURES")
        lines += [
            f" Output mode:             {self.output_mode}",
            f" Output path/directory:   {self.output_location or 'N/A'}",
            f" Structures written:      {stats['write_succeeded']}",
            "",
        ]
        written = [
            r
            for r in self.results
            if r.execution_status == STATUS_SUCCESS and r.output_path
        ]
        if self.output_mode == "merged":
            lines += [
                " Merged output mapping:",
                "",
                "   XYZ index  Combination",
                "   ---------  " + "-" * 48,
            ]
            if not written:
                lines.append(" None.")
            for r in written:
                lines.append(f"   {str(r.structure_index):>9}  {r.label}")
        else:
            lines += [
                " Separate output mapping:",
                "",
                "   Combination                               Output file",
                "   ----------------------------------------  " + "-" * 32,
            ]
            if not written:
                lines.append(" None.")
            for r in written:
                lines.append(f"   {r.label:<40}  {r.output_path}")
        return lines

    def _render_final_summary(
        self, stats: dict, error_codes: list, exit_code: int
    ) -> list:
        lines = self._section("FINAL SUMMARY")
        failed_numbers = _numbers(self.results, STATUS_FAILED)
        timed_numbers = _numbers(self.results, STATUS_TIMED_OUT)
        write_numbers = _numbers(self.results, STATUS_WRITE_FAILED)
        lines += [
            " Combination execution:",
            "",
            f"   Total combinations:          {stats['total']}",
            f"   Generated successfully:      {stats['generated']}",
            f"   Failed:                      {stats['failed']}",
            f"   Timed out:                   {stats['timed_out']}",
            "",
            " Output delivery:",
            "",
            f"   Structures intended:         {stats['generated']}",
            f"   Structures written:          {stats['write_succeeded']}",
            f"   Structures not delivered:    {stats['write_failed']}",
            "",
            " Problem combination numbers:",
            "",
            f"   Failed:                      {failed_numbers}",
            f"   Timed out:                   {timed_numbers}",
            f"   Write failed:                {write_numbers}",
            "",
            f"   Input load errors:           {stats['input_errors']}",
            "",
            " Error codes:",
            "",
        ]
        if error_codes:
            for code in error_codes:
                desc = ERROR_CODE_DESCRIPTIONS.get(code, "")
                lines.append(f"   {code:<20} {desc}".rstrip())
        else:
            lines.append("   None.")
        lines += [
            "",
            f" Exit code:                     {exit_code}",
            f" Total elapsed time:            "
            f"{_format_duration(self.duration_seconds)}",
        ]
        if self.error_message:
            lines.append(
                f" Error detail:                  {self.error_message}"
            )
        lines += ["", " " + self._termination_line(exit_code)]
        return lines

    def _termination_line(self, exit_code: int) -> str:
        """Return the fixed-format final termination line."""
        ts = self.finished_at.strftime("%Y-%m-%d %H:%M:%S")
        if exit_code == 0:
            return f"Normal termination of CHEMSMART Iterate at {ts}."
        return f"Error termination of CHEMSMART Iterate at {ts}."


def _numbers(results: list, status: str) -> str:
    nums = [
        r.combination_number for r in results if r.execution_status == status
    ]
    return ", ".join(str(n) for n in nums) if nums else "None"


def _format_list(value) -> str:
    if not value:
        return "N/A"
    if isinstance(value, (list, tuple)):
        return ", ".join(str(v) for v in value)
    return str(value)


def _format_sites(entry: dict) -> str:
    slots = entry.get("slots")
    if slots:
        parts = [
            f"group {s['group']}: [{_format_list(s['link_indices'])}]"
            for s in slots
        ]
        return "; ".join(parts)
    return _format_list(entry.get("link_index"))


def _format_duration(seconds: Optional[float]) -> str:
    if seconds is None:
        return "N/A"
    if seconds < 60:
        return f"{seconds:.2f} s"
    minutes, secs = divmod(seconds, 60)
    if minutes < 60:
        return f"{int(minutes)}m {secs:.1f}s"
    hours, minutes = divmod(int(minutes), 60)
    return f"{hours}h {int(minutes)}m {secs:.1f}s"


def write_report_atomically(path: str, text: str) -> None:
    """Write ``text`` to ``path`` atomically via a temp file + os.replace.

    A partial report is never left behind: the content is fully written to a
    temporary file in the same directory and then atomically renamed over the
    destination (overwriting any pre-existing file). Raises on failure so the
    caller can surface a clear terminal error.
    """
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    tmp_path = f"{path}.tmp.{os.getpid()}"
    try:
        with open(tmp_path, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, path)
    except Exception:
        if os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        raise
