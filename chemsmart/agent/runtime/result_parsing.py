"""Native calculation-result parsing with optional cclib enrichment."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from chemsmart.agent.runtime.cclib_parser import parse_cclib_output
from chemsmart.agent.runtime.output_io import read_output_window


def inspect_output(
    output_path: str | Path,
    *,
    program: str = "",
    kind: str = "",
) -> dict[str, Any]:
    path = Path(output_path) if output_path else Path()
    if not output_path or not path.is_file():
        return {}
    text = read_output_window(path)
    detected = program.lower() or _detect_program(text, path)
    result: dict[str, Any] = {
        "output_path": str(path.resolve()),
        "program": detected,
        "kind": kind,
    }
    if detected == "orca":
        energy_matches = _parse_orca(result, text)
    elif detected == "xtb":
        energy_matches = _parse_xtb(result, text)
    else:
        energy_matches = _parse_gaussian(result, text)
    if energy_matches:
        energies = [float(value) for value in energy_matches]
        result["energy_min"] = min(energies)
        result["energy_max"] = max(energies)
    result["path_direction"] = _path_direction(text)
    result["chemistry_elapsed_s"] = _program_runtime_seconds(text, detected)
    _merge_cclib_summary(result, path)
    return result


def _parse_orca(result: dict[str, Any], text: str) -> list[str]:
    energy_matches = re.findall(
        r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)",
        text,
    )
    scf_matches = re.findall(
        r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES",
        text,
        re.IGNORECASE,
    )
    optimization_cycles = re.findall(
        r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)",
        text,
    )
    scan_points = re.findall(
        r"(?:RELAXED SURFACE SCAN (?:STEP|POINT)|SCAN POINT)\s+(\d+)",
        text,
        re.IGNORECASE,
    )
    scan_totals = re.findall(
        r"SCAN POINT\s+\d+\s+(?:OF|/)\s+(\d+)",
        text,
        re.IGNORECASE,
    )
    neb_images = re.findall(
        r"Number of images \(incl\. end points\)\s+\.{2,}\s+(\d+)",
        text,
    )
    qmmm_total = re.findall(r"Size of QMMM System\s+\.{2,}\s+(\d+)", text)
    qmmm_mm = re.findall(
        r"Point charges in QM calc\. from MM atoms\.{2,}\s+(\d+)",
        text,
    )
    result.update(
        energy=float(energy_matches[-1]) if energy_matches else None,
        scf_cycles=int(scf_matches[-1]) if scf_matches else None,
        optimization_cycles=(
            int(optimization_cycles[-1]) if optimization_cycles else None
        ),
        optimization_converged=(
            "THE OPTIMIZATION HAS CONVERGED" in text
            if optimization_cycles
            else None
        ),
        scan_points_completed=(
            max(map(int, scan_points)) if scan_points else None
        ),
        scan_points_total=(int(scan_totals[-1]) if scan_totals else None),
        neb_images=int(neb_images[-1]) if neb_images else None,
        qmmm_total_atoms=int(qmmm_total[-1]) if qmmm_total else None,
        qmmm_mm_atoms=int(qmmm_mm[-1]) if qmmm_mm else None,
        normal_termination="ORCA TERMINATED NORMALLY" in text,
    )
    if (
        result.get("qmmm_total_atoms") is not None
        and result.get("qmmm_mm_atoms") is not None
    ):
        result["qmmm_qm_atoms"] = int(result["qmmm_total_atoms"]) - int(
            result["qmmm_mm_atoms"]
        )
    frequency_text = _last_section(text, "VIBRATIONAL FREQUENCIES")
    frequencies = [
        float(value)
        for value in re.findall(
            r"^\s*\d+:\s+(-?\d+(?:\.\d+)?)\s+cm\*\*-1",
            frequency_text,
            re.MULTILINE,
        )
        if abs(float(value)) > 1e-6
    ]
    _add_frequency_summary(result, frequencies)
    return energy_matches


def _parse_gaussian(result: dict[str, Any], text: str) -> list[str]:
    energy_matches = re.findall(
        r"SCF Done:\s+E\([^)]*\)\s*=\s*(-?\d+\.\d+)",
        text,
    )
    scf_matches = re.findall(r"SCF Done:.*?after\s+(\d+)\s+cycles", text)
    cycle_matches = re.findall(r"Step number\s+(\d+)", text)
    if not cycle_matches:
        cycle_matches = re.findall(r"Cycle\s+(\d+)", text)
    scan_points = re.findall(
        r"Scan point\s+(\d+)(?:\s+out of\s+(\d+))?",
        text,
        re.IGNORECASE,
    )
    result.update(
        energy=float(energy_matches[-1]) if energy_matches else None,
        scf_cycles=int(scf_matches[-1]) if scf_matches else None,
        optimization_cycles=int(cycle_matches[-1]) if cycle_matches else None,
        optimization_converged=(
            "Stationary point found" in text if cycle_matches else None
        ),
        scan_points_completed=(
            max(int(point) for point, _ in scan_points)
            if scan_points
            else None
        ),
        scan_points_total=(
            max(int(total) for _, total in scan_points if total)
            if any(total for _, total in scan_points)
            else None
        ),
        normal_termination="Normal termination of Gaussian" in text,
    )
    frequency_text = _last_section(text, "Harmonic frequencies")
    frequencies = [
        float(value)
        for line in re.findall(r"Frequencies\s+--\s+([^\n]+)", frequency_text)
        for value in re.findall(r"-?\d+(?:\.\d+)?", line)
    ]
    _add_frequency_summary(result, frequencies)
    return energy_matches


def _parse_xtb(result: dict[str, Any], text: str) -> list[str]:
    # xTB writes no Gaussian/ORCA termination banner; a completed run ends with
    # a "* finished run" line (the same marker XTBFile uses), and the fake
    # runner emits it too. Routing xtb through _parse_gaussian would report
    # normal_termination=False on every successful xtb job.
    energy_matches = re.findall(
        r"TOTAL ENERGY\s+(-?\d+\.\d+)\s*Eh",
        text,
    )
    opt_iterations = re.findall(
        r"GEOMETRY OPTIMIZATION CONVERGED AFTER\s+(\d+)\s+ITERATION",
        text,
        re.IGNORECASE,
    )
    result.update(
        energy=float(energy_matches[-1]) if energy_matches else None,
        scf_cycles=None,
        optimization_cycles=(
            int(opt_iterations[-1]) if opt_iterations else None
        ),
        optimization_converged=(bool(opt_iterations) or None),
        scan_points_completed=None,
        scan_points_total=None,
        normal_termination="* finished run" in text,
    )
    _add_frequency_summary(result, _xtb_frequencies(text))
    return energy_matches


def _xtb_frequencies(text: str) -> list[float]:
    # xTB prints all 3N mass-weighted-Hessian eigenvalues on "eigval :" lines
    # (near-zero translational/rotational modes first, then the real
    # vibrations), not one frequency per numbered line as Gaussian/ORCA do.
    section = _last_section(text, "projected vibrational frequencies")
    values: list[float] = []
    started = False
    for line in section.splitlines():
        stripped = line.strip()
        if stripped.startswith("eigval"):
            started = True
            values.extend(
                float(match) for match in re.findall(r"-?\d+\.\d+", stripped)
            )
        elif started:
            break
    return [value for value in values if abs(value) > 1e-6]


def _last_section(text: str, marker: str) -> str:
    start = text.rfind(marker)
    return text[start:] if start >= 0 else text


def _add_frequency_summary(
    result: dict[str, Any],
    frequencies: list[float],
) -> None:
    if not frequencies:
        return
    imaginary = [value for value in frequencies if value < 0]
    result["frequencies"] = frequencies[:12]
    result["imag_freqs"] = imaginary[:6]
    result["frequency_count"] = len(frequencies)
    result["imaginary_frequency_count"] = len(imaginary)


def _merge_cclib_summary(result: dict[str, Any], path: Path) -> None:
    cclib_result = parse_cclib_output(path)
    if cclib_result.get("status") != "ok":
        result["parser_backend"] = "native"
        warning = cclib_result.get("warning")
        result["parser_warnings"] = [str(warning)] if warning else []
        return

    result["parser_backend"] = "native+cclib"
    result["parser_warnings"] = []
    for key in (
        "energy",
        "scf_cycles",
        "optimization_converged",
        "frequency_count",
        "imaginary_frequency_count",
    ):
        if result.get(key) is None and cclib_result.get(key) is not None:
            result[key] = cclib_result[key]
    if not result.get("frequencies") and cclib_result.get("frequencies"):
        frequencies = list(cclib_result["frequencies"])
        result["frequencies"] = frequencies[:12]
        result["imag_freqs"] = [value for value in frequencies if value < 0][
            :6
        ]
    result["atom_count"] = cclib_result.get("atom_count")
    result["parsed_charge"] = cclib_result.get("charge")
    result["parsed_multiplicity"] = cclib_result.get("multiplicity")


def _path_direction(text: str) -> str:
    lowered = text.lower()
    has_forward = bool(
        re.search(r"\b(?:irc|qrc)?\s*forward direction\b", lowered)
    )
    has_reverse = bool(
        re.search(
            r"\b(?:irc|qrc)?\s*(?:reverse|backward) direction\b",
            lowered,
        )
    )
    if has_forward and has_reverse:
        return "both"
    if has_forward:
        return "forward"
    if has_reverse:
        return "reverse"
    return ""


def _program_runtime_seconds(text: str, program: str) -> float | None:
    if program == "orca":
        matches = re.findall(
            r"TOTAL RUN TIME:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+"
            r"(\d+) seconds\s+(\d+) msec",
            text,
        )
        if not matches:
            return None
        days, hours, minutes, seconds, milliseconds = map(int, matches[-1])
        return (
            days * 86_400
            + hours * 3_600
            + minutes * 60
            + seconds
            + milliseconds / 1_000
        )
    matches = re.findall(
        r"Elapsed time:\s+(\d+) days\s+(\d+) hours\s+(\d+) minutes\s+"
        r"([\d.]+) seconds",
        text,
    )
    if not matches:
        return None
    days, hours, minutes, seconds = matches[-1]
    return (
        int(days) * 86_400
        + int(hours) * 3_600
        + int(minutes) * 60
        + float(seconds)
    )


def _detect_program(text: str, path: Path) -> str:
    if "O   R   C   A" in text or "ORCA TERMINATED" in text:
        return "orca"
    if "Gaussian" in text or path.suffix.lower() == ".log":
        return "gaussian"
    return ""


__all__ = ["inspect_output"]
