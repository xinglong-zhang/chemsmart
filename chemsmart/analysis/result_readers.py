"""One result reader per program, behind a single semantic vocabulary.

Typed quantity extraction was written against PySCF because ChemSmart controls
its structured HDF5 schema, and its contracts were deliberately left
program-neutral "so that Gaussian, ORCA, and xTB readers can be registered
later without changing the expression evaluator".  Until they are, every
observable from those programs has to be read by eye out of a log, which is
exactly the channel the project-YAML hub exists to close: a number a model
typed is not a number ChemSmart measured.

Registering them here needs no new parser.  Each program's reader already
exposes the same few accessors under the same names -- ``energies``,
``gibbs_free_energy``, ``vibrational_frequencies``, ``molecule`` -- so a
selector maps to an attribute rather than to a per-program branch ladder.  A
program is supported exactly when it appears in ``RESULT_READERS``; nothing
else in the extraction path needs to know which programs exist.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

__all__ = [
    "RESULT_READERS",
    "MissingQuantityError",
    "ResultReaderV1",
    "SELECTOR_UNITS",
    "reader_for",
    "registered_reader_programs",
    "registered_reader_selectors",
]


class MissingQuantityError(LookupError):
    """The run did not produce this quantity.

    Distinct from an unsupported selector: a single point genuinely has no
    Gibbs energy, and reporting that honestly is what lets a caller decide to
    run a frequency job rather than treat a parse failure as a bug.
    """


#: The canonical unit each semantic selector is reported in.  A reader that
#: cannot produce a quantity in its declared unit must fail rather than guess:
#: a silently mis-scaled energy is worse than an absent one.
SELECTOR_UNITS = {
    "absorption_wavelengths": "nm",
    "energy": "Eh",
    "energies": "Eh",
    "entropy_times_temperature": "Eh",
    "excitation_energies": "eV",
    "singlet_excitation_energies": "eV",
    "triplet_excitation_energies": "eV",
    "excited_state_indices": "",
    "excited_state_manifold_roots": "",
    "excited_state_multiplicities": "",
    "excited_state_labels": "",
    "excited_state_spin_square": "",
    "gibbs_free_energy": "Eh",
    "oscillator_strengths": "",
    "singlet_oscillator_strengths": "",
    "triplet_oscillator_strengths": "",
    "vibrational_frequencies": "cm^-1",
    "vpt2_harmonic_frequencies": "cm^-1",
    "vpt2_fundamental_frequencies": "cm^-1",
    "vpt2_zero_point_rovibrational_energy": "cm^-1",
    "positions": "Angstrom",
    "connectivity": "",
    "symbols": "",
    "charge": "",
    "multiplicity": "",
    "scf_energy": "Eh",
    "correlation_energy": "Eh",
    "dipole_moment": "Debye",
    "dipole_moment_magnitude": "Debye",
    "homo": "eV",
    "lumo": "eV",
    "gap": "eV",
    "spin_square": "",
    "spin_square_after_annihilation": "",
    "spin_square_target": "",
    "spin_square_deviation": "",
    "effective_multiplicity": "",
    "wavefunction_stability_verdict": "",
    "wavefunction_stability_history": "",
    "trajectory_frame_count": "",
    "trajectory_start_positions": "Angstrom",
    "trajectory_end_positions": "Angstrom",
    "trajectory_start_connectivity": "",
    "trajectory_end_connectivity": "",
    "trajectory_connectivity_changed": "",
    "irc_direction": "",
}


def _last_energy(output: Any) -> float:
    """Return the converged energy, preferring an explicit final value."""

    final = getattr(output, "final_energy", None)
    if final is not None:
        return float(final)
    energies = getattr(output, "energies", None)
    if not energies:
        raise ValueError("result exposes no energy")
    return float(energies[-1])


def _positions(output: Any) -> list[list[float]]:
    molecule = output.molecule
    return [[float(value) for value in row] for row in molecule.positions]


def _symbols(output: Any) -> list[str]:
    return [str(item) for item in output.molecule.chemical_symbols]


def _required_record_values(
    records: list[dict[str, Any]], key: str, description: str
) -> list[Any]:
    if not records or any(record.get(key) is None for record in records):
        raise MissingQuantityError(
            f"result does not establish {description} for every excited state"
        )
    return [record[key] for record in records]


def _orca_excitation_energies(output: Any) -> list[float]:
    records = list(output.excited_state_records or ())
    if records:
        return [float(item["energy_eV"]) for item in records]
    # ORCA 6 spectrum-only fragments remain useful for the legacy aggregate
    # selector, but they cannot support multiplicity-specific selectors.
    return [float(item) for item in output.excitation_energies_eV]


def _spin_square_target(output: Any) -> float:
    multiplicity = getattr(output, "multiplicity", None)
    if not isinstance(multiplicity, int) or multiplicity <= 0:
        raise MissingQuantityError(
            "result does not establish a positive integer multiplicity"
        )
    return (float(multiplicity) ** 2 - 1.0) / 4.0


def _effective_multiplicity(spin_square: float) -> float:
    if spin_square < 0.0:
        raise MissingQuantityError("negative <S^2> cannot define multiplicity")
    return float((1.0 + 4.0 * spin_square) ** 0.5)


def _connectivity_matrix(molecule: Any) -> list[list[int]]:
    """Return geometry-perceived covalent connectivity in source atom order.

    This is intentionally a binary connectivity observation rather than an
    asserted electronic bond order.  IRC endpoints are not fully optimized,
    and covalent-radius perception is appropriate for deciding whether two
    path ends have different molecular graphs while leaving chemical
    interpretation to the scientist.
    """

    graph = molecule.to_graph(bond_cutoff_buffer=0.05, adjust_H=True)
    size = int(molecule.num_atoms)
    matrix = [[0 for _ in range(size)] for _ in range(size)]
    for first, second in graph.edges:
        matrix[int(first)][int(second)] = 1
        matrix[int(second)][int(first)] = 1
    return matrix


def _irc_structures(output: Any) -> list[Any]:
    jobtype = str(getattr(output, "jobtype", "") or "").casefold()
    if jobtype not in {"irc", "ircf", "ircr"}:
        raise MissingQuantityError(
            "trajectory selectors require an IRC program result"
        )
    structures = list(getattr(output, "all_structures", ()) or ())
    if len(structures) < 2:
        raise MissingQuantityError(
            "IRC output does not contain at least two parsed trajectory frames"
        )
    symbols = tuple(structures[0].chemical_symbols)
    if any(tuple(item.chemical_symbols) != symbols for item in structures[1:]):
        raise MissingQuantityError(
            "IRC trajectory changes atom identity or order"
        )
    return structures


def _irc_direction(output: Any) -> str:
    jobtype = str(getattr(output, "jobtype", "") or "").casefold()
    directions = {"ircf": "forward", "ircr": "reverse", "irc": "combined"}
    if jobtype not in directions:
        raise MissingQuantityError(
            "result does not establish an IRC direction"
        )
    return directions[jobtype]


def _orca_irc_direction(output: Any) -> str:
    """Return the direction explicitly echoed by an ORCA IRC input block."""

    if str(getattr(output, "jobtype", "") or "").casefold() != "irc":
        raise MissingQuantityError("result is not an ORCA IRC calculation")
    observed = getattr(output, "irc_direction", None)
    if observed is None:
        # Preserve the result-reader protocol for light-weight parser fixtures
        # while production ORCAInput/ORCAOutput objects use the shared typed
        # ``irc_direction`` property above.
        pattern = re.compile(
            r"\bDirection\s+(both|forward|backward|down)\b",
            re.IGNORECASE,
        )
        matches = [
            match.group(1).casefold()
            for line in getattr(output, "contents", ())
            if (match := pattern.search(str(line))) is not None
        ]
        observed = matches[-1] if matches else None
    if observed not in {"backward", "both", "down", "forward"}:
        raise MissingQuantityError(
            "ORCA output does not explicitly establish the IRC direction"
        )
    return observed


def _trajectory_connectivity_changed(output: Any) -> int:
    structures = _irc_structures(output)
    return int(
        _connectivity_matrix(structures[0])
        != _connectivity_matrix(structures[-1])
    )


def _xyz_trajectory_structures(output: Any) -> list[Any]:
    """Return a multi-frame XYZ trajectory without assigning program meaning.

    ORCA writes IRC paths to XYZ sidecars, which already enter Runtime V2 as
    ``geometry_xyz`` artifacts.  The file itself establishes ordered frames,
    but not whether a path is forward, reverse, or the combined IRC.  Generic
    trajectory selectors therefore expose only observations available from
    the registered bytes; ``irc_direction`` remains a program-log selector.
    """

    structures = list(output.get_molecules(index=":", return_list=True) or ())
    if len(structures) < 2:
        raise MissingQuantityError(
            "XYZ artifact does not contain at least two trajectory frames"
        )
    symbols = tuple(structures[0].chemical_symbols)
    if any(tuple(item.chemical_symbols) != symbols for item in structures[1:]):
        raise MissingQuantityError(
            "XYZ trajectory changes atom identity or order"
        )
    return structures


def _xyz_trajectory_connectivity_changed(output: Any) -> int:
    structures = _xyz_trajectory_structures(output)
    return int(
        _connectivity_matrix(structures[0])
        != _connectivity_matrix(structures[-1])
    )


@dataclass(frozen=True)
class ResultReaderV1:
    """How one program's result file answers the shared selector vocabulary."""

    program: str
    #: Host artifact kind this reader consumes; extraction refuses any other.
    artifact_kind: str
    parser_id: str
    #: Constructs the parser from an already host-verified path.
    open_output: Callable[[Path], Any]
    #: Selector name to a callable reading it from the parser object.
    accessors: dict[str, Callable[[Any], Any]]
    #: Native program outputs must prove normal termination.  A standalone
    #: geometry artifact is data rather than an engine run, so that format can
    #: opt out while retaining the same typed quantity path.
    requires_normal_termination: bool = True

    @property
    def selectors(self) -> frozenset[str]:
        return frozenset(self.accessors)

    def read(self, output: Any, selector: str) -> tuple[Any, str]:
        """Return ``(value, unit)`` for ``selector``.

        A selector this reader does not implement is refused by name, and a
        quantity this particular run did not produce -- a single point has no
        Gibbs energy and no optimised geometry -- is refused as absent.  Both
        are stated gaps rather than a silent ``None`` or a parser traceback
        leaking out of the extraction path.
        """

        if selector not in self.accessors:
            raise ValueError(
                f"{self.program} result reader does not provide {selector!r}; "
                f"it provides {sorted(self.selectors)}"
            )
        try:
            value = self.accessors[selector](output)
        except MissingQuantityError:
            raise
        except Exception as exc:
            # The parsers raise IndexError/TypeError when a block the run never
            # wrote is requested. That is an absent quantity, not a defect.
            raise MissingQuantityError(
                f"{self.program} result contains no {selector!r} value "
                f"({type(exc).__name__})"
            ) from exc
        if value is None or (isinstance(value, (list, tuple)) and not value):
            raise MissingQuantityError(
                f"{self.program} result contains no {selector!r} value"
            )
        return value, SELECTOR_UNITS[selector]


def _text_output_accessors(
    *, thermochemistry: bool = True
) -> dict[str, Callable[[Any], Any]]:
    """Accessors shared by the log-parsing programs."""

    accessors: dict[str, Callable[[Any], Any]] = {
        "energy": _last_energy,
        "energies": lambda output: [float(item) for item in output.energies],
        "vibrational_frequencies": lambda output: [
            float(item) for item in output.vibrational_frequencies
        ],
        "positions": _positions,
        "connectivity": lambda output: _connectivity_matrix(output.molecule),
        "symbols": _symbols,
        "charge": lambda output: int(output.charge),
        "multiplicity": lambda output: int(output.multiplicity),
    }
    if thermochemistry:
        accessors["gibbs_free_energy"] = lambda output: float(
            output.gibbs_free_energy
        )
    return accessors


def _orca_output(path: Path) -> Any:
    from chemsmart.io.orca.output import ORCAOutput

    return ORCAOutput(filename=str(path))


def _orca_total_energy(output: Any) -> float:
    """Return the final total ORCA energy, including post-SCF correlation.

    ``ORCAOutput.final_energy`` historically prefers ``final_scf_energy``.
    That is suitable for DFT but silently drops the correlation contribution
    from MP2 and coupled-cluster results.  The last explicit ORCA final-energy
    record is the program's total result and remains correct for DFT jobs too.
    """

    values: list[float] = []
    for line in getattr(output, "contents", ()):
        if "FINAL SINGLE POINT ENERGY" not in str(line).upper():
            continue
        try:
            values.append(float(str(line).split()[-1]))
        except (TypeError, ValueError):
            continue
    if values:
        return values[-1]
    return _last_energy(output)


def _orca_scf_energy(output: Any) -> float:
    # ``ORCAOutput.final_scf_energy`` currently treats an empty optimization
    # slice as an optimization result and can therefore return ``None`` for a
    # perfectly valid correlated single point.  Read the last explicit TOTAL
    # SCF ENERGY value first; for post-HF output that is the reference energy
    # paired with the final correlated total below it.
    values: list[float] = []
    for line in getattr(output, "contents", ()):
        text = str(line)
        if "Total Energy" not in text or ":" not in text:
            continue
        fields = text.split()
        try:
            colon = fields.index(":")
            values.append(float(fields[colon + 1]))
        except (ValueError, IndexError):
            continue
    if values:
        return values[-1]
    value = getattr(output, "final_scf_energy", None)
    if value is None:
        raise MissingQuantityError("ORCA result contains no final SCF energy")
    return float(value)


def _orca_correlation_energy(output: Any) -> float:
    """Return total post-SCF correlation as ``E(total) - E(SCF)``."""

    return _orca_total_energy(output) - _orca_scf_energy(output)


def _orca_positions(output: Any) -> list[list[float]]:
    molecule = output.thermochemistry_molecule
    return [[float(value) for value in row] for row in molecule.positions]


def _orca_symbols(output: Any) -> list[str]:
    return [
        str(item) for item in output.thermochemistry_molecule.chemical_symbols
    ]


def _orca_accessors() -> dict[str, Callable[[Any], Any]]:
    accessors = _text_output_accessors()
    accessors.update(
        {
            "absorption_wavelengths": lambda output: [
                float(item) for item in output.absorption_wavelengths
            ],
            "excitation_energies": _orca_excitation_energies,
            "oscillator_strengths": lambda output: [
                float(item) for item in output.oscillator_strengths
            ],
            "excited_state_indices": lambda output: [
                int(item["state_index"])
                for item in output.excited_state_records
            ],
            "excited_state_manifold_roots": lambda output: [
                int(item["manifold_root"])
                for item in output.excited_state_records
            ],
            "excited_state_multiplicities": lambda output: [
                int(item)
                for item in _required_record_values(
                    output.excited_state_records,
                    "multiplicity",
                    "a spin multiplicity",
                )
            ],
            "excited_state_spin_square": lambda output: [
                float(item["spin_square"])
                for item in output.excited_state_records
            ],
            "singlet_excitation_energies": lambda output: [
                float(item["energy_eV"])
                for item in output.excited_state_records
                if item["multiplicity"] == 1
            ],
            "triplet_excitation_energies": lambda output: [
                float(item["energy_eV"])
                for item in output.excited_state_records
                if item["multiplicity"] == 3
            ],
            "singlet_oscillator_strengths": lambda output: [
                float(item["oscillator_strength"])
                for item in output.electronic_absorption_transition_records
                if item["multiplicity"] == 1
            ],
            "triplet_oscillator_strengths": lambda output: [
                float(item["oscillator_strength"])
                for item in output.electronic_absorption_transition_records
                if item["multiplicity"] == 3
            ],
            "energy": _orca_total_energy,
            "entropy_times_temperature": lambda output: float(
                output.entropy_times_temperature
            ),
            "positions": _orca_positions,
            "connectivity": lambda output: _connectivity_matrix(
                output.thermochemistry_molecule
            ),
            "symbols": _orca_symbols,
            "charge": lambda output: int(output.thermochemistry_charge),
            "multiplicity": lambda output: int(
                output.thermochemistry_multiplicity
            ),
            "scf_energy": _orca_scf_energy,
            "correlation_energy": _orca_correlation_energy,
            "dipole_moment": lambda output: [
                float(item)
                for item in output.dipole_moment_in_debye.reshape(-1)
            ],
            "dipole_moment_magnitude": lambda output: float(
                output.dipole_moment_magnitude_in_debye
            ),
            "vpt2_harmonic_frequencies": lambda output: [
                float(item) for item in output.vpt2_harmonic_frequencies
            ],
            "vpt2_fundamental_frequencies": lambda output: [
                float(item) for item in output.vpt2_fundamental_frequencies
            ],
            "vpt2_zero_point_rovibrational_energy": lambda output: float(
                output.vpt2_zero_point_rovibrational_energy
            ),
            "spin_square": lambda output: float(
                output.spin_square_history[-1]
            ),
            "spin_square_target": _spin_square_target,
            "spin_square_deviation": lambda output: float(
                output.spin_square_history[-1] - _spin_square_target(output)
            ),
            "effective_multiplicity": lambda output: _effective_multiplicity(
                float(output.spin_square_history[-1])
            ),
            "trajectory_frame_count": lambda output: len(
                _irc_structures(output)
            ),
            "trajectory_start_positions": lambda output: [
                [float(value) for value in row]
                for row in _irc_structures(output)[0].positions
            ],
            "trajectory_end_positions": lambda output: [
                [float(value) for value in row]
                for row in _irc_structures(output)[-1].positions
            ],
            "trajectory_start_connectivity": lambda output: (
                _connectivity_matrix(_irc_structures(output)[0])
            ),
            "trajectory_end_connectivity": lambda output: _connectivity_matrix(
                _irc_structures(output)[-1]
            ),
            "trajectory_connectivity_changed": (
                _trajectory_connectivity_changed
            ),
            "irc_direction": _orca_irc_direction,
        }
    )
    return accessors


def _gaussian_output(path: Path) -> Any:
    from chemsmart.io.gaussian.output import Gaussian16Output

    return Gaussian16Output(filename=str(path))


def _gaussian_accessors() -> dict[str, Callable[[Any], Any]]:
    accessors = _text_output_accessors()
    accessors.update(
        {
            "absorption_wavelengths": lambda output: [
                float(item) for item in output.absorptions_in_nm
            ],
            "excitation_energies": lambda output: [
                float(item["energy_eV"])
                for item in output.excited_state_records
            ],
            "oscillator_strengths": lambda output: [
                float(item) for item in output.oscillatory_strengths
            ],
            "excited_state_indices": lambda output: [
                int(item["state_index"])
                for item in output.excited_state_records
            ],
            "excited_state_manifold_roots": lambda output: [
                int(item)
                for item in _required_record_values(
                    output.excited_state_records,
                    "manifold_root",
                    "a manifold-local root index",
                )
            ],
            "excited_state_multiplicities": lambda output: [
                int(item)
                for item in _required_record_values(
                    output.excited_state_records,
                    "multiplicity",
                    "a source-labelled spin multiplicity",
                )
            ],
            "excited_state_labels": lambda output: [
                str(item["state_label"])
                for item in output.excited_state_records
            ],
            "excited_state_spin_square": lambda output: [
                float(item)
                for item in _required_record_values(
                    output.excited_state_records,
                    "spin_square",
                    "a printed excited-state <S^2>",
                )
            ],
            "singlet_excitation_energies": lambda output: [
                float(item["energy_eV"])
                for item in output.excited_state_records
                if item["multiplicity"] == 1
            ],
            "triplet_excitation_energies": lambda output: [
                float(item["energy_eV"])
                for item in output.excited_state_records
                if item["multiplicity"] == 3
            ],
            "singlet_oscillator_strengths": lambda output: [
                float(item["oscillator_strength"])
                for item in output.excited_state_records
                if item["multiplicity"] == 1
                and item["oscillator_strength"] is not None
            ],
            "triplet_oscillator_strengths": lambda output: [
                float(item["oscillator_strength"])
                for item in output.excited_state_records
                if item["multiplicity"] == 3
                and item["oscillator_strength"] is not None
            ],
            "dipole_moment": lambda output: [
                float(item) for item in output.all_dipole_moments[-1]
            ],
            "dipole_moment_magnitude": lambda output: float(
                output.all_dipole_moment_magnitudes[-1]
            ),
            "spin_square": lambda output: float(
                output.spin_square_history[-1]["before_annihilation"]
            ),
            "spin_square_after_annihilation": lambda output: float(
                output.spin_square_history[-1]["after_annihilation"]
            ),
            "spin_square_target": _spin_square_target,
            "spin_square_deviation": lambda output: float(
                output.spin_square_history[-1]["before_annihilation"]
                - _spin_square_target(output)
            ),
            "effective_multiplicity": lambda output: _effective_multiplicity(
                float(output.spin_square_history[-1]["before_annihilation"])
            ),
            "wavefunction_stability_verdict": lambda output: str(
                output.wavefunction_stability_history[-1]
            ),
            "wavefunction_stability_history": lambda output: [
                str(item) for item in output.wavefunction_stability_history
            ],
            "trajectory_frame_count": lambda output: len(
                _irc_structures(output)
            ),
            "trajectory_start_positions": lambda output: [
                [float(value) for value in row]
                for row in _irc_structures(output)[0].positions
            ],
            "trajectory_end_positions": lambda output: [
                [float(value) for value in row]
                for row in _irc_structures(output)[-1].positions
            ],
            "trajectory_start_connectivity": lambda output: (
                _connectivity_matrix(_irc_structures(output)[0])
            ),
            "trajectory_end_connectivity": lambda output: _connectivity_matrix(
                _irc_structures(output)[-1]
            ),
            "trajectory_connectivity_changed": (
                _trajectory_connectivity_changed
            ),
            "irc_direction": _irc_direction,
        }
    )
    return accessors


def _xtb_output(path: Path) -> Any:
    from chemsmart.io.xtb.output import XTBOutput

    # Unlike the Gaussian and ORCA readers, XTBOutput represents the complete
    # calculation directory because xTB can distribute one result across the
    # main log, geometry, Hessian, and vibrational-spectrum files.  The host
    # has already resolved and verified the exact main output artifact; its
    # parent is therefore the corresponding calculation directory.
    return XTBOutput(folder=str(path.parent))


def _xtb_accessors() -> dict[str, Callable[[Any], Any]]:
    # xTB's reader exposes no charge/multiplicity or energy series, so those
    # selectors stay unregistered rather than being faked from the geometry.
    return {
        "energy": _last_energy,
        "vibrational_frequencies": lambda output: [
            float(item) for item in output.vibrational_frequencies
        ],
        "positions": _positions,
        "connectivity": lambda output: _connectivity_matrix(output.molecule),
        "symbols": _symbols,
        "dipole_moment": lambda output: [
            float(item) for item in output.molecular_dipole_full
        ],
        "dipole_moment_magnitude": lambda output: float(
            output.total_molecular_dipole_moment
        ),
        "homo": lambda output: float(output.homo_energy),
        "lumo": lambda output: float(output.lumo_energy),
        "gap": lambda output: float(output.fmo_gap),
    }


def _xyz_output(path: Path) -> Any:
    """Open a registered XYZ as molecular data, not as a program run."""

    from chemsmart.io.xyz.xyzfile import XYZFile

    return XYZFile(filename=str(path))


_XYZ_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
_XYZ_HARTREE_PATTERNS = (
    re.compile(
        rf"\bEnergy\s*\(\s*Hartree\s*\)\s*:\s*({_XYZ_NUMBER})\b",
        re.IGNORECASE,
    ),
    re.compile(
        rf"^Coordinates from ORCA-job\b.*\sE\s+({_XYZ_NUMBER})\s*$",
        re.IGNORECASE,
    ),
)


def _xyz_energy(output: Any) -> float:
    """Read only an explicitly Hartree-grounded XYZ comment energy."""

    comment = str(output.comments or "").strip()
    for pattern in _XYZ_HARTREE_PATTERNS:
        match = pattern.search(comment)
        if match is not None:
            return float(match.group(1).replace("D", "E").replace("d", "e"))
    raise MissingQuantityError(
        "XYZ energy requires an explicit Energy(Hartree) label or the "
        "ChemSmart/ORCA 'Coordinates from ORCA-job ... E <value>' form"
    )


def _xyz_accessors() -> dict[str, Callable[[Any], Any]]:
    """Expose the quantities ChemSmart's XYZ parser actually establishes."""

    return {
        "energy": _xyz_energy,
        "positions": _positions,
        "connectivity": lambda output: _connectivity_matrix(output.molecule),
        "symbols": _symbols,
        "trajectory_frame_count": lambda output: len(
            _xyz_trajectory_structures(output)
        ),
        "trajectory_start_positions": lambda output: [
            [float(value) for value in row]
            for row in _xyz_trajectory_structures(output)[0].positions
        ],
        "trajectory_end_positions": lambda output: [
            [float(value) for value in row]
            for row in _xyz_trajectory_structures(output)[-1].positions
        ],
        "trajectory_start_connectivity": lambda output: _connectivity_matrix(
            _xyz_trajectory_structures(output)[0]
        ),
        "trajectory_end_connectivity": lambda output: _connectivity_matrix(
            _xyz_trajectory_structures(output)[-1]
        ),
        "trajectory_connectivity_changed": (
            _xyz_trajectory_connectivity_changed
        ),
    }


#: Programs whose results can be read into typed quantities.  PySCF keeps its
#: dedicated structured-HDF5 path in ``result_quantities``; the entries here are
#: the log-parsing programs that had none.
RESULT_READERS: dict[str, ResultReaderV1] = {
    "orca": ResultReaderV1(
        program="orca",
        artifact_kind="orca_output",
        parser_id="chemsmart.io.orca.output.ORCAOutput",
        open_output=_orca_output,
        accessors=_orca_accessors(),
    ),
    "gaussian": ResultReaderV1(
        program="gaussian",
        artifact_kind="gaussian_output",
        parser_id="chemsmart.io.gaussian.output.Gaussian16Output",
        open_output=_gaussian_output,
        accessors=_gaussian_accessors(),
    ),
    "xtb": ResultReaderV1(
        program="xtb",
        artifact_kind="xtb_output",
        parser_id="chemsmart.io.xtb.output.XTBOutput",
        open_output=_xtb_output,
        accessors=_xtb_accessors(),
    ),
    "xyz": ResultReaderV1(
        program="xyz",
        artifact_kind="geometry_xyz",
        parser_id="chemsmart.io.xyz.xyzfile.XYZFile",
        open_output=_xyz_output,
        accessors=_xyz_accessors(),
        requires_normal_termination=False,
    ),
}


#: Physical dimension of each selector, in the shared quantity vocabulary.
_SELECTOR_DIMENSIONS = {
    "absorption_wavelengths": "LENGTH",
    "energy": "ENERGY",
    "energies": "ENERGY",
    "entropy_times_temperature": "ENERGY",
    "excitation_energies": "ENERGY",
    "singlet_excitation_energies": "ENERGY",
    "triplet_excitation_energies": "ENERGY",
    "excited_state_indices": "DIMENSIONLESS",
    "excited_state_manifold_roots": "DIMENSIONLESS",
    "excited_state_multiplicities": "DIMENSIONLESS",
    "excited_state_labels": "DIMENSIONLESS",
    "excited_state_spin_square": "DIMENSIONLESS",
    "gibbs_free_energy": "ENERGY",
    "oscillator_strengths": "DIMENSIONLESS",
    "singlet_oscillator_strengths": "DIMENSIONLESS",
    "triplet_oscillator_strengths": "DIMENSIONLESS",
    "vibrational_frequencies": "FREQUENCY",
    "vpt2_harmonic_frequencies": "FREQUENCY",
    "vpt2_fundamental_frequencies": "FREQUENCY",
    "vpt2_zero_point_rovibrational_energy": "FREQUENCY",
    "positions": "LENGTH",
    "connectivity": "DIMENSIONLESS",
    "symbols": "DIMENSIONLESS",
    "charge": "DIMENSIONLESS",
    "multiplicity": "DIMENSIONLESS",
    "scf_energy": "ENERGY",
    "correlation_energy": "ENERGY",
    "dipole_moment": "DIPOLE_MOMENT",
    "dipole_moment_magnitude": "DIPOLE_MOMENT",
    "homo": "ENERGY",
    "lumo": "ENERGY",
    "gap": "ENERGY",
    "spin_square": "DIMENSIONLESS",
    "spin_square_after_annihilation": "DIMENSIONLESS",
    "spin_square_target": "DIMENSIONLESS",
    "spin_square_deviation": "DIMENSIONLESS",
    "effective_multiplicity": "DIMENSIONLESS",
    "wavefunction_stability_verdict": "DIMENSIONLESS",
    "wavefunction_stability_history": "DIMENSIONLESS",
    "trajectory_frame_count": "DIMENSIONLESS",
    "trajectory_start_positions": "LENGTH",
    "trajectory_end_positions": "LENGTH",
    "trajectory_start_connectivity": "DIMENSIONLESS",
    "trajectory_end_connectivity": "DIMENSIONLESS",
    "trajectory_connectivity_changed": "DIMENSIONLESS",
    "irc_direction": "DIMENSIONLESS",
}

_TEXT_SELECTORS = frozenset(
    {"irc_direction", "wavefunction_stability_verdict"}
)
_TEXT_VECTOR_SELECTORS = frozenset(
    {"excited_state_labels", "wavefunction_stability_history"}
)
_INTEGER_SELECTORS = frozenset(
    {
        "charge",
        "multiplicity",
        "trajectory_frame_count",
        "trajectory_connectivity_changed",
    }
)


def extract_logged_quantities(
    *,
    request: Any,
    artifact_path: str | Path,
) -> Any:
    """Extract selected quantities from a log-parsing program's result.

    Mirrors the structured-PySCF path: the artifact bytes are verified before
    and after parsing, so a file that changes mid-extraction is refused rather
    than silently mixing two runs.
    """

    from chemsmart.analysis import result_quantities as rq

    reader = reader_for(request.program)
    if reader is None:
        raise rq.QuantityContractError(
            f"no result reader is registered for {request.program!r}; "
            f"registered: {list(registered_reader_programs())}"
        )
    artifact = rq._verify_artifact(artifact_path, request.artifact_sha256)
    output = reader.open_output(artifact)
    if (
        reader.requires_normal_termination
        and getattr(output, "normal_termination", None) is not True
    ):
        raise rq.QuantityExtractionError(
            f"{request.program} scientific quantities require a normally "
            "terminated program result"
        )
    evidence_ref = f"artifact:{request.artifact_id}#{request.artifact_sha256}"
    quantities = []
    for selector in request.selectors:
        try:
            source_value, source_unit = reader.read(output, selector.selector)
        except MissingQuantityError as exc:
            raise rq.QuantityExtractionError(str(exc)) from exc
        dimension = getattr(rq, _SELECTOR_DIMENSIONS[selector.selector])
        if selector.selector in {"symbols", *_TEXT_VECTOR_SELECTORS}:
            value = source_value
            unit = source_unit
            data_kind = "text_vector"
        elif selector.selector in _TEXT_SELECTORS:
            value = source_value
            unit = source_unit
            data_kind = "text"
        elif selector.selector in _INTEGER_SELECTORS:
            value = int(source_value)
            unit = "1"
            data_kind = "integer"
        else:
            from chemsmart.analysis.quantity_expressions import (
                normalize_numeric_value,
            )

            value, unit, observed_dimension = normalize_numeric_value(
                source_value, source_unit
            )
            if observed_dimension != dimension:
                raise rq.QuantityExtractionError(
                    f"selector {selector.selector!r} produced an "
                    "incompatible unit"
                )
            data_kind = None
        quantities.append(
            rq.make_quantity_value(
                quantity_id=selector.quantity_id,
                source_value=source_value,
                source_unit=source_unit,
                value=value,
                unit=unit,
                dimension=dimension,
                evidence_ref=evidence_ref,
                data_kind=data_kind,
            )
        )
    if rq.result_file_sha256(artifact) != request.artifact_sha256:
        raise rq.QuantityExtractionError(
            "result artifact changed during extraction"
        )
    body = {
        "schema_version": "chemsmart.quantity-extraction-receipt.v1",
        "artifact_id": request.artifact_id,
        "artifact_sha256": request.artifact_sha256,
        "program": request.program,
        "parser_id": reader.parser_id,
        "quantities": tuple(quantities),
        "status": "extracted",
    }
    return rq.QuantityExtractionReceiptV1(
        **body, receipt_sha256=rq.canonical_quantity_sha256(body)
    )


def reader_for(program: str) -> ResultReaderV1 | None:
    """Return the registered reader for ``program``, or ``None``."""

    return RESULT_READERS.get(str(program).strip().lower())


def registered_reader_programs() -> tuple[str, ...]:
    """Return every program with a registered log reader, sorted."""

    return tuple(sorted(RESULT_READERS))


def registered_reader_selectors() -> dict[str, tuple[str, ...]]:
    """Return the model-discoverable selector inventory for each reader."""

    return {
        program: tuple(sorted(reader.selectors))
        for program, reader in sorted(RESULT_READERS.items())
    }
