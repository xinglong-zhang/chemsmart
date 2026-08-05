"""Reader for PySCF calculation results.

This is a typed reader over the structured ``label.h5`` file that ChemSmart's
own generated driver script writes -- not a log parser. ``io/orca/output.py``
is ~4,000 lines of regex because ORCA's printed text is the only interface it
offers; here ChemSmart controls both ends, so the numbers come back as exact
float64 arrays and the sibling ``label.out`` (PySCF's own log) is never
parsed. It exists for humans and for program sniffing only.

Thermochemistry is **delegated** to ``chemsmart/analysis/thermochemistry.py``
rather than recomputed. PySCF ships ``pyscf.hessian.thermo.thermo()``, but
using it would give free energies that silently disagree with the Gaussian
and ORCA numbers for the same molecule, since ChemSmart's engine has its own
rotational-symmetry and quasi-harmonic conventions. Cross-program consistency
is the whole point of a unified toolkit, so there is exactly one such engine.
"""

import json
import logging
import os
from datetime import datetime
from functools import cached_property

import numpy as np
from ase import units

from chemsmart.utils.constants import joule_per_mol_to_hartree
from chemsmart.utils.mixins import FileMixin

logger = logging.getLogger(__name__)

#: Layouts this reader understands. Version 1.0 used JSON string datasets for
#: three sections; 2.0 uses navigable groups and typed datasets throughout.
SUPPORTED_SCHEMA_VERSIONS = ("1.0", "2.0")
H5_NULL_ATTRIBUTE = "chemsmart_is_null"

#: Conditions at which the thermochemistry delegates are evaluated. A PySCF
#: Hessian carries no temperature of its own, so the standard state is the
#: only defensible default.
STANDARD_TEMPERATURE_K = 298.15
STANDARD_PRESSURE_ATM = 1.0


class PySCFOutput(FileMixin):
    """Typed reader for a completed (or failed) PySCF job.

    Constructed with the path to ``label.out`` so that ``FileMixin`` supplies
    ``contents``, ``content_lines_string`` and the file-identity helpers the
    database provenance needs. All numeric data is read from the sibling
    ``label.h5``.

    Attributes:
        filename (str): Path to the PySCF log file.
        resultsfile (str): Path to the structured results file.
    """

    def __init__(self, filename):
        self.filename = os.path.abspath(filename)
        stem, _ = os.path.splitext(self.filename)
        self.resultsfile = stem + ".h5"

    # ------------------------------------------------------------------
    # raw payload
    # ------------------------------------------------------------------

    @cached_property
    def _payload(self):
        """Return ``(spec, provenance, status, results)`` from the h5 file."""
        if not os.path.exists(self.resultsfile):
            raise FileNotFoundError(
                f"PySCF results file not found: {self.resultsfile}. "
                f"The calculation did not reach the point of writing results."
            )
        return read_pyscf_h5(self.resultsfile)

    @property
    def spec(self):
        """Settings as they were actually applied."""
        return self._payload[0]

    @property
    def provenance(self):
        """Versions, host, engine, timings and the settings digest."""
        return self._payload[1]

    @property
    def status(self):
        """Per-stage convergence, termination flag and any typed failure."""
        return self._payload[2]

    @property
    def results(self):
        """Numeric result arrays."""
        return self._payload[3]

    # ------------------------------------------------------------------
    # provenance / identity
    # ------------------------------------------------------------------

    @property
    def normal_termination(self):
        """Return whether every requested stage converged.

        Keyed on the structured status rather than a log string, so a run
        that crashed after writing a partial log cannot read as complete.
        """
        try:
            return bool(self.status.get("normal_termination", False))
        except (FileNotFoundError, ValueError, KeyError, OSError) as e:
            logger.debug(f"Could not read PySCF status: {e}")
            return False

    @property
    def failure(self):
        """Return the typed failure record, or None."""
        return self.status.get("failure")

    @property
    def version(self):
        """Return the PySCF version that produced this result."""
        return self.provenance.get("pyscf_version")

    @property
    def engine(self):
        """Return ``cpu`` or ``gpu``.

        Recorded because GPU and CPU results are not bit-identical, so a
        comparison set must never mix them.
        """
        return self.provenance.get("engine")

    @property
    def settings_digest(self):
        """Digest of the scientifically meaningful settings."""
        return self.provenance.get("settings_digest")

    @property
    def project_yaml_digest(self):
        """Digest of the source project YAML, or None when unavailable."""
        return self.provenance.get("project_yaml_digest")

    @property
    def file_date(self):
        """Return the results-file modification time.

        Overrides the Gaussian-flavoured default on ``FileMixin``, which
        scans the log tail for a Gaussian date banner.
        """
        if not os.path.exists(self.resultsfile):
            return None
        stamp = os.path.getmtime(self.resultsfile)
        return datetime.fromtimestamp(stamp).strftime("%Y-%m-%d %H:%M:%S")

    @property
    def total_elapsed_walltime(self):
        return self.provenance.get("wall_seconds")

    @property
    def total_core_hours(self):
        seconds = self.provenance.get("core_seconds")
        return seconds / 3600.0 if seconds is not None else None

    # ------------------------------------------------------------------
    # method / basis metadata
    # ------------------------------------------------------------------

    @property
    def method(self):
        return self.spec.get("method")

    @property
    def basis(self):
        return self.spec.get("basis")

    @property
    def jobtype(self):
        return self.spec.get("jobtype")

    @property
    def charge(self):
        return self.spec.get("charge")

    @property
    def multiplicity(self):
        return self.spec.get("multiplicity")

    @property
    def spin(self):
        """Return PySCF's ``mol.spin`` (2S), not the multiplicity."""
        return self.spec.get("spin")

    @property
    def num_basis_functions(self):
        return self.spec.get("num_basis_functions")

    @property
    def num_shells(self):
        return self.spec.get("num_shells")

    @property
    def point_group(self):
        # Version 1.0 placed this derived result in spec. Version 2.0 stores it
        # under results while preserving the public property for old artifacts.
        if "point_group" in self.spec:
            return self.spec.get("point_group")
        return self.results.get("point_group")

    @property
    def freq(self):
        """Return whether this run produced frequencies."""
        return "hess" in self.spec.get("stages", [])

    @property
    def rotational_symmetry_number(self):
        """Return the rotational symmetry number from the point group.

        Gaussian and ORCA print this; PySCF does not, so it is derived from
        the point group PySCF detected. It divides the rotational partition
        function, so getting it wrong shifts the rotational entropy by
        R*ln(sigma) -- about 1.4 kJ/mol at 298 K for sigma = 2.
        """
        return rotational_symmetry_number_from_point_group(self.point_group)

    @property
    def solvent_on(self):
        return self.spec.get("solvent_model") is not None

    @property
    def solvent_model(self):
        return self.spec.get("solvent_model")

    @property
    def solvent_id(self):
        return self.spec.get("solvent_id")

    @property
    def custom_solvent(self):
        """PySCF jobs never carry a custom solvent block."""
        return None

    @property
    def route_string(self):
        """Return a canonical one-line description of the calculation.

        There is no native route string for a library backend, so this is
        synthesised from the applied settings. It is stable for a given
        settings digest, which is what the database's route hashing needs.
        """
        spec = self.spec
        parts = [str(spec.get("method") or ""), str(spec.get("basis") or "")]
        if spec.get("density_fit"):
            parts.append("df")
        if spec.get("defgrid"):
            parts.append(str(spec["defgrid"]))
        if spec.get("dispersion"):
            parts.append(str(spec["dispersion"]))
        if spec.get("solvent_model"):
            parts.append(
                f"{spec['solvent_model']}({spec.get('solvent_id') or ''})"
            )
        parts.append(str(spec.get("jobtype") or ""))
        parts.append(f"engine={spec.get('engine')}")
        return " ".join(p for p in parts if p)

    # ------------------------------------------------------------------
    # energies and orbitals
    # ------------------------------------------------------------------

    @cached_property
    def energies(self):
        """Return every SCF energy in Hartree, in stage order."""
        values = self.results.get("energies")
        return [float(v) for v in values] if values is not None else []

    @cached_property
    def _eigenvalues(self):
        """Split MO energies into occupied/virtual, per spin channel.

        ``FileMixin`` derives HOMO, LUMO, SOMO and the FMO gaps from these
        four lists, so PySCF inherits exactly the definitions Gaussian and
        ORCA use instead of re-deriving them.

        **Converted from Hartree to eV.** Both ``Gaussian16Output`` and
        ``ORCAOutput`` document these lists as eV, and the database stores
        whatever the output object returns, so returning PySCF's native
        Hartree here would make orbital energies incomparable across
        programs by a factor of 27.2.
        """
        energy = self.results.get("mo_energy")
        occupancy = self.results.get("mo_occ")
        empty = ([], [], [], [])
        if energy is None or occupancy is None:
            return empty

        energy = np.asarray(energy)
        occupancy = np.asarray(occupancy)

        if energy.ndim == 2:  # unrestricted: shape (2, nmo)
            alpha_e, beta_e = energy[0], energy[1]
            alpha_o, beta_o = occupancy[0], occupancy[1]
        else:  # restricted: occupancies are 2.0 / 0.0
            alpha_e = beta_e = energy
            alpha_o = beta_o = occupancy

        ev = units.Hartree  # Hartree -> eV
        return (
            [float(e) * ev for e in alpha_e[alpha_o > 0]],
            [float(e) * ev for e in beta_e[beta_o > 0]],
            [float(e) * ev for e in alpha_e[alpha_o == 0]],
            [float(e) * ev for e in beta_e[beta_o == 0]],
        )

    @property
    def alpha_occ_eigenvalues(self):
        return self._eigenvalues[0]

    @property
    def beta_occ_eigenvalues(self):
        return self._eigenvalues[1]

    @property
    def alpha_virtual_eigenvalues(self):
        return self._eigenvalues[2]

    @property
    def beta_virtual_eigenvalues(self):
        return self._eigenvalues[3]

    # ------------------------------------------------------------------
    # structure and vibrations
    # ------------------------------------------------------------------

    @cached_property
    def positions(self):
        """Return the final geometry in Angstrom."""
        values = self.results.get("positions")
        return np.asarray(values, dtype=float) if values is not None else None

    @cached_property
    def forces(self):
        """Return final-geometry forces in Hartree/Bohr, or None."""
        values = self.results.get("forces")
        return np.asarray(values, dtype=float) if values is not None else None

    @cached_property
    def chemical_symbols(self):
        return list(self.spec.get("symbols", []))

    @cached_property
    def vibrational_frequencies(self):
        """Return harmonic frequencies in cm^-1, or None."""
        values = self.results.get("vibrational_frequencies")
        return [float(v) for v in values] if values is not None else None

    @cached_property
    def vibrational_modes(self):
        values = self.results.get("normal_modes")
        return np.asarray(values, dtype=float) if values is not None else None

    @cached_property
    def mulliken_atomic_charges(self):
        values = self.results.get("mulliken_charges")
        return [float(v) for v in values] if values is not None else None

    @cached_property
    def dipole_moment(self):
        values = self.results.get("dipole_moment")
        return [float(v) for v in values] if values is not None else None

    def get_molecule(self, index="-1", return_list=False):
        """Build a :class:`Molecule` from the stored geometry and results."""
        from chemsmart.io.molecules.structure import Molecule

        molecule = Molecule(
            symbols=self.chemical_symbols,
            positions=self.positions,
            charge=self.charge,
            multiplicity=self.multiplicity,
            energy=self.energies[-1] if self.energies else None,
            forces=self.forces,
            vibrational_frequencies=self.vibrational_frequencies,
            vibrational_modes=self.vibrational_modes,
            mulliken_atomic_charges=self.mulliken_atomic_charges,
            rotational_symmetry_number=self.rotational_symmetry_number,
            is_optimized_structure=self.jobtype in ("opt",),
        )
        if return_list:
            return [molecule]
        return molecule

    # ------------------------------------------------------------------
    # thermochemistry -- delegated, never recomputed
    # ------------------------------------------------------------------

    @cached_property
    def _thermochemistry(self):
        """Return ChemSmart's thermochemistry engine for this result.

        Evaluated at standard conditions. A Hessian is temperature
        independent -- the temperature enters only when the partition
        functions are evaluated -- and PySCF records no temperature, so
        these properties report the standard state. For any other condition,
        construct ``Thermochemistry`` directly (or use
        ``chemsmart run thermochemistry -T ...``), which is also how the
        Gaussian and ORCA paths handle non-standard conditions.

        Safe against recursion: ``Thermochemistry`` reads only ``jobtype``,
        ``normal_termination`` and the vibrational data from this object, and
        never touches a thermochemical property.
        """
        if not self.freq:
            return None
        from chemsmart.analysis.thermochemistry import Thermochemistry

        return Thermochemistry(
            filename=self.filename,
            temperature=STANDARD_TEMPERATURE_K,
            pressure=STANDARD_PRESSURE_ATM,
        )

    def _thermo(self, attribute):
        """Return a raw engine value in the engine's SI units."""
        engine = self._thermochemistry
        return getattr(engine, attribute, None) if engine else None

    def _thermo_energy(self, attribute):
        """Return an engine energy converted from J/mol to Hartree.

        The thermochemistry engine works in SI, but ``Gaussian16Output`` and
        ``ORCAOutput`` expose these quantities in **Hartree**, and the
        database stores whatever the output object returns. Converting here
        is what keeps a PySCF record numerically comparable with a Gaussian
        or ORCA record instead of silently off by a factor of 2.6e6.
        """
        value = self._thermo(attribute)
        return value * joule_per_mol_to_hartree if value is not None else None

    def _thermo_entropy(self, attribute):
        """Return an engine entropy converted from J/mol/K to Hartree/K.

        Matches ``ORCAOutput.entropy``, which is documented as Hartree/K.
        """
        value = self._thermo(attribute)
        return value * joule_per_mol_to_hartree if value is not None else None

    @property
    def temperature_in_K(self):
        return self._thermo("temperature")

    @property
    def pressure_in_atm(self):
        return self._thermo("pressure")

    @property
    def zero_point_energy(self):
        return self._thermo_energy("zero_point_energy")

    @property
    def internal_energy(self):
        """Return the absolute internal energy U.

        The engine's ``total_internal_energy`` is the thermal *correction*
        (it excludes the electronic energy), unlike its ``enthalpy`` and
        ``gibbs_free_energy``, which are absolute. Verified: electronic +
        correction reproduces H - RT.
        """
        correction = self._thermo_energy("total_internal_energy")
        electronic = self._thermo_energy("electronic_energy")
        if correction is None or electronic is None:
            return None
        return electronic + correction

    @property
    def enthalpy(self):
        return self._thermo_energy("enthalpy")

    @property
    def entropy(self):
        return self._thermo_entropy("total_entropy")

    @property
    def entropy_times_temperature(self):
        return self._thermo_energy("entropy_times_temperature")

    @property
    def gibbs_free_energy(self):
        return self._thermo_energy("gibbs_free_energy")

    @property
    def electronic_entropy(self):
        return self._thermo_entropy("electronic_entropy")

    @property
    def vibrational_entropy(self):
        return self._thermo_entropy("vibrational_entropy")

    @property
    def rotational_entropy(self):
        return self._thermo_entropy("rotational_entropy")

    @property
    def translational_entropy(self):
        return self._thermo_entropy("translational_entropy")

    @property
    def rotational_temperatures(self):
        return self._thermo("effective_rotational_temperatures")

    @property
    def rotational_constants_in_Hz(self):
        return self._thermo("effective_rotational_constants_in_Hz")

    # The "thermal correction to ..." quantities are Gaussian's naming. The
    # engine exposes the absolute terms, so the corrections are derived the
    # same way Gaussian defines them: the thermal quantity minus the
    # electronic energy. Kept in the engine's SI units, like its siblings.

    @property
    def thermal_vibration_correction(self):
        return self._thermo_energy("vibrational_internal_energy")

    @property
    def thermal_rotation_correction(self):
        return self._thermo_energy("rotational_internal_energy")

    @property
    def thermal_translation_correction(self):
        return self._thermo_energy("translational_internal_energy")

    def _correction_over_electronic(self, attribute):
        """Return (absolute quantity - electronic energy), in Hartree."""
        total = self._thermo_energy(attribute)
        electronic = self._thermo_energy("electronic_energy")
        if total is None or electronic is None:
            return None
        return total - electronic

    @property
    def thermal_energy_correction(self):
        # Already a correction in the engine -- see internal_energy.
        return self._thermo_energy("total_internal_energy")

    @property
    def thermal_enthalpy_correction(self):
        return self._correction_over_electronic("enthalpy")

    @property
    def thermal_gibbs_free_energy_correction(self):
        return self._correction_over_electronic("gibbs_free_energy")

    def __repr__(self):
        return (
            f"{type(self).__name__}(label={self.spec.get('label')!r}, "
            f"jobtype={self.jobtype!r}, engine={self.engine!r}, "
            f"normal_termination={self.normal_termination})"
        )


def read_pyscf_h5(filename):
    """Read a v1 or v2 PySCF HDF5 artifact.

    Version 1.0 is retained for Stage A artifacts and stores three mappings as
    JSON strings. Version 2.0 requires four top-level groups. Mapping arrays
    are converted to ordinary Python lists to preserve the historical public
    payload shape; values below ``results/`` remain exact numpy arrays.
    """
    import h5py

    with h5py.File(filename, "r") as handle:
        raw_version = handle.attrs.get("schema_version")
        version = _text(raw_version) if raw_version is not None else None
        if version not in SUPPORTED_SCHEMA_VERSIONS:
            raise ValueError(
                f"Unsupported PySCF results schema {version!r} in "
                f"{filename}; this ChemSmart understands "
                f"{SUPPORTED_SCHEMA_VERSIONS}."
            )

        if version == "1.0":
            _require_nodes(handle, ("spec", "provenance", "status", "results"))
            for name in ("spec", "provenance", "status"):
                if not isinstance(handle[name], h5py.Dataset):
                    raise ValueError(
                        f"PySCF schema 1.0 requires /{name} to be a JSON "
                        f"dataset in {filename}"
                    )
            if not isinstance(handle["results"], h5py.Group):
                raise ValueError(
                    f"PySCF schema 1.0 requires /results to be a group in "
                    f"{filename}"
                )
            spec = json.loads(_text(handle["spec"][()]))
            provenance = json.loads(_text(handle["provenance"][()]))
            status = json.loads(_text(handle["status"][()]))
            results = _read_group(handle["results"], preserve_arrays=True)
            return spec, provenance, status, results

        _require_nodes(handle, ("spec", "provenance", "status", "results"))
        for name in ("spec", "provenance", "status", "results"):
            if not isinstance(handle[name], h5py.Group):
                raise ValueError(
                    f"PySCF schema 2.0 requires /{name} to be a group in "
                    f"{filename}"
                )
        spec = _read_group(handle["spec"])
        provenance = _read_group(handle["provenance"])
        status = _read_group(handle["status"])
        results = _read_group(handle["results"], preserve_arrays=True)
        return spec, provenance, status, results


def _require_nodes(handle, names):
    """Raise a useful schema error for a missing top-level node."""
    missing = [name for name in names if name not in handle]
    if missing:
        raise ValueError(
            f"PySCF results artifact {handle.filename} is missing required "
            f"nodes: {missing}"
        )


def _read_group(group, preserve_arrays=False):
    """Decode a mapping group recursively."""
    import h5py

    decoded = {}
    for key in group:
        node = group[key]
        if isinstance(node, h5py.Group):
            decoded[key] = _read_group(node, preserve_arrays=preserve_arrays)
        else:
            decoded[key] = _read_dataset(node, preserve_arrays=preserve_arrays)
    return decoded


def _read_dataset(dataset, preserve_arrays=False):
    """Decode one typed dataset, including the explicit null marker."""
    if bool(dataset.attrs.get(H5_NULL_ATTRIBUTE, False)):
        return None

    value = dataset[()]
    if preserve_arrays:
        if isinstance(value, bytes):
            return value.decode("utf-8")
        return np.asarray(value)
    return _python_value(value)


def _python_value(value):
    """Convert an h5py scalar/array into JSON-compatible Python values."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.ndarray):
        return _python_value(value.tolist())
    if isinstance(value, np.generic):
        return _python_value(value.item())
    if isinstance(value, list):
        return [_python_value(item) for item in value]
    return value


def _text(value):
    """Decode an h5py string dataset to ``str``."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


#: Point groups whose rotational symmetry number is not derived from an
#: order n. PySCF spells the linear groups ``Coov`` and ``Dooh``.
_FIXED_SYMMETRY_NUMBERS = {
    "C1": 1,
    "CS": 1,
    "CI": 1,
    "COOV": 1,
    "DOOH": 2,
    "T": 12,
    "TD": 12,
    "TH": 12,
    "O": 24,
    "OH": 24,
    "I": 60,
    "IH": 60,
}


def rotational_symmetry_number_from_point_group(point_group):
    """Return the rotational symmetry number sigma for a point group.

    Standard result: sigma is the order of the rotational subgroup -- n for
    Cn/Cnv/Cnh, 2n for Dn/Dnh/Dnd, n/2 for Sn, and the fixed values above for
    the special groups.

    Returns ``1`` for an unrecognised or missing label, with a warning.
    That matches what Gaussian reports when symmetry is switched off, and is
    the conservative choice in the sense that it never *removes* entropy that
    a symmetric molecule really has -- but it does overestimate the entropy
    of a symmetric molecule, so the warning matters. Check it before
    comparing free energies across programs.
    """
    if not point_group:
        logger.warning(
            "No point group available; using rotational symmetry number 1. "
            "Rotational entropy may be overestimated."
        )
        return 1

    label = str(point_group).strip().upper()
    if label in _FIXED_SYMMETRY_NUMBERS:
        return _FIXED_SYMMETRY_NUMBERS[label]

    family = label[0]
    digits = "".join(ch for ch in label if ch.isdigit())
    if digits and family in ("C", "D", "S"):
        order = int(digits)
        if family == "C":
            return order
        if family == "D":
            return 2 * order
        if family == "S":
            # Sn only exists for even n; its rotational subgroup is C(n/2).
            return max(order // 2, 1)

    logger.warning(
        f"Unrecognised point group {point_group!r}; using rotational "
        f"symmetry number 1. Rotational entropy may be overestimated."
    )
    return 1
