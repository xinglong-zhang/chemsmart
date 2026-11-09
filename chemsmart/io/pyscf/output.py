"""Reader for PySCF calculation results.

This is a typed reader over the structured ``label.h5`` file that CHEMSMART's
own generated driver script writes -- not a log parser. ``io/orca/output.py``
is ~4,000 lines of regex because ORCA's printed text is the only interface it
offers; here CHEMSMART controls both ends, so the numbers come back as exact
float64 arrays and the sibling ``label.out`` (PySCF's own log) is never
parsed. It exists for humans and for program sniffing only.

Thermochemistry is **delegated** to ``chemsmart/analysis/thermochemistry.py``
rather than recomputed. PySCF ships ``pyscf.hessian.thermo.thermo()``, but
using it would give free energies that silently disagree with the Gaussian
and ORCA numbers for the same molecule, since CHEMSMART's engine has its own
rotational-symmetry and quasi-harmonic conventions. Cross-program consistency
is the whole point of a unified toolkit, so there is exactly one such engine.
"""

import logging
import os
from datetime import datetime
from functools import cached_property

import numpy as np
from ase import units

from chemsmart.io.pyscf.driver import APPLIED_SPEC_FIELDS, H5_NULL_ATTRIBUTE
from chemsmart.utils.constants import joule_per_mol_to_hartree
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.utils import parse_index_specification

logger = logging.getLogger(__name__)

#: Layouts this reader understands.
SUPPORTED_SCHEMA_VERSIONS = ("2.0",)

#: Conditions at which the thermochemistry delegates are evaluated. A PySCF
#: Hessian carries no temperature of its own, so the standard state is the
#: only defensible default.
STANDARD_TEMPERATURE_K = 298.15
STANDARD_PRESSURE_ATM = 1.0

#: Applied settings that do not change what a calculation computed, so two
#: results differing only in these are the same database record. Every other
#: field in ``APPLIED_SPEC_FIELDS`` reaches the record identity through
#: ``PySCFOutput.route_string``; a setting is therefore part of the identity
#: unless it is deliberately listed here.
#:
#: Only execution bookkeeping qualifies: the label, thread count and memory
#: change neither the numbers nor which of them were computed. The settings
#: digest is a hash of the other fields, so including it would count them
#: twice and re-key every result whenever the driver gains a setting.
#: Iteration limits, ``forces`` and ``stages`` stay in, as Gaussian's
#: ``maxcycles``, ``force`` and ``freq`` route keywords do: a run that stopped
#: at its limit, or stored fewer properties, must not replace a more
#: complete result.
ROUTE_IDENTITY_EXCLUDED = frozenset(
    {"label", "num_threads", "max_memory_mb", "settings_digest"}
)

#: Facts the driver derives that change the physical model, so they join
#: the record identity although no setting names them: the same basis name
#: with and without its core potential is a different calculation.
ROUTE_IDENTITY_DERIVED = ("ecp_core_electrons",)


class PySCFOutput(FileMixin):
    """Typed reader for a completed (or failed) PySCF job.

    Constructed with the path to ``label.out`` so that ``FileMixin`` supplies
    ``contents``, ``content_lines_string`` and the file-identity helpers the
    database provenance needs. All numeric data is read from the sibling
    ``label.h5``.

    Attributes:
        filename (str): Path to the PySCF log file.
        resultsfile (str): Path to the structured results file.
        h5file (PySCFH5File): Reader for that file.
    """

    def __init__(self, filename):
        # A log may be absent or empty after an early failure, so the results
        # file can be given directly; only its sibling path is derived here.
        self.filename = os.path.abspath(filename)
        stem, _ = os.path.splitext(self.filename)
        self.resultsfile = stem + ".h5"
        self.h5file = PySCFH5File(self.resultsfile)

    # ------------------------------------------------------------------
    # raw payload
    # ------------------------------------------------------------------

    @property
    def spec(self):
        """Settings as they were actually applied."""
        return self.h5file.spec

    @property
    def provenance(self):
        """Versions, host, engine, timings and the settings digest."""
        return self.h5file.provenance

    @property
    def status(self):
        """Per-stage convergence, termination flag and any typed failure."""
        return self.h5file.status

    @property
    def results(self):
        """Numeric result arrays."""
        return self.h5file.results

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
    def failure_reason(self):
        """Return why this run did not finish normally, or None if it did."""
        if self.normal_termination:
            return None
        try:
            status = self.status
        except (FileNotFoundError, ValueError) as e:
            return str(e)
        failure = status.get("failure")
        if failure:
            return (
                f"{failure.get('type')} at stage {failure.get('stage')!r}: "
                f"{failure.get('message')}"
            )
        unfinished = [
            (
                f"{stage} ({entry['not_a_transition_state']})"
                if "not_a_transition_state" in entry
                else stage
            )
            for stage, entry in status.get("stages", {}).items()
            if not entry.get("converged", entry.get("completed", False))
        ]
        return f"these stages did not finish: {', '.join(unfinished)}"

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
    def reference(self):
        """Return the SCF reference PySCF used, e.g. ``RKS`` or ``UHF``.

        None for results written before the reference was recorded; it is
        not reconstructed for them.
        """
        return self.spec.get("reference")

    @property
    def spin(self):
        """Return ``restricted`` or ``unrestricted``, as the Gaussian and ORCA
        readers do. PySCF's ``mol.spin`` (2S) is ``spec["spin"]``.

        Results written before the reference was recorded used
        ``scf.HF``/``dft.KS``, which choose a restricted reference exactly
        when 2S is zero.
        """
        reference = self.reference
        if reference is None:
            return (
                "restricted" if self.spec.get("spin") == 0 else "unrestricted"
            )
        return "unrestricted" if reference.startswith("U") else "restricted"

    @property
    def xc_applied(self):
        """Return the exchange-correlation functional PySCF evaluated.

        A mapping of backend component ids and coefficients, the short- and
        long-range exact-exchange fractions, the range-separation parameter
        and whether non-local correlation was included; None for Hartree-Fock
        and for results written before it was recorded.
        """
        return self.spec.get("xc_applied")

    @property
    def ecp_core_electrons(self):
        """Return core electrons replaced by a potential, by element.

        None for an all-electron calculation, and for results written before
        core potentials were attached (those ran every element all-electron).
        """
        return self.spec.get("ecp_core_electrons")

    @property
    def numerics_applied(self):
        """Return the numerical settings PySCF ran with, defaults resolved.

        SCF ``conv_tol`` and ``max_cycle``; for DFT the integration grid
        (level, per-atom radial/angular points, pruning, radial and Becke
        schemes, point count); and the density-fitting basis, or None.
        None for results written before it was recorded.
        """
        return self.spec.get("numerics_applied")

    @property
    def spin_square(self):
        """Return <S^2> of the final determinant, in units of hbar^2.

        For a correlated method this is the HF reference determinant.
        """
        return self._scalar_result("spin_square")

    @property
    def spin_contamination(self):
        """Return <S^2> - S(S+1) for the requested multiplicity, or None."""
        spin_square = self.spin_square
        if spin_square is None:
            return None
        total_spin = (self.multiplicity - 1) / 2
        return spin_square - total_spin * (total_spin + 1)

    @property
    def num_basis_functions(self):
        return self.spec.get("num_basis_functions")

    @property
    def num_shells(self):
        return self.spec.get("num_shells")

    @property
    def point_group(self):
        return self.results.get("point_group")

    @property
    def level_of_theory(self):
        """Return the applied level of theory, by ``PySCFJobSettings`` field."""
        spec = self.spec
        return {
            "functional": spec.get("xc"),
            "ab_initio": spec.get("ab_initio"),
            "basis": spec.get("basis"),
            "density_fit": spec.get("density_fit"),
            "aux_basis": spec.get("aux_basis"),
            "defgrid": spec.get("defgrid"),
            "dispersion": spec.get("dispersion"),
            "solvent_model": spec.get("solvent_model"),
            "solvent_id": spec.get("solvent_id"),
        }

    @property
    def freq(self):
        """Return whether this run produced frequencies."""
        return self.results.get("vibrational_frequencies") is not None

    @property
    def rotational_symmetry_number(self):
        """Return the rotational symmetry number from the point group.

        Gaussian and ORCA print this; PySCF does not, so it is derived from
        the point group PySCF detected. It divides the rotational partition
        function, so getting it wrong shifts the rotational entropy by
        R*ln(sigma) -- about 1.7 kJ/mol at 298 K for sigma = 2.
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
        synthesised from the applied settings. The readable part names the
        level of theory; the trailing ``spec=`` hash covers every applied
        setting outside ``ROUTE_IDENTITY_EXCLUDED``, so a setting that changes
        the result cannot be left out of the database record identity. Unset
        (``None``) fields are omitted, so a field added to a later driver does
        not change the identity of results written before it existed.
        """
        # Deferred: the database package imports this reader.
        from chemsmart.database.utils import canonical_json_hash

        spec = self.spec
        parts = [str(spec.get("method") or ""), str(spec.get("basis") or "")]
        if spec.get("density_fit"):
            parts.append("df")
        # These change the numbers but have no structured database field, so
        # they must reach the route hash or the records would overwrite.
        for name in ("aux_basis", "scf_tol", "opt_solver"):
            if spec.get(name) is not None:
                parts.append(f"{name}={spec[name]}")
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
        identity = {
            name: spec[name]
            for name in APPLIED_SPEC_FIELDS + ROUTE_IDENTITY_DERIVED
            if name not in ROUTE_IDENTITY_EXCLUDED
            and spec.get(name) is not None
        }
        parts.append(f"spec={canonical_json_hash(identity)}")
        return " ".join(p for p in parts if p)

    # ------------------------------------------------------------------
    # energies and orbitals
    # ------------------------------------------------------------------

    @cached_property
    def energies(self):
        """Return the total energies in Hartree at the requested level.

        As in the Gaussian reader: the correlated total for MP2, CCSD or
        CCSD(T), which is empty when its iterations did not converge, and
        the SCF energies, in stage order, otherwise.
        """
        values = self.results.get("energies")
        return [float(v) for v in values] if values is not None else []

    @property
    def scf_energies(self):
        """Return the SCF (reference) energies in Hartree, in stage order."""
        values = self.results.get("scf_energies")
        if values is None:
            # Written before correlated methods, when they were the same.
            return self.energies
        return [float(v) for v in values]

    @property
    def correlation_energy(self):
        """Return the correlation energy in Hartree, or None.

        For CCSD(T) it includes the perturbative triples correction.
        """
        return self._scalar_result("correlation_energy")

    @property
    def ccsd_correlation_energy(self):
        """Return the CCSD part of a CCSD(T) correlation energy, or None."""
        return self._scalar_result("ccsd_correlation_energy")

    @property
    def triples_correction(self):
        """Return the perturbative (T) correction in Hartree, or None."""
        return self._scalar_result("triples_correction")

    @property
    def saddle_order(self):
        """Return the number of imaginary modes found by a TS search.

        One for a transition state. None for other jobs, and for a search
        whose result had no Hessian, so it is not known to be a saddle point.
        """
        return self.status.get("stages", {}).get("ts", {}).get("saddle_order")

    @property
    def scf_stability(self):
        """Return the stability analysis of the converged SCF reference.

        A mapping with the reference and, per question, whether it is stable
        and its lowest Hessian eigenvalue in Eh: ``internal``,
        ``real_to_complex`` and ``external`` (its ``question`` names the
        direction, e.g. RHF/RKS -> UHF/UKS; ROHF has neither of the last two).
        A ``stable`` of None means PySCF gave no answer. None when stability
        was not requested, never an implicit "stable".
        """
        record = self.results.get("scf_stability")
        return None if record is None else PySCFH5File.python_value(record)

    def _scalar_result(self, name):
        value = self.results.get(name)
        return None if value is None else float(value)

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

    @property
    def optimizer_converged(self):
        """Return whether a geometry optimisation met its criteria.

        For a transition-state search this is the saddle-point search's own
        convergence; whether the result is a transition state is
        :attr:`saddle_order`.
        """
        stages = self.status.get("stages", {})
        stage = stages.get("opt", stages.get("ts", {}))
        return bool(stage.get("optimizer_converged", False))

    @cached_property
    def all_structures(self):
        """Return the stored structures: a PySCF job keeps only the final one.

        Raises:
            ValueError: If the run failed before producing a geometry.
        """
        from chemsmart.io.molecules.structure import Molecule

        if self.positions is None:
            failure = self.failure or {}
            raise ValueError(
                f"No structure in {self.resultsfile}; the calculation failed "
                f"at stage {failure.get('stage')!r}: "
                f"{failure.get('type')}: {failure.get('message')}"
            )
        dipole = self.dipole_moment
        point_group = self.point_group
        return [
            Molecule(
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
                is_optimized_structure=self.optimizer_converged,
                dipole_moment=(
                    np.asarray(dipole) if dipole is not None else None
                ),
                dipole_moment_magnitude=(
                    float(np.linalg.norm(dipole))
                    if dipole is not None
                    else None
                ),
                point_group=(
                    str(point_group).upper() if point_group else None
                ),
            )
        ]

    @property
    def optimized_structure(self):
        """Return the optimised structure, or None if not converged."""
        if not self.optimizer_converged:
            return None
        return self.all_structures[-1]

    def get_molecule(self, index="-1", return_list=False):
        """Return the stored structures selected by a 1-based index."""
        selector = parse_index_specification(str(index))
        structures = self.all_structures
        if isinstance(selector, list):
            return [structures[i] for i in selector]
        selected = structures[selector]
        if return_list and not isinstance(selected, list):
            return [selected]
        return selected

    # ------------------------------------------------------------------
    # thermochemistry -- delegated, never recomputed
    # ------------------------------------------------------------------

    @cached_property
    def _thermochemistry(self):
        """Return CHEMSMART's thermochemistry engine for this result.

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
        if not self.freq or not self.normal_termination:
            return None
        from chemsmart.analysis.thermochemistry import Thermochemistry

        # The results file always exists here, whereas the log may not.
        return Thermochemistry(
            filename=self.resultsfile,
            temperature=STANDARD_TEMPERATURE_K,
            pressure=STANDARD_PRESSURE_ATM,
        )

    def _thermo(self, attribute):
        """Return a raw engine value in the engine's SI units."""
        engine = self._thermochemistry
        return getattr(engine, attribute) if engine else None

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
    # electronic energy.

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


class PySCFH5File(FileMixin):
    """Structured results file written by CHEMSMART's PySCF driver.

    The file is read on first access and its layout is checked before any
    value is returned. Mapping groups decode to ordinary Python values;
    values below ``results/`` remain exact numpy arrays.

    Attributes:
        filename (str): Path to the ``label.h5`` file.
    """

    #: Top-level groups every supported layout carries.
    SECTIONS = ("spec", "provenance", "status", "results")

    def __init__(self, filename):
        self.filename = filename

    @property
    def spec(self):
        """Settings as they were actually applied."""
        return self._sections["spec"]

    @property
    def provenance(self):
        """Versions, host, engine, timings and the settings digest."""
        return self._sections["provenance"]

    @property
    def status(self):
        """Per-stage convergence, termination flag and any typed failure."""
        return self._sections["status"]

    @property
    def results(self):
        """Numeric result arrays."""
        return self._sections["results"]

    @cached_property
    def _sections(self):
        """Return every top-level group, decoded, keyed by name."""
        if not os.path.exists(self.filename):
            raise FileNotFoundError(
                f"PySCF results file not found: {self.filename}. "
                f"The calculation did not reach the point of writing results."
            )
        try:
            return self._read()
        except OSError as e:
            raise ValueError(
                f"{self.filename} is not a readable PySCF results file: {e}"
            ) from e

    def _read(self):
        """Check the layout and decode every section in one pass."""
        import h5py

        with h5py.File(self.filename, "r") as handle:
            raw_version = handle.attrs.get("schema_version")
            version = (
                self._text(raw_version) if raw_version is not None else None
            )
            if version not in SUPPORTED_SCHEMA_VERSIONS:
                raise ValueError(
                    f"Unsupported PySCF results schema {version!r} in "
                    f"{self.filename}; this CHEMSMART understands "
                    f"{SUPPORTED_SCHEMA_VERSIONS}."
                )

            missing = [name for name in self.SECTIONS if name not in handle]
            if missing:
                raise ValueError(
                    f"PySCF results artifact {self.filename} is missing "
                    f"required nodes: {missing}"
                )
            for name in self.SECTIONS:
                if not isinstance(handle[name], h5py.Group):
                    raise ValueError(
                        f"PySCF schema {version} requires /{name} to be a "
                        f"group in {self.filename}"
                    )
            # Another driver writes the same schema version with its own
            # result contract; its energies sit under other names and would
            # be misread, e.g. an MP2 record reporting its HF reference
            # energy.
            if "result_contract_version" in handle["spec"]:
                contract = self._text(
                    handle["spec/result_contract_version"][()]
                )
                raise ValueError(
                    f"{self.filename} declares the PySCF result contract "
                    f"{contract!r}, a layout this CHEMSMART does not read. "
                    "Its quantities are stored under other names and would "
                    "be misread; rerun the calculation with this CHEMSMART."
                )
            return {
                name: self._read_group(
                    handle[name], preserve_arrays=name == "results"
                )
                for name in self.SECTIONS
            }

    @classmethod
    def _read_group(cls, group, preserve_arrays=False):
        """Decode a mapping group recursively."""
        import h5py

        decoded = {}
        for key in group:
            node = group[key]
            if isinstance(node, h5py.Group):
                decoded[key] = cls._read_group(
                    node, preserve_arrays=preserve_arrays
                )
            else:
                decoded[key] = cls._read_dataset(
                    node, preserve_arrays=preserve_arrays
                )
        return decoded

    @classmethod
    def _read_dataset(cls, dataset, preserve_arrays=False):
        """Decode one typed dataset, including the explicit null marker."""
        if bool(dataset.attrs.get(H5_NULL_ATTRIBUTE, False)):
            return None

        value = dataset[()]
        if preserve_arrays:
            if isinstance(value, bytes):
                return value.decode("utf-8")
            return np.asarray(value)
        return cls.python_value(value)

    @classmethod
    def python_value(cls, value):
        """Convert decoded h5py values, mappings included, to plain Python."""
        if isinstance(value, dict):
            return {key: cls.python_value(item) for key, item in value.items()}
        if isinstance(value, bytes):
            return value.decode("utf-8")
        if isinstance(value, np.ndarray):
            return cls.python_value(value.tolist())
        if isinstance(value, np.generic):
            return cls.python_value(value.item())
        if isinstance(value, list):
            return [cls.python_value(item) for item in value]
        return value

    @staticmethod
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
