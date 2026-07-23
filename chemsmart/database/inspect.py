"""
Database inspection module for viewing database overview, record details,
and structure details.
"""

import logging
import os
import sqlite3

from chemsmart.database.database import Database
from chemsmart.database.utils import (
    bool_to_str,
    format_energy,
    format_float,
    format_kv,
    human_size,
    open_connection,
    separator,
    sort_structure_dicts_by_energy,
    truncate_iso,
)
from chemsmart.utils.constants import energy_conversion

logger = logging.getLogger(__name__)


class DatabaseInspector:
    """Inspect a chemsmart database.

    Provides five views:
    Overview – database-level metadata and statistics
    Record detail – full detail for one record
    Structure detail – full detail for one structure within a record
    Molecule detail – full detail for a molecule (chemical species)
    Standalone structure detail – full detail for a structure (conformer)
    """

    def __init__(
        self,
        db_file,
        index=None,
        record_id=None,
        structure_index=None,
        molecule_id=None,
        structure_id=None,
    ):
        self.db_file = db_file
        self.db = Database(db_file)
        self.index = index
        self.record_id = record_id
        self.structure_index = structure_index
        self.molecule_id = molecule_id
        self.structure_id = structure_id

    def resolve_id(self):
        """Return the full record_id from either index or partial ID."""
        if self.index is not None:
            record = self.db.get_record(record_index=self.index)
            if record is None:
                raise ValueError(f"No record found with index {self.index}.")
            return record["record_id"]
        if self.record_id is not None:
            return self.db.get_record_by_partial_id(self.record_id)
        raise ValueError("Either index or record_id must be provided.")

    def overview(self):
        """Return database-level statistics as a dictionary."""
        conn = open_connection(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            stats = {}

            # File info
            stats["db_file"] = os.path.basename(self.db_file)
            stats["db_size"] = os.path.getsize(self.db_file)

            # Record counts
            stats["num_records"] = conn.execute(
                "SELECT COUNT(*) FROM records"
            ).fetchone()[0]
            stats["num_molecules"] = conn.execute(
                "SELECT COUNT(*) FROM molecules"
            ).fetchone()[0]
            stats["num_structures"] = conn.execute(
                "SELECT COUNT(*) FROM structures"
            ).fetchone()[0]
            stats["num_record_structures"] = conn.execute(
                "SELECT COUNT(*) FROM record_structures"
            ).fetchone()[0]

            # Program breakdown
            rows = conn.execute(
                "SELECT program, COUNT(*) AS cnt FROM records GROUP BY program ORDER BY cnt DESC"
            ).fetchall()
            stats["programs"] = [(r["program"], r["cnt"]) for r in rows]

            # Job type breakdown
            rows = conn.execute(
                "SELECT jobtype, COUNT(*) AS cnt FROM records GROUP BY jobtype ORDER BY cnt DESC"
            ).fetchall()
            stats["jobtypes"] = [(r["jobtype"], r["cnt"]) for r in rows]

            # Unique methods and basis sets
            rows = conn.execute(
                "SELECT method, COUNT(*) AS cnt FROM records WHERE method IS NOT NULL GROUP BY method ORDER BY cnt DESC"
            ).fetchall()
            stats["methods"] = [(r["method"], r["cnt"]) for r in rows]
            rows = conn.execute(
                "SELECT basis, COUNT(*) AS cnt FROM records WHERE basis IS NOT NULL GROUP BY basis ORDER BY cnt DESC"
            ).fetchall()
            stats["basis_sets"] = [(r["basis"], r["cnt"]) for r in rows]

            # Energy, thermochemistry, and solvation counts (from records table)
            row = conn.execute("""
                SELECT
                    COUNT(CASE WHEN total_energy IS NOT NULL THEN 1 END) AS energy_count,
                    COUNT(CASE WHEN gibbs_free_energy IS NOT NULL THEN 1 END) AS thermo_count,
                    COUNT(CASE WHEN solvent_on = 1 THEN 1 END) AS solvated_records,
                    COUNT(CASE WHEN solvent_on = 0 THEN 1 END) AS gas_phase_records
                FROM records
                """).fetchone()
            stats["energy_count"] = row["energy_count"]
            stats["thermo_count"] = row["thermo_count"]
            stats["solvated_records"] = row["solvated_records"]
            stats["gas_phase_records"] = row["gas_phase_records"]

            # Forces count (records that have at least one structure with forces)
            stats["forces_count"] = conn.execute(
                "SELECT COUNT(DISTINCT record_id) FROM record_structures WHERE forces_json IS NOT NULL"
            ).fetchone()[0]
            # Solvent count
            rows = conn.execute("""
                SELECT solvent_id, COUNT(*) AS cnt
                FROM records
                WHERE solvent_on = 1
                GROUP BY solvent_id
                ORDER BY cnt DESC
                """).fetchall()
            stats["solvents"] = [(r["solvent_id"], r["cnt"]) for r in rows]

            # Assembled-at date range
            row = conn.execute(
                "SELECT MIN(assembled_at) AS first, MAX(assembled_at) AS last FROM records"
            ).fetchone()
            stats["assembled_first"] = row["first"]
            stats["assembled_last"] = row["last"]

            return stats
        finally:
            conn.close()

    def record_detail(self):
        """Return full record data with structure summaries."""
        record_id = self.resolve_id()
        record = self.db.get_record(record_id=record_id)
        if record is None:
            raise ValueError(f"Record not found: {record_id}")
        return record

    def structure_detail(self):
        """Return full structure data for one structure in a record."""
        record = self.record_detail()
        molecules = record.get("molecules", [])
        for s in molecules:
            if s.get("index") == self.structure_index:
                return record, s
        raise ValueError(
            f"Structure index {self.structure_index} not found in record "
            f"(available: {[s.get('index') for s in molecules]})."
        )

    def format_overview(self):
        """Return a human-readable overview string."""
        stats = self.overview()
        lines = []
        lines.append("")
        lines.append(separator("Database Overview"))

        # General
        lines.append(format_kv("File", stats["db_file"]))
        lines.append(format_kv("Size", human_size(stats["db_size"])))
        lines.append(format_kv("Records", stats["num_records"]))
        lines.append(format_kv("Molecules", stats["num_molecules"]))
        lines.append(
            format_kv("Structures (Conformers)", stats["num_structures"])
        )
        assembled_first = truncate_iso(stats.get("assembled_first"))
        assembled_last = truncate_iso(stats.get("assembled_last"))
        lines.append(format_kv("Created", f"{assembled_first} (UTC)"))
        lines.append(format_kv("Last Updated", f"{assembled_last} (UTC)"))

        # Programs
        lines.append("")
        lines.append(separator("Programs"))
        for prog, cnt in stats["programs"]:
            lines.append(f"  {prog or '(unknown)':<30}: {cnt}")

        # Functionals & Basis sets
        lines.append("")
        lines.append(separator("Electronic Structure Methods"))
        if stats["methods"]:
            methods_str = ", ".join(
                f"{f[0]} ({f[1]})" for f in stats["methods"]
            )
            lines.append(format_kv("Methods", methods_str))
        else:
            lines.append(format_kv("Methods", "NULL"))
        if stats["basis_sets"]:
            basis_sets_str = ", ".join(
                f"{b[0]} ({b[1]})" for b in stats["basis_sets"]
            )
            lines.append(format_kv("Basis Sets", basis_sets_str))
        else:
            lines.append(format_kv("Basis Sets", "NULL"))

        # Solvation
        lines.append("")
        lines.append(separator("Solvation"))
        lines.append(
            format_kv("Gas Phase Records", stats["gas_phase_records"])
        )
        lines.append(format_kv("Solvated Records", stats["solvated_records"]))
        if stats["solvents"]:
            solvents_str = ", ".join(
                f"{s[0]} ({s[1]})" for s in stats["solvents"]
            )
            lines.append(format_kv("Solvents", solvents_str))

        # Job types
        lines.append("")
        lines.append(separator("Job Types"))
        for jt, cnt in stats["jobtypes"]:
            lines.append(f"  {jt or '(unknown)':<30}: {cnt}")

        # Calculation data
        lines.append("")
        lines.append(separator("Calculation Data"))
        lines.append(
            format_kv("Records with Energies", f"{stats['energy_count']}")
        )
        lines.append(
            format_kv("Records with Forces", f"{stats['forces_count']}")
        )
        lines.append(
            format_kv("Records with Frequencies", f"{stats['thermo_count']}")
        )

        lines.append("")
        lines.append(separator())
        lines.append("")
        return "\n".join(lines)

    def format_record_detail(self):
        """Return a human-readable record detail string."""
        record = self.record_detail()
        meta = record.get("meta", {})
        results = record.get("results", {})
        provenance = record.get("provenance", {})
        molecules = record.get("molecules", [])

        lines = []
        lines.append("")
        lines.append(separator("Record Detail"))

        # Identification
        lines.append(format_kv("Record Index", record.get("record_index")))
        lines.append(format_kv("Record ID", record.get("record_id")))
        lines.append(
            format_kv("Source", os.path.basename(provenance.get("source")))
        )

        # Meta
        lines.append("")
        lines.append(separator("Calculation"))
        lines.append(format_kv("Program", provenance.get("program")))
        lines.append(
            format_kv("Program Version", provenance.get("program_version"))
        )
        lines.append(format_kv("Method", meta.get("method")))
        lines.append(format_kv("Basis Set", meta.get("basis")))
        custom_basis = meta.get("custom_basis") or {}
        custom_solvent = meta.get("custom_solvent") or {}
        if meta.get("basis") == "customized_basis":
            has_detail = any(
                custom_basis.get(k) is not None
                for k in (
                    "light_elements",
                    "light_elements_basis",
                    "heavy_elements",
                    "heavy_elements_basis",
                )
            )
            if not has_detail:
                lines.append(
                    format_kv(
                        "  Basis Detail",
                        "unavailable (use #p in route to print expanded Gen/GenECP basis)",
                    )
                )
        if custom_basis:
            light_elems = custom_basis.get("light_elements")
            light_basis = custom_basis.get("light_elements_basis")
            heavy_elems = custom_basis.get("heavy_elements")
            heavy_basis = custom_basis.get("heavy_elements_basis")
            heavy_ecp = custom_basis.get("heavy_elements_ecp")
            if light_elems:
                light_str = ", ".join(light_elems)
                if light_basis:
                    light_str += f"  [{light_basis}]"
                lines.append(format_kv("  Light Elements", light_str))
            if heavy_elems:
                lines.append(
                    format_kv("  Heavy Elements", ", ".join(heavy_elems))
                )
            if heavy_basis:
                for elem, shells in heavy_basis.items():
                    shell_types = ", ".join(s["shell"] for s in shells)
                    lines.append(format_kv(f"    {elem} Shells", shell_types))
            if heavy_ecp:
                ecp_elems = ", ".join(heavy_ecp.keys())
                lines.append(format_kv("  Heavy Elements ECP", ecp_elems))
        lines.append(format_kv("Spin Type", meta.get("spin")))
        lines.append(format_kv("Job Type", meta.get("jobtype")))
        lines.append(format_kv("Solvent", bool_to_str(meta.get("solvent_on"))))
        if meta.get("solvent_on"):
            lines.append(format_kv("Solvent Model", meta.get("solvent_model")))
            lines.append(format_kv("Solvent ID", meta.get("solvent_id")))
            if custom_solvent:
                solvent_name = custom_solvent.get("SolventName")
                if solvent_name and solvent_name.strip().lower() != "generic":
                    lines.append(
                        format_kv("  Custom Solvent Name", solvent_name)
                    )
                for key in (
                    "Eps",
                    "EpsInf",
                    "HbondAcidity",
                    "HbondBasicity",
                    "SurfaceTensionAtInterface",
                    "CarbonAromaticity",
                    "ElectronegativeHalogenicity",
                ):
                    val = custom_solvent.get(key)
                    if val is not None:
                        lines.append(format_kv(f"  {key}", val))
        lines.append(format_kv("Route String", meta.get("route_string")))

        # Results — electronic
        lines.append("")
        lines.append(separator("Electronic Results"))
        lines.append(
            format_kv(
                "Total Energy (Eh)", format_energy(results.get("total_energy"))
            )
        )
        lines.append(
            format_kv(
                "Unpaired Electrons", results.get("num_unpaired_electrons")
            )
        )
        # Closed-shell orbitals
        if results.get("num_unpaired_electrons") == 0:
            lines.append(
                format_kv(
                    "HOMO Energy (eV)",
                    format_energy(results.get("homo_energy")),
                )
            )
            lines.append(
                format_kv(
                    "LUMO Energy (eV)",
                    format_energy(results.get("lumo_energy")),
                )
            )
        else:
            lines.append(
                format_kv(
                    "α-HOMO Energy (eV)",
                    format_energy(results.get("alpha_homo_energy")),
                )
            )
            lines.append(
                format_kv(
                    "β-HOMO Energy (eV)",
                    format_energy(results.get("beta_homo_energy")),
                )
            )
            lines.append(
                format_kv(
                    "α-LUMO Energy (eV)",
                    format_energy(results.get("alpha_lumo_energy")),
                )
            )
            lines.append(
                format_kv(
                    "β-LUMO Energy (eV)",
                    format_energy(results.get("beta_lumo_energy")),
                )
            )
            lines.append(
                format_kv(
                    "α-FMO Gap (eV)",
                    format_energy(results.get("alpha_fmo_gap")),
                )
            )
            lines.append(
                format_kv(
                    "β-FMO Gap (eV)",
                    format_energy(results.get("beta_fmo_gap")),
                )
            )
            if results.get("somo_energies") is not None:
                somo = results["somo_energies"]
                lines.append(format_kv("SOMO Energies (eV)", f"{somo}"))
        lines.append(
            format_kv("FMO Gap (eV)", format_energy(results.get("fmo_gap")))
        )

        # Performance
        if results.get("total_core_hours") is not None:
            lines.append("")
            lines.append(separator("Performance"))
            lines.append(
                format_kv(
                    "Core Hours",
                    format_float(results.get("total_core_hours"), 2),
                )
            )
            lines.append(
                format_kv(
                    "Elapsed Walltime (h)",
                    format_float(results.get("total_elapsed_walltime"), 2),
                )
            )

        # Thermochemistry (only if present)
        if results.get("gibbs_free_energy") is not None:
            lines.append("")
            lines.append(separator("Thermochemistry"))
            lines.append(
                format_kv("Temperature (K)", f"{meta.get('temperature_in_K')}")
            )
            lines.append(
                format_kv("Pressure (atm)", f"{meta.get('pressure_in_atm')}")
            )
            lines.append(
                format_kv(
                    "Zero Point Energy (Eh)",
                    format_energy(results.get("zero_point_energy")),
                )
            )
            lines.append(
                format_kv(
                    "Internal Energy (Eh)",
                    format_energy(results.get("internal_energy")),
                )
            )
            lines.append(
                format_kv(
                    "Enthalpy (Eh)", format_energy(results.get("enthalpy"))
                )
            )
            lines.append(
                format_kv(
                    "Entropy (Eh/K)", format_energy(results.get("entropy"))
                )
            )
            lines.append(
                format_kv(
                    "Gibbs Free Energy (Eh)",
                    format_energy(results.get("gibbs_free_energy")),
                )
            )

        # Provenance
        lines.append("")
        lines.append(separator("Provenance"))
        lines.append(format_kv("Parser", provenance.get("parser")))
        lines.append(
            format_kv("CHEMSMART Version", provenance.get("chemsmart_version"))
        )
        lines.append(format_kv("Source", provenance.get("source")))
        lines.append(
            format_kv("Source Hash", provenance.get("source_file_hash"))
        )
        lines.append(
            format_kv(
                "Source Size",
                human_size(provenance.get("source_size")),
            )
        )
        lines.append(format_kv("Source Date", provenance.get("source_date")))
        lines.append(format_kv("Assembled At", provenance.get("assembled_at")))
        if provenance.get("normal_termination") is not None:
            lines.append(
                format_kv(
                    "Normal Termination",
                    bool_to_str(provenance.get("normal_termination")),
                )
            )

        # Structure summary table
        if molecules:
            lines.append("")
            lines.append(separator(f"Structures ({len(molecules)})"))

            # Check whether identity fields are uniform across all molecules
            formulas = {(s.get("chemical_formula") or "") for s in molecules}
            charges = {s.get("charge") for s in molecules}
            mults = {s.get("multiplicity") for s in molecules}
            natoms = {s.get("number_of_atoms") for s in molecules}
            uniform = (
                len(formulas) == 1
                and len(charges) == 1
                and len(mults) == 1
                and len(natoms) == 1
            )

            if uniform:
                # Print shared identity once, then show differentiating columns
                formula = next(iter(formulas))
                charge = next(iter(charges))
                mult = next(iter(mults))
                lines.append(format_kv("Formula", formula or "N/A"))
                lines.append(format_kv("Charge", charge))
                lines.append(format_kv("Multiplicity", mult))
                lines.append("")
                lines.append(
                    f"  {'Idx':>4}  {'Structure ID':<12}  {'Energy (Eh)':>15}  {'Optimized':>9}"
                )
                lines.append(
                    f"  {'----':>4}  {'------------':<12}  {'---------------':>15}  {'---------':>9}"
                )
                for s in molecules:
                    idx = s.get("index", "-")
                    sid = str(s.get("structure_id", ""))[:12]
                    energy = s.get("energy")
                    energy_str = (
                        format_energy(energy) if energy is not None else "N/A"
                    )
                    is_opt = s.get("is_optimized_structure")
                    opt_str = bool_to_str(is_opt) if is_opt is not None else ""
                    lines.append(
                        f"  {idx:>4}  {sid:<12}  {energy_str:>15}  {opt_str:>9}"
                    )
            else:
                # Heterogeneous structures — show full info
                hdr = (
                    f"  {'Idx':>4}  {'Structure ID':<12}  {'Formula':<20}  {'Charge':>6}  "
                    f"{'Mult':>4}  {'Energy (Eh)':>15}"
                )
                lines.append(hdr)
                lines.append(
                    f"  {'----':>4}  {'------------':<12}  {'-------':<20}  {'------':>6}  "
                    f"{'----':>4}  {'---------------':>15}"
                )
                for s in molecules:
                    idx = s.get("index", "-")
                    sid = str(s.get("structure_id", ""))[:12]
                    formula = s.get("chemical_formula", "") or ""
                    charge = s.get("charge", "")
                    mult = s.get("multiplicity", "")
                    energy = s.get("energy")
                    energy_str = (
                        format_energy(energy) if energy is not None else ""
                    )
                    lines.append(
                        f"  {idx:>4}  {sid:<12}  {formula:<20}  {charge:>6}  "
                        f"{mult:>4}  {energy_str:>15}"
                    )

        db_name = os.path.basename(self.db_file)
        rid = str(record.get("record_id"))[:12]
        lines.append("")
        lines.append(
            f"  Tip: chemsmart run database inspect -f {db_name}"
            f" --rid {rid} --si <structure index>"
        )

        lines.append("")
        lines.append(separator())
        lines.append("")
        return "\n".join(lines)

    def format_structure_detail(self):
        """Return a human-readable structure detail string."""
        record, struct = self.structure_detail()
        provenance = record.get("provenance", {})
        lines = []
        lines.append("")
        lines.append(separator("Record Context"))

        # Context
        lines.append(format_kv("Record Index", record.get("record_index")))
        lines.append(format_kv("Record ID", record.get("record_id")))
        lines.append(
            format_kv("Source", os.path.basename(provenance.get("source")))
        )
        lines.append(
            format_kv(
                "Structure Index in File",
                struct.get("structure_index_in_file"),
            )
        )
        lines.append(
            format_kv("Energy (Eh)", format_energy(struct.get("energy")))
        )
        is_opt = struct.get("is_optimized_structure")
        if is_opt is not None:
            lines.append(format_kv("Optimized", bool_to_str(is_opt)))

        # Structure identity fields
        lines.append("")
        lines.append(separator("Structure Detail"))
        lines.append(format_kv("Structure ID", struct.get("structure_id")))
        lines.append(format_kv("Charge", struct.get("charge")))
        lines.append(format_kv("Multiplicity", struct.get("multiplicity")))

        # Molecule-level identity fields
        lines.append("")
        lines.append(separator("Parent Molecule"))
        lines.append(format_kv("Molecule ID", struct.get("molecule_id")))
        lines.append(
            format_kv("Chemical Formula", struct.get("chemical_formula"))
        )
        lines.append(
            format_kv("Mass (amu)", format_float(struct.get("mass"), 4))
        )
        lines.append(format_kv("SMILES", struct.get("smiles")))

        # Coordinates
        symbols = struct.get("chemical_symbols", [])
        positions = struct.get("positions", [])
        if symbols and positions:
            lines.append("")
            lines.append(separator(f"Coordinates ({len(symbols)} atoms, Å)"))
            lines.append(
                f"  {'Idx':>4}  {'Elem':<4}  "
                f"{'X':>14}  {'Y':>14}  {'Z':>14}"
            )
            lines.append(
                f"  {'----':>4}  {'----':<4}  "
                f"{'----------':>14}  {'----------':>14}  {'----------':>14}"
            )
            for i, (sym, pos) in enumerate(zip(symbols, positions)):
                x, y, z = pos
                lines.append(
                    f"  {i + 1:>4}  {sym:<4}  "
                    f"{x:>14.8f}  {y:>14.8f}  {z:>14.8f}"
                )

        # Frozen atoms
        frozen = struct.get("frozen_atoms")
        if frozen:
            lines.append("")
            lines.append(separator("Frozen Atoms"))
            lines.append(f"  {frozen}")

        # Physical properties
        com = struct.get("center_of_mass")
        moi = struct.get("moments_of_inertia")
        rot_sym = struct.get("rotational_symmetry_number")
        rot_consts = struct.get("rotational_constants")
        point_group = struct.get("point_group")
        if any(
            v is not None for v in (com, moi, rot_sym, rot_consts, point_group)
        ):
            lines.append("")
            lines.append(separator("Physical Properties"))
            if point_group is not None:
                lines.append(format_kv("Point Group", point_group))
            if rot_sym is not None:
                lines.append(format_kv("Rotational Symmetry Number", rot_sym))
            if com is not None:
                lines.append(
                    format_kv(
                        "Center of Mass (Å)",
                        f"[{', '.join(format_float(c, 6) for c in com)}]",
                    )
                )
            if moi is not None:
                lines.append(
                    format_kv(
                        "Moments of Inertia (amu·Å²)",
                        f"[{', '.join(format_float(v, 6) for v in moi)}]",
                    )
                )
            if rot_consts is not None:
                lines.append(
                    format_kv(
                        "Rotational Constants (Hz)",
                        f"[{', '.join(f'{v:.4e}' for v in rot_consts)}]",
                    )
                )

        # Dipole moment
        dipole = struct.get("dipole_moment")
        dipole_mag = struct.get("dipole_moment_magnitude")
        if dipole is not None:
            lines.append("")
            lines.append(separator("Dipole Moment (Debye)"))
            lines.append(f"  {'X':>10}  {'Y':>10}  {'Z':>10}  {'Tot':>10}")
            lines.append(
                f"  {'----------':>10}  {'----------':>10}  {'----------':>10}  {'----------':>10}"
            )
            x, y, z = dipole
            tot = dipole_mag if dipole_mag is not None else float("nan")
            lines.append(f"  {x:>10.6f}  {y:>10.6f}  {z:>10.6f}  {tot:>10.6f}")

        # Mulliken population analysis
        mulliken = struct.get("mulliken_atomic_charges")
        spin_densities = struct.get("mulliken_spin_densities")
        if mulliken is not None or spin_densities is not None:
            lines.append("")
            lines.append(separator("Mulliken Population Analysis"))
            if mulliken is not None and spin_densities is not None:
                lines.append(
                    f"  {'Idx':>4}  {'Elem':<4}  {'Charge':>10}  {'Spin':>10}"
                )
                lines.append(
                    f"  {'----':>4}  {'----':<4}  {'----------':>10}  {'----------':>10}"
                )
                for i, (atom, q) in enumerate(mulliken.items(), 1):
                    s = spin_densities.get(atom, float("nan"))
                    lines.append(
                        f"  {i:>4}  {atom:<4}  {q:>10.6f}  {s:>10.6f}"
                    )
            elif mulliken is not None:
                lines.append(f"  {'Idx':>4}  {'Elem':<4}  {'Charge':>10}")
                lines.append(f"  {'----':>4}  {'----':<4}  {'----------':>10}")
                for i, (atom, q) in enumerate(mulliken.items(), 1):
                    lines.append(f"  {i:>4}  {atom:<4}  {q:>10.6f}")
            else:
                lines.append(f"  {'Idx':>4}  {'Elem':<4}  {'Spin':>10}")
                lines.append(f"  {'----':>4}  {'----':<4}  {'----------':>10}")
                for i, (atom, s) in enumerate(spin_densities.items(), 1):
                    lines.append(f"  {i:>4}  {atom:<4}  {s:>10.6f}")

        # Vibrational data
        num_vib = struct.get("num_vibrational_modes")
        if num_vib is not None:
            freqs = struct.get("vibrational_frequencies", [])
            lines.append("")
            lines.append(
                separator(f"Vibrational Frequencies ({num_vib} modes)")
            )
            lines.append(f"  {'Mode':>6}  {'Frequency (cm⁻¹)':>18}")
            lines.append(f"  {'------':>6}  {'------------------':>18}")
            if freqs:
                for i, freq in enumerate(freqs):
                    lines.append(f"  {i + 1:>6}  {freq:>18.4f}")

        lines.append("")
        lines.append(separator())
        lines.append("")
        return "\n".join(lines)

    def molecule_detail(self):
        """Return molecule data with related structures and records.

        Returns:
            tuple: (molecule_dict, structures_list, records_list)
        """
        full_id = self.db.get_molecule_by_partial_id(self.molecule_id)
        molecule = self.db.get_molecule(full_id)
        if molecule is None:
            raise ValueError(f"Molecule not found: {full_id}")
        structures = self.db.get_structures_for_molecule(full_id)
        records = self.db.get_records_for_molecule(full_id)
        return molecule, structures, records

    def format_molecule_detail(self):
        """Return a human-readable molecule detail string."""
        molecule, structures, records = self.molecule_detail()
        lines = []
        lines.append("")
        lines.append(separator("Molecule Detail"))

        # Core attributes
        lines.append(format_kv("Molecule ID", molecule.get("molecule_id")))
        lines.append(
            format_kv("Chemical Formula", molecule.get("chemical_formula"))
        )
        lines.append(
            format_kv("Mass (amu)", format_float(molecule.get("mass"), 4))
        )

        # Compact boolean properties as tags
        tag_map = {
            "aromatic": (molecule.get("is_aromatic"), "nonaromatic"),
            "chiral": (molecule.get("is_chiral"), "achiral"),
            "linear": (molecule.get("is_linear"), "nonlinear"),
            "multicomponent": (
                molecule.get("is_multicomponent"),
                "monocomponent",
            ),
        }
        active_tags = sorted(k for k, (v, _) in tag_map.items() if v)
        inactive_tags = sorted(neg for k, (v, neg) in tag_map.items() if not v)
        tags_str = ", ".join(active_tags + inactive_tags) or "(none)"
        lines.append(format_kv("Tags", tags_str))

        # Representation
        lines.append("")
        lines.append(separator("Representation"))
        lines.append(format_kv("SMILES", molecule.get("smiles")))
        lines.append(format_kv("InChI", molecule.get("inchi")))

        # Compact element counts
        element_counts = molecule.get("element_counts")
        if element_counts:
            n_atoms = molecule.get("number_of_atoms", "")
            n_elems = len(element_counts)
            lines.append("")
            lines.append(
                separator(f"Composition: {n_elems} Elements, {n_atoms} Atoms")
            )
            counts_str = "  " + "  ".join(
                f"{elem} {cnt}" for elem, cnt in sorted(element_counts.items())
            )
            lines.append(counts_str)

        # Structures (conformers) table
        lines.append("")
        lines.append(separator(f"Structures: {len(structures)}"))
        if structures:
            sorted_structures = sort_structure_dicts_by_energy(
                self.db_file, structures
            )
            # Determine the primary (method, basis) label for the column header
            primary_mb = None
            if sorted_structures:
                primary_mb = sorted_structures[0].get("primary_method_basis")
            if primary_mb is not None:
                mb_label = "/".join(x for x in primary_mb if x) or "energy"
                energy_header = f"Energy[{mb_label}] (Eh)"
            else:
                energy_header = "Energy (Eh)"

            # Check if all conformers share the same electronic state
            charges = {s.get("charge") for s in sorted_structures}
            multiplicities = {s.get("multiplicity") for s in sorted_structures}
            mixed_states = len(charges) > 1 or len(multiplicities) > 1

            # Reference energy for ΔE: only meaningful for uniform electronic states
            if not mixed_states:
                primary_energies = [
                    s.get("primary_energy")
                    for s in sorted_structures
                    if s.get("primary_energy") is not None
                ]
                e_ref = min(primary_energies) if primary_energies else None
            else:
                e_ref = None  # cannot compare energies across different charge/spin states

            delta_header = "ΔE (—)" if mixed_states else "ΔE (kcal/mol)"
            lines.append(
                f"  {'Structure ID':<12}  {'Charge':>6}  {'Mult':>4}"
                f"  {energy_header:>32}  {delta_header:>14}"
            )
            lines.append(
                f"  {'------------':<12}  {'------':>6}  {'----':>4}"
                f"  {'--------------------------------':>32}  {'--------------':>14}"
            )
            for s in sorted_structures:
                sid = str(s.get("structure_id", ""))[:12]
                charge = s.get("charge", "")
                mult = s.get("multiplicity", "")
                primary_e = s.get("primary_energy")
                energy_str = (
                    format_energy(primary_e)
                    if primary_e is not None
                    else "N/A"
                )
                if primary_e is not None and e_ref is not None:
                    delta_str = f"{energy_conversion('hartree', 'kcal/mol', primary_e - e_ref):.2f}"
                elif mixed_states:
                    delta_str = (
                        "—"  # mixed charge/multiplicity -> ΔE is meaningless
                    )
                else:
                    delta_str = "N/A"  # uniform state but missing energy
                lines.append(
                    f"  {sid:<12}  {charge:>6}  {mult:>4}"
                    f"  {energy_str:>32}  {delta_str:>14}"
                )

            if mixed_states:
                lines.append("")
                lines.append(
                    "  Note: Mixed electronic states detected "
                    "(different charge or multiplicity). ΔE is not applicable."
                )
        else:
            lines.append("  (none)")

        db_name = os.path.basename(self.db_file)
        lines.append("")
        lines.append(
            f"  Tip: chemsmart run database inspect -f {db_name}"
            f" --sid <structure id>"
        )

        # Related records table
        lines.append("")
        lines.extend(self._format_related_records_table(records))

        lines.append("")
        lines.append(
            f"  Tip: chemsmart run database inspect -f {db_name}"
            f" --rid <record id>"
        )

        lines.append("")
        lines.append(separator())
        lines.append("")
        return "\n".join(lines)

    def standalone_structure_detail(self):
        """Return structure data with related records (standalone view).

        Unlike structure_detail() which shows a structure in record context
        (including per-calculation data like vibrational frequencies), this
        method presents the structure from the structures table directly.

        Returns:
            tuple: (structure_dict, records_list)
        """
        full_id = self.db.get_structure_by_partial_id(self.structure_id)
        structure = self.db.get_structure(full_id)
        if structure is None:
            raise ValueError(f"Structure not found: {full_id}")
        records = self.db.get_records_for_structure(full_id)
        return structure, records

    def format_standalone_structure_detail(self):
        """Return a human-readable standalone structure detail string."""
        structure, records = self.standalone_structure_detail()
        lines = []
        lines.append("")
        lines.append(separator("Structure Detail"))
        lines.append(format_kv("Structure ID", structure.get("structure_id")))
        lines.append(format_kv("Charge", structure.get("charge")))
        lines.append(format_kv("Multiplicity", structure.get("multiplicity")))

        # Molecule-level identity fields
        lines.append("")
        lines.append(separator("Parent Molecule"))
        lines.append(format_kv("Molecule ID", structure.get("molecule_id")))
        lines.append(
            format_kv("Chemical Formula", structure.get("chemical_formula"))
        )
        lines.append(
            format_kv("Mass (amu)", format_float(structure.get("mass"), 4))
        )
        lines.append(format_kv("SMILES", structure.get("smiles")))

        # Coordinates
        symbols = structure.get("chemical_symbols", [])
        positions = structure.get("positions", [])
        if symbols and positions:
            lines.append("")
            lines.append(separator(f"Coordinates ({len(symbols)} atoms, Å)"))
            lines.append(
                f"  {'Idx':>4}  {'Elem':<4}  "
                f"{'X':>14}  {'Y':>14}  {'Z':>14}"
            )
            lines.append(
                f"  {'----':>4}  {'----':<4}  "
                f"{'----------':>14}  {'----------':>14}  {'----------':>14}"
            )
            for i, (sym, pos) in enumerate(zip(symbols, positions)):
                x, y, z = pos
                lines.append(
                    f"  {i + 1:>4}  {sym:<4}  "
                    f"{x:>14.8f}  {y:>14.8f}  {z:>14.8f}"
                )

        # Physical properties
        com = structure.get("center_of_mass")
        moi = structure.get("moments_of_inertia")
        if com is not None or moi is not None:
            lines.append("")
            lines.append(separator("Physical Properties"))
            if com is not None:
                lines.append(
                    format_kv(
                        "Center of Mass (Å)",
                        f"[{', '.join(format_float(c, 6) for c in com)}]",
                    )
                )
            if moi is not None:
                lines.append(
                    format_kv(
                        "Moments of Inertia (amu·Å²)",
                        f"[{', '.join(format_float(v, 6) for v in moi)}]",
                    )
                )

        db_name = os.path.basename(self.db_file)
        # Related records table
        lines.append("")
        lines.extend(self._format_related_records_table(records))

        lines.append("")
        lines.append(
            f"  Tip: chemsmart run database inspect -f {db_name}"
            f" --rid <record id>"
        )

        lines.append("")
        lines.append(separator())
        lines.append("")
        return "\n".join(lines)

    def _format_related_records_table(self, records):
        """Format a Related Records table section.

        Args:
            records: list of record summary dicts from
                db.get_records_for_molecule() or db.get_records_for_structure()

        Returns:
            list[str]: formatted lines
        """
        lines = []
        lines.append(separator(f"Related Records: {len(records)}"))

        if not records:
            lines.append("  (none)")
            return lines

        lines.append(
            f"  {'Record ID':<12}  {'Job':<8}  {'Program':<8}  "
            f"{'Method':<12}  {'Basis':<14}  {'Energy (Eh)':<15}"
        )
        lines.append(
            f"  {'------------':<12}  {'--------':<8}  {'--------':<8}  "
            f"{'------------':<12}  {'--------------':<14}  {'---------------':<15}"
        )
        for r in records:
            rid = str(r.get("record_id", ""))[:12]
            job = r.get("jobtype", "") or ""
            prog = r.get("program", "") or ""
            func = r.get("method", "") or ""
            basis = r.get("basis", "") or ""
            energy = r.get("total_energy")
            energy_str = format_energy(energy) if energy is not None else ""
            lines.append(
                f"  {rid:<12}  {job:<8}  {prog:<8}  "
                f"{func:<12}  {basis:<14}  {energy_str:<15}"
            )
        return lines
