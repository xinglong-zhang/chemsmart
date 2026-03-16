"""
Database inspection module for viewing database overview, record details,
and molecule details.
"""

import logging
import os
import sqlite3

from chemsmart.assembler.database import Database
from chemsmart.assembler.utils import (
    bool_to_str,
    format_energy,
    format_float,
    format_kv,
    human_size,
    separator,
    truncate_iso,
)

logger = logging.getLogger(__name__)


class DatabaseInspector:
    """Inspect a chemsmart database.

    Provides three views:
    Overview – database-level metadata and statistics
    Record detail – full detail for one record
    Molecule detail – full detail for one molecule within a record
    """

    def __init__(self, db_file, index=None, record_id=None, molecule=None):
        self.db_file = db_file
        self.db = Database(db_file)
        self.index = index
        self.record_id = record_id
        self.molecule = molecule

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
        conn = sqlite3.connect(self.db_file)
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

            # Unique functionals and basis sets
            rows = conn.execute(
                "SELECT functional, COUNT(*) AS cnt FROM records WHERE functional IS NOT NULL GROUP BY functional ORDER BY cnt DESC"
            ).fetchall()
            stats["functionals"] = [(r["functional"], r["cnt"]) for r in rows]
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

            # Frequency and solvent count (from molecules table)
            stats["freq_count"] = conn.execute(
                "SELECT COUNT(DISTINCT record_id) FROM molecules WHERE vibrational_frequencies_json IS NOT NULL"
            ).fetchone()[0]
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
        """Return full record data with molecule summaries."""
        record_id = self.resolve_id()
        record = self.db.get_record(record_id=record_id)
        if record is None:
            raise ValueError(f"Record not found: {record_id}")
        return record

    def molecule_detail(self):
        """Return full molecule data for one molecule in a record."""
        record = self.record_detail()
        molecules = record.get("molecules", [])
        for mol in molecules:
            if mol.get("index") == self.molecule:
                return record, mol
        raise ValueError(
            f"Molecule index {self.molecule} not found in record "
            f"(available: {[m.get('index') for m in molecules]})."
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
        assembled_first = truncate_iso(stats.get("assembled_first"))
        assembled_last = truncate_iso(stats.get("assembled_last"))
        lines.append(format_kv("Created", f"{assembled_first} (UTC)"))
        lines.append(format_kv("Last Updated", f"{assembled_last} (UTC)"))

        # Programs
        lines.append("")
        lines.append(separator("Programs"))
        for prog, cnt in stats["programs"]:
            lines.append(f"  {prog or '(unknown)':<28}: {cnt}")

        # Functionals & Basis sets
        lines.append("")
        lines.append(separator("Electronic Structure Methods"))
        if stats["functionals"]:
            functionals_str = ", ".join(
                f"{f[0]} ({f[1]})" for f in stats["functionals"]
            )
            lines.append(format_kv("Methods", functionals_str))
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
            lines.append(f"  {jt or '(unknown)':<28}: {cnt}")

        # Calculation data
        lines.append("")
        lines.append(separator("Calculation Data"))
        lines.append(
            format_kv("Records with Energies", f"{stats['energy_count']}")
        )
        lines.append(
            format_kv("Records with Frequencies", f"{stats['freq_count']}")
        )
        lines.append(
            format_kv(
                "Records with Thermochemistry", f"{stats['thermo_count']}"
            )
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
            format_kv(
                "Filename", os.path.basename(provenance.get("source_file"))
            )
        )

        # Meta
        lines.append("")
        lines.append(separator("Calculation"))
        lines.append(format_kv("Program", provenance.get("program")))
        lines.append(
            format_kv("Program Version", provenance.get("program_version"))
        )
        lines.append(format_kv("Method", meta.get("functional")))
        lines.append(format_kv("Basis Set", meta.get("basis")))
        lines.append(format_kv("Spin", meta.get("spin")))
        lines.append(format_kv("Job Type", meta.get("jobtype")))
        lines.append(format_kv("Solvent", bool_to_str(meta.get("solvent_on"))))
        if meta.get("solvent_on"):
            lines.append(format_kv("Solvent Model", meta.get("solvent_model")))
            lines.append(format_kv("Solvent ID", meta.get("solvent_id")))
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
        lines.append(format_kv("Source File", provenance.get("source_file")))
        lines.append(
            format_kv("Source File Hash", provenance.get("source_file_hash"))
        )
        lines.append(
            format_kv(
                "Source File Size",
                human_size(provenance.get("source_file_size")),
            )
        )
        lines.append(
            format_kv("Source File Date", provenance.get("source_file_date"))
        )
        lines.append(format_kv("Assembled At", provenance.get("assembled_at")))

        # Molecule summary table
        if molecules:
            lines.append("")
            lines.append(separator(f"Molecules ({len(molecules)})"))

            # Check whether identity fields are uniform across all molecules
            formulas = {
                (mol.get("chemical_formula") or "") for mol in molecules
            }
            charges = {mol.get("charge") for mol in molecules}
            mults = {mol.get("multiplicity") for mol in molecules}
            natoms = {mol.get("number_of_atoms") for mol in molecules}
            smiles_set = {(mol.get("smiles") or "") for mol in molecules}
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
                natom = next(iter(natoms))
                smiles = next(iter(smiles_set))
                lines.append(format_kv("Formula", formula or "N/A"))
                lines.append(format_kv("Charge", charge))
                lines.append(format_kv("Multiplicity", mult))
                lines.append(format_kv("Atoms", natom))
                if smiles:
                    lines.append(format_kv("SMILES", smiles))
                lines.append("")
                lines.append(
                    f"  {'Idx':>4}  {'Energy (Eh)':>20}  {'Optimized':>9}"
                )
                lines.append(
                    f"  {'----':>4}  {'--------------------':>20}  {'---------':>9}"
                )
                for mol in molecules:
                    idx = mol.get("index", "-")
                    energy = mol.get("energy")
                    energy_str = (
                        format_energy(energy) if energy is not None else "N/A"
                    )
                    is_opt = mol.get("is_optimized_structure")
                    opt_str = bool_to_str(is_opt) if is_opt is not None else ""
                    lines.append(f"  {idx:>4}  {energy_str:>20}  {opt_str:>9}")
            else:
                # Heterogeneous molecules — show full info
                hdr = (
                    f"  {'Idx':>4}  {'Formula':<20}  {'Charge':>6}  "
                    f"{'Mult':>4}  {'Atoms':>5}  {'Energy (Eh)':>20}  {'SMILES':<30}"
                )
                lines.append(hdr)
                lines.append(
                    f"  {'----':>4}  {'-------':<20}  {'------':>6}  "
                    f"{'----':>4}  {'-----':>5}  {'--------------------':>20}  {'------':<30}"
                )
                for mol in molecules:
                    idx = mol.get("index", "-")
                    formula = mol.get("chemical_formula", "") or ""
                    charge = mol.get("charge", "")
                    mult = mol.get("multiplicity", "")
                    natom = mol.get("number_of_atoms", "")
                    energy = mol.get("energy")
                    energy_str = (
                        format_energy(energy) if energy is not None else ""
                    )
                    smiles = mol.get("smiles", "") or ""
                    lines.append(
                        f"  {idx:>4}  {formula:<20}  {charge:>6}  "
                        f"{mult:>4}  {natom:>5}  {energy_str:>20}  {smiles:<30}"
                    )

        lines.append("")
        lines.append(separator())
        lines.append("")
        return "\n".join(lines)

    def format_molecule_detail(self):
        """Return a human-readable molecule detail string."""
        record, mol = self.molecule_detail()
        lines = []
        lines.append("")
        lines.append(separator("Molecule Detail"))

        # Context
        lines.append(format_kv("Record Index", record.get("record_index")))
        lines.append(format_kv("Record ID", record.get("record_id")))
        lines.append(
            format_kv(
                "Structure Index in File", mol.get("structure_index_in_file")
            )
        )

        # Basic info
        lines.append("")
        lines.append(separator("Identity"))
        lines.append(format_kv("Molecule Index", mol.get("index")))
        lines.append(
            format_kv("Chemical Formula", mol.get("chemical_formula"))
        )
        lines.append(format_kv("SMILES", mol.get("smiles")))
        lines.append(format_kv("Charge", mol.get("charge")))
        lines.append(format_kv("Multiplicity", mol.get("multiplicity")))
        lines.append(format_kv("Number of Atoms", mol.get("number_of_atoms")))
        lines.append(format_kv("Mass (amu)", format_float(mol.get("mass"), 4)))
        lines.append(
            format_kv("Energy (Eh)", format_energy(mol.get("energy")))
        )

        # Coordinates
        symbols = mol.get("chemical_symbols", [])
        positions = mol.get("positions", [])
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
        frozen = mol.get("frozen_atoms")
        if frozen:
            lines.append("")
            lines.append(separator("Frozen Atoms"))
            lines.append(f"  {frozen}")

        # Structure properties
        lines.append("")
        lines.append(separator("Structure Properties"))
        lines.append(
            format_kv(
                "Is Optimized Structure",
                bool_to_str(mol.get("is_optimized_structure")),
            )
        )
        lines.append(format_kv("Is Chiral", bool_to_str(mol.get("is_chiral"))))
        lines.append(format_kv("Is Ring", bool_to_str(mol.get("is_ring"))))
        lines.append(
            format_kv("Is Aromatic", bool_to_str(mol.get("is_aromatic")))
        )
        lines.append(
            format_kv("Is Monoatomic", bool_to_str(mol.get("is_monoatomic")))
        )
        lines.append(
            format_kv("Is Diatomic", bool_to_str(mol.get("is_diatomic")))
        )
        lines.append(format_kv("Is Linear", bool_to_str(mol.get("is_linear"))))

        # Element counts
        element_counts = mol.get("element_counts")
        if element_counts:
            lines.append("")
            lines.append(separator("Element Counts"))
            for elem, cnt in sorted(element_counts.items()):
                lines.append(f"  {elem:<4}: {cnt}")

        # Center of mass / moments of inertia
        com = mol.get("center_of_mass")
        if com is not None:
            lines.append("")
            lines.append(separator("Physical Properties"))
            lines.append(
                format_kv(
                    "Center of Mass (Å)",
                    f"[{', '.join(format_float(c, 6) for c in com)}]",
                )
            )
        moi = mol.get("moments_of_inertia")
        if moi is not None:
            lines.append(
                format_kv(
                    "Moments of Inertia (amu·Å²)",
                    f"[{', '.join(format_float(v, 6) for v in moi)}]",
                )
            )
        rot_sym = mol.get("rotational_symmetry_number")
        if rot_sym is not None:
            lines.append(format_kv("Rotational Symmetry Number", rot_sym))

        # Mulliken charges
        mulliken = mol.get("mulliken_atomic_charges")
        if mulliken is not None:
            lines.append("")
            lines.append(separator("Mulliken Atomic Charges"))
            lines.append(f"  {'Idx':>4}  {'Elem':<4}  {'Charge':>12}")
            lines.append(f"  {'----':>4}  {'----':<4}  {'----------':>12}")
            for i, (atom, q) in enumerate(mulliken.items(), 1):
                lines.append(f"  {i:>4}  {atom:<6}  {q:>10.6f}")

        # Vibrational data
        num_vib = mol.get("num_vibrational_modes")
        if num_vib is not None:
            freqs = mol.get("vibrational_frequencies", [])
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
