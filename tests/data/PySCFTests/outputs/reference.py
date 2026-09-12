#!/usr/bin/env python
"""PySCF's own numbers for an archived ChemSmart PySCF artifact.

Runs in the compute env, imports pyscf, never chemsmart. Rebuilds the
Mole from the artifact's spec, recomputes the harmonic analysis and
thermochemistry from the stored Hessian with pyscf.hessian.thermo, and
re-detects the point group -- an independent representation of the same
bytes for the differential oracles in tests/. Writes reference.json
beside the .h5 with the command that produced it.

Usage: /opt/miniforge3/envs/chemsmart-pyscf/bin/python reference.py FILE.h5
"""

import json
import sys

import h5py
import numpy as np
import pyscf
from pyscf import gto, symm
from pyscf.hessian import thermo


def _read(handle, name):
    node = handle[name]
    if bool(node.attrs.get("chemsmart_is_null", False)):
        return None
    value = node[()]
    if isinstance(value, bytes):
        return value.decode()
    return value.tolist() if hasattr(value, "tolist") else value


def _plain(value):
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if hasattr(value, "tolist"):
        return value.tolist()
    return value


class _Model:
    """thermo.thermo reads .mol and .e_tot only."""


def main(path):
    with h5py.File(path, "r") as handle:
        spec = {
            key: _read(handle, f"spec/{key}")
            for key in (
                "symbols",
                "basis",
                "charge",
                "spin",
                "multiplicity",
                "jobtype",
            )
        }
        results = {
            key: (
                np.asarray(handle[f"results/{key}"][()])
                if key in handle["results"]
                else None
            )
            for key in ("hessian", "vibrational_frequencies", "energies")
        }
        positions = np.asarray(handle["results/positions"][()])
    symbols = [
        s.decode() if isinstance(s, bytes) else s for s in spec["symbols"]
    ]
    mol = gto.M(
        atom=[(s, tuple(p)) for s, p in zip(symbols, positions.tolist())],
        basis=spec["basis"],
        charge=int(spec["charge"]),
        spin=int(spec["spin"]),
        unit="Angstrom",
        verbose=0,
    )
    name = path.split("/")[-1]
    out = {
        "pyscf_version": pyscf.__version__,
        "artifact": name,
        "jobtype": spec["jobtype"],
        "point_group_detect_symm": symm.detect_symm(mol._atom)[0],
        "symm_geom_tol_bohr": float(symm.geom.TOLERANCE),
        "masses_isotope_averaged_amu": mol.atom_mass_list(
            isotope_avg=True
        ).tolist(),
        "command": "reference.py " + name,
    }
    if results["hessian"] is not None:
        analysis = thermo.harmonic_analysis(
            mol, results["hessian"], imaginary_freq=False
        )
        recomputed = np.asarray(analysis["freq_wavenumber"], dtype=float)
        out["harmonic_analysis_freq_wavenumber"] = recomputed.tolist()
        out["reduced_mass_amu"] = np.asarray(
            analysis["reduced_mass"], dtype=float
        ).tolist()
        model = _Model()
        model.mol = mol
        model.e_tot = float(results["energies"][-1])
        thermo_values = thermo.thermo(
            model, analysis["freq_au"], 298.15, 101325
        )
        out["thermo_298K_1atm"] = {
            k: _plain(v) for k, v in thermo_values.items()
        }
        stored = np.sort(np.asarray(results["vibrational_frequencies"]))
        out["stored_frequencies_match_recomputed_max_abs_cm1"] = float(
            np.max(np.abs(stored - np.sort(recomputed)))
        )
    target = path.rsplit(".", 1)[0] + ".reference.json"
    json.dump(out, open(target, "w"), indent=1, sort_keys=True)
    print("reference written for", name)


if __name__ == "__main__":
    main(sys.argv[1])
