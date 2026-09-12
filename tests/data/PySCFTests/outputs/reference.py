#!/usr/bin/env python
"""PySCF's own numbers for an archived ChemSmart PySCF artifact.

Runs in the compute env, imports pyscf, never chemsmart. Rebuilds the
Mole from the artifact's spec, recomputes the harmonic analysis and
thermochemistry from the stored Hessian with pyscf.hessian.thermo, and
re-detects the point group -- an independent representation of the same
bytes for the differential oracles in tests/. Writes reference.json
beside the .h5 with the command that produced it."""
import json, sys
import h5py, numpy as np, pyscf
from pyscf import gto, symm
from pyscf.hessian import thermo

path = sys.argv[1]
with h5py.File(path, "r") as h:
    def read(name):
        node = h[name]
        if bool(node.attrs.get("chemsmart_is_null", False)): return None
        v = node[()]
        return v.decode() if isinstance(v, bytes) else (v.tolist() if hasattr(v, "tolist") else v)
    spec = {k: read(f"spec/{k}") for k in ("symbols","positions","unit","basis","charge","spin","multiplicity","jobtype","stages","xc")}
    results = {k: (np.asarray(h[f"results/{k}"][()]) if k in h["results"] else None) for k in ("hessian","vibrational_frequencies","energies","positions")}
symbols = [s.decode() if isinstance(s, bytes) else s for s in spec["symbols"]]
mol = gto.M(atom=[(s, tuple(p)) for s, p in zip(symbols, results["positions"].tolist())], basis=spec["basis"], charge=int(spec["charge"]), spin=int(spec["spin"]), unit="Angstrom", verbose=0)
out = {"pyscf_version": pyscf.__version__, "artifact": path.split("/")[-1], "jobtype": spec["jobtype"], "point_group_detect_symm": symm.detect_symm(mol._atom)[0], "symm_geom_tol_bohr": float(symm.geom.TOLERANCE), "masses_isotope_averaged_amu": mol.atom_mass_list(isotope_avg=True).tolist(), "command": "reference.py " + path.split("/")[-1]}
if results["hessian"] is not None:
    analysis = thermo.harmonic_analysis(mol, results["hessian"], imaginary_freq=False)
    out["harmonic_analysis_freq_wavenumber"] = np.asarray(analysis["freq_wavenumber"], dtype=float).tolist()
    out["reduced_mass_amu"] = np.asarray(analysis["reduced_mass"], dtype=float).tolist()
    class _Model:  # thermo.thermo reads .mol and .e_tot only
        pass
    model = _Model(); model.mol = mol; model.e_tot = float(results["energies"][-1])
    th = thermo.thermo(model, analysis["freq_au"], 298.15, 101325)
    def _plain(value):
        if isinstance(value, (list, tuple)):
            return [_plain(item) for item in value]
        if hasattr(value, "tolist"):
            return value.tolist()
        return value
    out["thermo_298K_1atm"] = {k: _plain(v) for k, v in th.items()}
    out["stored_frequencies_match_recomputed_max_abs_cm1"] = float(np.max(np.abs(np.sort(np.asarray(results["vibrational_frequencies"])) - np.sort(np.asarray(analysis["freq_wavenumber"], dtype=float)))))
json.dump(out, open(path.rsplit(".", 1)[0] + ".reference.json", "w"), indent=1, sort_keys=True)
print("reference written for", out["artifact"], "| group", out["point_group_detect_symm"], "| freqs", out.get("harmonic_analysis_freq_wavenumber"))
