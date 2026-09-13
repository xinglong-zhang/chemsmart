#!/usr/bin/env python
"""PySCF's own numbers for an archived ChemSmart PySCF artifact.

Runs in the compute env, imports pyscf, never chemsmart. Rebuilds the
Mole and the mean field from the artifact's spec at the artifact's own
(final) geometry, then recomputes what the artifact claims with PySCF's
public API: the harmonic analysis and thermochemistry from the stored
Hessian, the TDA/TDDFT roots, oscillator strengths and transition dipoles
of a response stage, the MP2 / CCSD / CCSD(T) components of a correlated
stage, and the point group -- an independent representation of the same
bytes for the differential oracles in tests/. Writes reference.json
beside the .h5 with the command that produced it.

Usage: /opt/miniforge3/envs/chemsmart-pyscf/bin/python reference.py FILE.h5
"""

import json
import sys

import h5py
import numpy as np
import pyscf
from pyscf import cc, dft, gto, mp, scf, symm, tdscf
from pyscf.data import nist
from pyscf.hessian import thermo


def _read(handle, name):
    if name not in handle:
        return None
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


def _mean_field(spec, mol):
    """The reference the driver built, from the applied spec alone."""
    xc = spec.get("xc")
    if xc is None:
        mf = scf.HF(mol)
    else:
        mf = dft.KS(mol, xc=xc)
        grid = spec.get("atom_grid")
        if grid:
            mf.grids.atom_grid = tuple(int(item) for item in grid)
    if spec.get("scf_tol") is not None:
        mf.conv_tol = float(spec["scf_tol"])
    if spec.get("scf_maxiter") is not None:
        mf.max_cycle = int(spec["scf_maxiter"])
    call = spec.get("solvent_call")
    if call == "PCM":
        mf = mf.PCM()
        mf.with_solvent.method = spec["solvent_method"]
        mf.with_solvent.eps = float(spec["solvent_eps"])
    elif call == "SMD":
        mf = mf.SMD()
        mf.with_solvent.solvent = spec["solvent_id"]
    return mf


def _response(spec, mf):
    method = str(spec["response_method"]).strip().lower()
    manifold = str(spec["state_manifold"]).strip().lower()
    td = {"tda": tdscf.TDA, "tddft": tdscf.TDDFT}[method](mf)
    if manifold in ("singlet", "triplet"):
        td.singlet = manifold == "singlet"
    td.nstates = int(spec["nstates"])
    if spec.get("td_max_cycle") is not None:
        td.max_cycle = int(spec["td_max_cycle"])
    td.kernel()
    return td


def _correlated(spec, mf):
    method = str(spec["ab_initio"]).strip().lower()
    obj = mp.MP2(mf) if method == "mp2" else cc.CCSD(mf)
    if method != "mp2" and spec.get("cc_max_cycle") is not None:
        obj.max_cycle = int(spec["cc_max_cycle"])
    frozen = spec.get("frozen_core")
    if isinstance(frozen, str) and frozen.strip().lower() == "auto":
        obj.set_frozen(method="auto")
    elif frozen is not None:
        obj.frozen = int(frozen)
    obj.kernel()
    return method, obj


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
                "xc",
                "ab_initio",
                "atom_grid",
                "scf_tol",
                "scf_maxiter",
                "solvent_call",
                "solvent_method",
                "solvent_id",
                "solvent_eps",
                "response_method",
                "state_manifold",
                "nstates",
                "excited_state_root",
                "td_max_cycle",
                "frozen_core",
                "cc_max_cycle",
            )
        }
        results = {
            key: (
                np.asarray(handle[f"results/{key}"][()])
                if key in handle["results"]
                else None
            )
            for key in (
                "hessian",
                "vibrational_frequencies",
                "energies",
                "excitation_energies",
                "correlation_energy",
                "total_energy",
            )
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
    needs_reference = (
        results["excitation_energies"] is not None
        or results["correlation_energy"] is not None
    )
    if needs_reference:
        mf = _mean_field(spec, mol)
        mf.kernel()
        out["scf_energy_eh"] = float(mf.e_tot)
        out["scf_converged"] = bool(mf.converged)
    if results["excitation_energies"] is not None:
        td = _response(spec, mf)
        excitations = np.asarray(td.e, dtype=float).reshape(-1)
        out["td_excitation_energies_eh"] = excitations.tolist()
        out["td_converged_per_root"] = [
            bool(item)
            for item in np.asarray(td.converged, dtype=bool).reshape(-1)
        ]
        out["td_positive_eig_threshold_eh"] = float(td.positive_eig_threshold)
        out["td_nstates_requested"] = int(spec["nstates"])
        out["td_singlet_flag"] = (
            None if td.singlet is None else bool(td.singlet)
        )
        try:
            out["td_oscillator_strengths"] = (
                np.asarray(td.oscillator_strength(), dtype=float)
                .reshape(-1)[: excitations.size]
                .tolist()
            )
        except Exception as exc:  # noqa: BLE001 - recorded, never hidden
            out["td_oscillator_strengths_error"] = type(exc).__name__
        try:
            dipoles = np.asarray(td.transition_dipole(), dtype=float)
            out["td_transition_dipole_moments_debye"] = (
                dipoles.reshape(excitations.size, -1)[:, :3]
                * float(nist.AU2DEBYE)
            ).tolist()
        except Exception as exc:  # noqa: BLE001
            out["td_transition_dipole_moments_error"] = type(exc).__name__
        response_solvent = getattr(td, "with_solvent", None)
        if response_solvent is not None:
            out["td_response_eps"] = float(response_solvent.eps)
            out["td_equilibrium_solvation"] = bool(
                getattr(response_solvent, "equilibrium_solvation", False)
            )
        stored = np.asarray(results["excitation_energies"], dtype=float)
        common = min(stored.size, excitations.size)
        out["stored_excitations_match_recomputed_max_abs_eh"] = (
            float(np.max(np.abs(stored[:common] - excitations[:common])))
            if common
            else None
        )
        root = spec.get("excited_state_root")
        if root is not None and excitations.size >= int(root):
            out["followed_root_total_energy_eh"] = float(
                td.e_tot[int(root) - 1]
            )
    if results["correlation_energy"] is not None:
        method, obj = _correlated(spec, mf)
        out["corr_method"] = method
        out["corr_reference_energy_eh"] = float(mf.e_tot)
        correlation = float(obj.e_corr)
        out["corr_converged"] = bool(getattr(obj, "converged", True))
        frozen = getattr(obj, "frozen", None)
        out["corr_frozen_orbitals"] = (
            0
            if frozen is None
            else (
                int(frozen)
                if isinstance(frozen, (int, np.integer))
                else len(frozen)
            )
        )
        if method in ("ccsd", "ccsd(t)"):
            out["corr_ccsd_correlation_energy_eh"] = correlation
        if method == "ccsd(t)":
            triples = float(obj.ccsd_t())
            out["corr_triples_correction_eh"] = triples
            correlation += triples
        out["corr_correlation_energy_eh"] = correlation
        out["corr_total_energy_eh"] = float(mf.e_tot) + correlation
        out["stored_total_energy_match_recomputed_abs_eh"] = float(
            abs(float(results["total_energy"]) - out["corr_total_energy_eh"])
        )
    target = path.rsplit(".", 1)[0] + ".reference.json"
    json.dump(out, open(target, "w"), indent=1, sort_keys=True)
    print("reference written for", name)


if __name__ == "__main__":
    main(sys.argv[1])
