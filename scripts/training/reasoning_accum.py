"""Accumulate verified reasoning-synthesis episodes from DashScope teachers.

This is a live-data collection helper for v17 reasoning distillation. It runs
the real ChemSmart unified agent loop against Alibaba/DashScope models, stores
raw episodes in a model-specific repo-local training directory, and appends a
small result ledger for quick triage.

The raw training ledger remains append-only. Use ``export_sft.py`` afterwards
to split trusted ``synthesis.reasoning`` from review-only reasoning candidates.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from dotenv import dotenv_values
from goal_batches import BATCHES

REPO = Path(__file__).resolve().parents[2]
BASE_URL = "https://dashscope-intl.aliyuncs.com/compatible-mode/v1"

# provider -> (api.env key name, base_url). Always pin an explicit base_url:
# OpenAIProvider falls back to a third-party gateway when base_url is None, and
# the official key must never leave api.openai.com.
PROVIDERS: dict[str, tuple[str, str]] = {
    "dashscope": ("alibaba-cloud-api-key", BASE_URL),
    "openai": ("OPEN-AI-api-key", "https://api.openai.com/v1"),
    "deepseek": ("DEEPSEEK-api-key", "https://api.deepseek.com/v1"),
}

# USD per 1M tokens, keyed by model-id prefix. Cached input is billed at the
# discounted rate; OpenAI reports it under prompt_tokens_details.cached_tokens.
PRICING: dict[str, tuple[float, float, float]] = {
    # model prefix: (input, cached_input, output)
    "gpt-5.4-mini": (0.75, 0.075, 4.50),
    "gpt-5.4-nano": (0.15, 0.015, 1.20),
}


class BudgetExceeded(RuntimeError):
    """Raised when accumulated spend crosses the configured cap."""


class CostMeter:
    """Meter token spend for a priced model and stop the run at the cap."""

    def __init__(self, model: str, budget_usd: float) -> None:
        self.model = model
        self.budget_usd = budget_usd
        self.calls = 0
        self.prompt_tokens = 0
        self.cached_tokens = 0
        self.completion_tokens = 0
        rate = None
        for prefix, value in PRICING.items():
            if model.startswith(prefix):
                rate = value
                break
        self.rate = rate

    @property
    def cost_usd(self) -> float:
        if self.rate is None:
            return 0.0
        rate_in, rate_cached, rate_out = self.rate
        fresh = max(self.prompt_tokens - self.cached_tokens, 0)
        return (
            fresh * rate_in
            + self.cached_tokens * rate_cached
            + self.completion_tokens * rate_out
        ) / 1_000_000

    def record(self, usage: Any) -> None:
        if not usage:
            return
        get = (
            usage.get
            if isinstance(usage, dict)
            else lambda k, d=None: getattr(usage, k, d)
        )
        self.calls += 1
        self.prompt_tokens += int(get("prompt_tokens", 0) or 0)
        self.completion_tokens += int(get("completion_tokens", 0) or 0)
        details = get("prompt_tokens_details", None) or {}
        cached = (
            details.get("cached_tokens")
            if isinstance(details, dict)
            else getattr(details, "cached_tokens", 0)
        )
        self.cached_tokens += int(cached or 0)

    def check(self) -> None:
        if self.rate is not None and self.cost_usd >= self.budget_usd:
            raise BudgetExceeded(
                f"spend ${self.cost_usd:.4f} reached cap ${self.budget_usd:.2f}"
            )

    def summary(self) -> dict[str, Any]:
        return {
            "model": self.model,
            "calls": self.calls,
            "prompt_tokens": self.prompt_tokens,
            "cached_tokens": self.cached_tokens,
            "completion_tokens": self.completion_tokens,
            "cost_usd": round(self.cost_usd, 4),
            "budget_usd": self.budget_usd,
            "priced": self.rate is not None,
        }


H2O_XYZ = (
    "3\n"
    "water-like placeholder\n"
    "O 0.000000 0.000000 0.117000\n"
    "H 0.000000 0.757000 -0.464000\n"
    "H 0.000000 -0.757000 -0.464000\n"
)


@dataclass(frozen=True)
class CorpusItem:
    kind: str
    program: str
    request: str
    workflow: str
    difficulty: str = "hard"


AGENT_DOMAIN_CORPUS: list[CorpusItem] = [
    CorpusItem(
        "sp",
        "orca",
        "Run a DLPNO-CCSD(T) single-point energy calculation on the optimized iron complex. I need a complete basis set (CBS) extrapolation using def2-TZVPP and def2-QZVPP, along with the corresponding def2-TZVPP/C and def2-QZVPP/C auxiliary basis sets.",
        "dlpno_ccsd_t_cbs",
    ),
    CorpusItem(
        "ts",
        "gaussian",
        "Locate the transition state for the concerted [4+2] cycloaddition using the Berny algorithm with a pre-calculated analytical Hessian (calcall).",
        "ts_opt_calcall",
    ),
    CorpusItem(
        "irc",
        "gaussian",
        "Wait, freeze the core backbone atoms using modredundant. After finding the saddle point, confirm it has exactly one imaginary frequency corresponding to the proton transfer, then run a MaxPoints=30 IRC.",
        "irc_validation",
    ),
    CorpusItem(
        "opt",
        "gaussian",
        "Optimize this fluorinated intermediate using the SMD solvation model for 1,2-dichloroethane.",
        "implicit_solvation_smd",
    ),
    CorpusItem(
        "nci",
        "gaussian",
        "Run a Hirshfeld population analysis alongside a standard Mulliken charge analysis on this transition metal complex, then perform NCI analysis.",
        "nci_population",
    ),
    CorpusItem(
        "sp",
        "gaussian",
        "Evaluate the spin-state energetics of this octahedral Fe(II) complex. I need geometry optimizations for the low-spin (S=0), intermediate-spin (S=1), and high-spin (S=2) states.",
        "spin_state_energetics",
    ),
    CorpusItem(
        "opt",
        "orca",
        "Set up an open-shell broken-symmetry DFT calculation for a dinuclear copper complex, charge +2, multiplicity 1 (antiferromagnetic coupling). Use B3LYP with def2-SVP.",
        "broken_symmetry_dft",
    ),
    CorpusItem(
        "modred",
        "orca",
        "Perform a relaxed surface scan on the C-C bond length from 1.5 to 2.5 Angstroms in 10 steps using ORCA, neutral singlet.",
        "relaxed_surface_scan",
    ),
    CorpusItem(
        "sp",
        "gaussian",
        "Calculate the NMR shielding tensors for this organic molecule using the GIAO method at the m062x/def2-tzvp level.",
        "nmr_giao",
    ),
    CorpusItem(
        "ts",
        "orca",
        "Optimize the transition state for the bulky chiral phosphoric acid catalyzed step. Use Cartesian coordinates with the GEDIIS algorithm.",
        "ts_opt_gediis",
    ),
]

WORKFLOW_HARD_CORPUS: list[CorpusItem] = [
    CorpusItem(
        kind="td",
        program="gaussian",
        workflow="excited_state_spectroscopy",
        request=(
            "Set up Gaussian TD-DFT for chromophore_a.xyz, singlets, 12 states, "
            "root 2, eqsolv, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="scan",
        program="gaussian",
        workflow="pes_scan_with_constraint",
        request=(
            "Prepare a relaxed Gaussian scan of bond 1-2 in alcohol_scan.xyz for "
            "10 steps of 0.05 A while keeping bond 1-3 frozen, neutral singlet."
        ),
    ),
    CorpusItem(
        kind="modred",
        program="orca",
        workflow="constrained_optimization",
        request=(
            "Use ORCA to optimize cage.xyz while freezing the distance between "
            "atoms 4 and 9, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="neb",
        program="orca",
        workflow="reaction_path_neb",
        request=(
            "Set up an ORCA NEB-TS path from reactant_neb.xyz to product_neb.xyz "
            "with 8 images and XTB2 pre-optimization, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="qmmm",
        program="orca",
        workflow="multiscale_qmmm",
        request=(
            "For enzyme_cluster.xyz, prepare an ORCA QM/MM optimization with "
            "atoms 1-8 as the high-level region, total charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="dias",
        program="gaussian",
        workflow="activation_strain_dias",
        request=(
            "Run Gaussian DIAS on cycloadd_ts.xyz in TS mode with fragment 1 "
            "atoms 1-5 and fragment 2 atoms 6-10, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="wbi",
        program="gaussian",
        workflow="bond_order_analysis",
        request=(
            "I need Wiberg bond indices for cationic_intermediate.xyz, Gaussian, "
            "charge +1 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="irc",
        program="gaussian",
        workflow="ts_connectivity_irc",
        request=(
            "Take barrier_guess.xyz and follow the Gaussian IRC in both directions "
            "to confirm connectivity, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="ts",
        program="orca",
        workflow="scan_guided_ts",
        request=(
            "ORCA ScanTS for ring_opening.xyz along bond 2-7 from 1.8 to 2.6 A "
            "with 9 points, neutral singlet."
        ),
    ),
    CorpusItem(
        kind="resp",
        program="gaussian",
        workflow="force_field_charge_derivation",
        request=(
            "Generate RESP charge input for ligand_resp.xyz using Gaussian, "
            "charge -1 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="ask",
        program="gaussian",
        workflow="clarification_missing_charge_spin",
        request="Run the calculation for molecule.xyz with the usual charge and spin.",
    ),
    CorpusItem(
        kind="opt",
        program="orca",
        workflow="multilingual_geometry_optimization",
        request="用 ORCA 优化 nickel_complex.xyz，电荷 0，多重度 3。",
    ),
]

COMMAND_HARD_CORPUS: list[CorpusItem] = [
    CorpusItem(
        kind="td",
        program="gaussian",
        workflow="excited_state_spectroscopy",
        request=(
            "Using the loaded demo project settings, make the ChemSmart CLI "
            "command for Gaussian TD-DFT on chromophore_cmd.xyz: singlets, "
            "12 states, root 2, eqsolv, charge 0 multiplicity 1. Do not build "
            "a new project YAML."
        ),
    ),
    CorpusItem(
        kind="scan",
        program="gaussian",
        workflow="pes_scan_with_constraint",
        request=(
            "Using the active Gaussian demo project, prepare the CLI command for "
            "a relaxed scan of bond 1-2 in alcohol_cmd.xyz for 10 steps of "
            "0.05 A while freezing bond 1-3, neutral singlet. Do not initialize "
            "or edit project YAML."
        ),
    ),
    CorpusItem(
        kind="modred",
        program="orca",
        workflow="constrained_optimization",
        request=(
            "With the active ORCA demo project already loaded, create the "
            "ChemSmart run command to optimize cage_cmd.xyz while freezing the "
            "distance between atoms 4 and 9, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="neb",
        program="orca",
        workflow="reaction_path_neb",
        request=(
            "Use the loaded ORCA demo project and synthesize only the CLI command "
            "for an NEB-TS path from reactant_cmd.xyz to product_cmd.xyz with "
            "8 images and XTB2 pre-optimization, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="qmmm",
        program="orca",
        workflow="multiscale_qmmm",
        request=(
            "Using current ORCA project settings, make the CLI command for a "
            "QM/MM optimization of enzyme_cmd.xyz with atoms 1-8 as high-level "
            "atoms, total charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="dias",
        program="gaussian",
        workflow="activation_strain_dias",
        request=(
            "Using active Gaussian demo project settings, synthesize the command "
            "for DIAS in TS mode on cycloadd_cmd.xyz with fragment 1 atoms 1-5 "
            "and fragment 2 atoms 6-10, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="wbi",
        program="gaussian",
        workflow="bond_order_analysis",
        request=(
            "Using the loaded Gaussian project, create the CLI command for Wiberg "
            "bond indices on cationic_cmd.xyz, charge +1 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="irc",
        program="gaussian",
        workflow="ts_connectivity_irc",
        request=(
            "Using the active Gaussian project, prepare the CLI command to follow "
            "the IRC from barrier_cmd.xyz in both directions, charge 0 "
            "multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="crest",
        program="gaussian",
        workflow="conformer_sampling_to_dft",
        request=(
            "Using the loaded Gaussian demo project, synthesize only the CLI "
            "command for a CREST conformer search on flexible_cmd.xyz followed "
            "by Gaussian opt evaluation, neutral singlet."
        ),
    ),
    CorpusItem(
        kind="nci",
        program="gaussian",
        workflow="noncovalent_interaction_analysis",
        request=(
            "With the active Gaussian demo project, create the command for NCI "
            "analysis of pi_stack_cmd.xyz, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="resp",
        program="gaussian",
        workflow="force_field_charge_derivation",
        request=(
            "Using the loaded Gaussian project, create the command for RESP "
            "charge fitting on ligand_cmd.xyz, charge -1 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="ts",
        program="gaussian",
        workflow="transition_state_frequency_confirmation",
        request=(
            "Using the current Gaussian project, create a TS optimization "
            "command for proton_transfer_cmd.xyz with maxstep 8, charge 0 "
            "multiplicity 1. Do not create YAML."
        ),
    ),
]

# Targets the coverage holes measured on 2026-07-10: qrc/traj had zero gated
# positives, and the kinds below sat under ten rows each. Every request uses
# command-hard framing (a loaded project, "command only") because plain job
# phrasing gets misrouted into project-YAML authoring.

AGENT_QMMM_CORPUS: list[CorpusItem] = [
    CorpusItem(
        "qmmm",
        "orca",
        "Run a QM/MM optimization of the enzyme-substrate complex in enzyme.xyz. Treat atoms 1-25 at the QM level (B3LYP/def2-SVP) and the rest at the MM level. I have the force field file prepared.",
        "orca_qmmm_basic",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Set up an ORCA QM/MM frequency calculation for protein_opt.xyz. The high-level region is atoms 15-40. Charge is -1 and multiplicity is 2, but wait, the QM region itself has a charge of -1 and multiplicity 2, while the MM region is neutral. Figure out how to specify this.",
        "orca_qmmm_charge_spin",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "gaussian",
        "I want a 2-layer ONIOM (B3LYP/6-31G*:UFF) optimization on complex.xyz. The real system has charge +1 and multiplicity 1, and the model system is neutral singlet. Make sure the charge/multiplicity scope is correctly mapped.",
        "gaussian_oniom_opt",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "gaussian",
        "Gaussian QM/MM (ONIOM) single point energy on the snapshot from my MD simulation (snapshot_100.xyz). High level: M06-2X/def2-TZVP. Low level: AMBER. The total system charge is -2, multiplicity 1. Model system charge is -1, multiplicity 1.",
        "gaussian_oniom_sp",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Do a QM/MM MD simulation step using ORCA for solvated_box.xyz. High layer atoms are 1-15. Ensure you use the current project settings and don't forget the total charge/multiplicity and QM region charge/multiplicity. It's all neutral singlets.",
        "orca_qmmm_md",
        "hard",
    ),
]

AGENT_RARE_CORPUS: list[CorpusItem] = [
    CorpusItem(
        "td",
        "gaussian",
        "Calculate the first 15 excited singlet and triplet states of OLED_emitter.xyz using TD-DFT in Gaussian. Also account for solvation in THF.",
        "gaussian_tddft_oled",
        "hard",
    ),
    CorpusItem(
        "modred",
        "orca",
        "I need to do a relaxed surface scan using ORCA modredundant for conformer_scan.xyz. Scan the dihedral angle between atoms 4 5 6 7 from 0 to 180 degrees in 18 steps.",
        "orca_modred_scan",
        "hard",
    ),
    CorpusItem(
        "qrc",
        "orca",
        "The TS optimization finished (ts_final.xyz). Now I need to run an ORCA QRC (Quasi-Reaction Coordinate) calculation to verify it connects the correct minimums. Neutral singlet.",
        "orca_qrc_verify",
        "hard",
    ),
    CorpusItem(
        "crest",
        "gaussian",
        "Run a CREST conformational sampling for this highly flexible macrocycle (macrocycle.xyz) before I do DFT optimizations. Use the semiempirical GFN2-xTB level.",
        "gaussian_crest_sampling",
        "hard",
    ),
    CorpusItem(
        "modred",
        "orca",
        "Freeze the distance between the metal center (atom 1) and the ligand nitrogen (atom 15) at 2.1 Angstroms in ORCA, then optimize the rest of the structure. Charge +2, spin 1.",
        "orca_modred_freeze",
        "hard",
    ),
]


AGENT_REPAIR_FODDER: list[CorpusItem] = [
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys1.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_1",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys2.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_2",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys3.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_3",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys4.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_4",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys5.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_5",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys6.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_6",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys7.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_7",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys8.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_8",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys9.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_9",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys10.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_10",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys11.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_11",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys12.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_12",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys13.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_13",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys14.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_14",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys15.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_15",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys16.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_16",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys17.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_17",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys18.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_18",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys19.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_19",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys20.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_20",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys21.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_21",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys22.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_22",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys23.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_23",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys24.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_24",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys25.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_25",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys26.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_26",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys27.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_27",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys28.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_28",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys29.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_29",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys30.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_30",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys31.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_31",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys32.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_32",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys33.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_33",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys34.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_34",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys35.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_35",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys36.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_36",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys37.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_37",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys38.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_38",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys39.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_39",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys40.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_40",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys41.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_41",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys42.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_42",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys43.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_43",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys44.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_44",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys45.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_45",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys46.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_46",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys47.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_47",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys48.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_48",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys49.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 1 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_49",
        "hard",
    ),
    CorpusItem(
        "qmmm",
        "orca",
        "Run ORCA QM/MM optimization on sys50.xyz with QM atoms 1-10. Wait, I forgot to tell you the total charge is 0 and total multiplicity is 1. Figure out the flags.",
        "orca_qmmm_repair_50",
        "hard",
    ),
]

GOAL_BATCH_CORPORA = {}
for batch_name, prompts in BATCHES.items():
    items = []
    for i, p in enumerate(prompts):
        items.append(
            CorpusItem(
                kind="mixed",
                program="mixed",
                request=p,
                workflow=f"{batch_name}_{i}",
                difficulty="hard",
            )
        )
    GOAL_BATCH_CORPORA[batch_name] = items

GAP_FILL_CORPUS: list[CorpusItem] = [
    CorpusItem(
        kind="qrc",
        program="gaussian",
        workflow="qrc_from_ts_log",
        request=(
            "The Gaussian project is already loaded. Give me only the CLI "
            "command to run a quasi-reaction coordinate job from the saddle "
            "point in ts_run.log, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="qrc",
        program="gaussian",
        workflow="qrc_forward_reverse",
        request=(
            "Using the current project settings, synthesize the command that "
            "displaces along the imaginary mode of barrier_qrc.log in both "
            "directions. Command only, no YAML edits. Neutral singlet."
        ),
    ),
    CorpusItem(
        kind="qrc",
        program="orca",
        workflow="orca_qrc_from_out",
        request=(
            "ORCA project is set. Produce only the CLI command for a QRC run "
            "seeded from sn2_saddle.out, charge -1 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="traj",
        program="gaussian",
        workflow="traj_frame_extraction",
        request=(
            "With the loaded Gaussian project, give me just the command to "
            "process the trajectory file md_traj.xyz, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="traj",
        program="gaussian",
        workflow="traj_batch_singlepoint",
        request=(
            "Project YAML is already correct — do not touch it. I need the CLI "
            "command that runs single-point jobs over the last 5 structures of "
            "the trajectory water_traj.xyz, neutral singlet."
        ),
    ),
    # DIAS reads geometries along a reaction coordinate (IRC log by default)
    # and -i/--fragment-indices is a required option.
    CorpusItem(
        kind="dias",
        program="gaussian",
        workflow="distortion_interaction_irc",
        request=(
            "Using the active Gaussian project, synthesize only the command for "
            "a distortion/interaction analysis along cycloaddition_irc.log, "
            "fragment one is atoms 1-6, sample every 5th point, charge 0 "
            "multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="dias",
        program="gaussian",
        workflow="dias_ts_mode",
        request=(
            "Command only, current project: run DIAS in TS mode on "
            "da_adduct_ts.log where fragment one is atoms 1-6, charge 0 "
            "multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="crest",
        program="gaussian",
        workflow="crest_conformer_search",
        request=(
            "The project is loaded. Emit only the CLI command for a CREST "
            "conformer search on flexible_chain.xyz, neutral singlet."
        ),
    ),
    CorpusItem(
        kind="nci",
        program="gaussian",
        workflow="nci_dimer",
        request=(
            "Do not write YAML. Using the loaded project, give the command for "
            "NCI analysis of the benzene dimer in dimer_nci.xyz, charge 0 "
            "multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="irc",
        program="orca",
        workflow="orca_irc_from_hessian",
        request=(
            "ORCA project already configured. Only the command please: IRC from "
            "orca_ts.xyz in both directions, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="modred",
        program="orca",
        workflow="orca_frozen_distance",
        request=(
            "Use the loaded ORCA project. Command only: optimize "
            "constrained.xyz while freezing the distance between atoms 2 and 5. "
            "Charge 0, multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="scan",
        program="orca",
        workflow="orca_relaxed_scan",
        request=(
            "ORCA settings are loaded. Synthesize only the CLI command for a "
            "relaxed scan of the 1-4 distance in scan_target.xyz from 1.4 to "
            "2.4 angstrom in 11 steps, neutral singlet."
        ),
    ),
    CorpusItem(
        kind="sp",
        program="orca",
        workflow="orca_sp_on_anion",
        request=(
            "Current ORCA project is fine. Give me the single-point energy "
            "command for anion_sp.xyz, charge -1 multiplicity 1. Command only."
        ),
    ),
    CorpusItem(
        kind="ts",
        program="orca",
        workflow="orca_ts_optimization",
        request=(
            "Project already set up. Emit only the ORCA transition-state "
            "optimization command for saddle_guess.xyz, charge 0 mult 1."
        ),
    ),
    CorpusItem(
        kind="neb",
        program="orca",
        workflow="orca_neb_endpoints",
        request=(
            "Do not create project YAML. Using current ORCA settings, give the "
            "NEB command from reactant_neb.xyz to product_neb.xyz with 8 "
            "images, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="resp",
        program="gaussian",
        workflow="resp_cationic_ligand",
        request=(
            "Loaded Gaussian project. Command only: RESP charge fitting for "
            "cation_resp.xyz, charge +1 multiplicity 1."
        ),
    ),
]

# Second wave: same coverage targets as GAP_FILL_CORPUS, but each request is
# written in a different register (methods-section prose, lab shorthand,
# Chinese, Korean, submit-instead-of-run) so the rows do not collapse onto one
# query skeleton. Authored by hand; never generate these with a template loop.
GAP_FILL_WIDE_CORPUS: list[CorpusItem] = [
    CorpusItem(
        kind="qrc",
        program="gaussian",
        workflow="qrc_methods_section",
        request=(
            "From the SI: 'Saddle points were verified by displacing along the "
            "imaginary mode and reoptimizing to the adjacent minima.' The "
            "project is loaded; give me only the command for verify_ts.log, "
            "neutral singlet."
        ),
    ),
    CorpusItem(
        kind="qrc",
        program="gaussian",
        workflow="qrc_lab_shorthand",
        request=("qrc off ts_a.log, chg 0 mult 1, current proj, cmd only pls"),
    ),
    CorpusItem(
        kind="qrc",
        program="orca",
        workflow="qrc_chinese",
        request=(
            "ORCA 项目已经配置好了。请只给出命令行：以 saddle_cn.out 为起点做 "
            "quick reaction coordinate 计算，电荷 0，多重度 1。不要修改 YAML。"
        ),
    ),
    CorpusItem(
        kind="qrc",
        program="gaussian",
        workflow="qrc_reviewer_request",
        request=(
            "A reviewer asked us to confirm the saddle point connects the "
            "intended minima. Project is loaded; command only, for "
            "cluster_ts.log, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="traj",
        program="gaussian",
        workflow="traj_conformer_grouping",
        request=(
            "Loaded project. Command only: group the conformers in "
            "md_run_traj.xyz ignoring hydrogens, then run opt on the 3 lowest "
            "structures, neutral singlet."
        ),
    ),
    CorpusItem(
        kind="traj",
        program="gaussian",
        workflow="traj_korean",
        request=(
            "현재 프로젝트 설정 그대로 쓰고 YAML은 건드리지 마. "
            "dynamics_traj.xyz 마지막 10개 구조에 대해 single point 돌리는 "
            "명령어만 알려줘. 전하 0, 다중도 1."
        ),
    ),
    CorpusItem(
        kind="dias",
        program="gaussian",
        workflow="dias_every_n_prose",
        request=(
            "Reviewer wants an activation-strain decomposition. Using the "
            "active project, emit only the command: analyze da_irc.log, "
            "fragment 1 = atoms 1-10, take every 2nd point, charge 0 mult 1."
        ),
    ),
    CorpusItem(
        kind="dias",
        program="gaussian",
        workflow="dias_charged_fragments",
        request=(
            "Command only, project already loaded: distortion/interaction along "
            "anionic_irc.log with fragment one atoms 1-4 carrying charge -1, "
            "fragment two neutral. Overall charge -1, multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="crest",
        program="gaussian",
        workflow="crest_then_opt",
        request=(
            "Use the loaded Gaussian project. Only the CLI command: CREST "
            "conformer search on sugar_flex.xyz, then optimize the 5 "
            "lowest-energy conformers. Charge 0, multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="crest",
        program="gaussian",
        workflow="crest_chinese",
        request=(
            "项目已加载，不要写 YAML。只需要命令：对 peptide_flex.xyz 做构象搜索，"
            "取能量最低的 10 个构象跑单点，电荷 0，多重度 1。"
        ),
    ),
    CorpusItem(
        kind="nci",
        program="gaussian",
        workflow="nci_host_guest",
        request=(
            "Host-guest complex in hg_complex.xyz. Project settings are final. "
            "Give the noncovalent interaction analysis command only, charge 0 "
            "multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="nci",
        program="gaussian",
        workflow="nci_pi_stacking",
        request=(
            "Quantify the pi-stacking contacts in stack_pi.xyz with the current "
            "Gaussian project. Neutral singlet. Command only, no YAML edits."
        ),
    ),
    CorpusItem(
        kind="resp",
        program="gaussian",
        workflow="resp_methods_prose",
        request=(
            "Methods: 'Atomic charges were derived by RESP fitting to the "
            "electrostatic potential.' Project loaded. Command only for "
            "drug_like.xyz, charge 0, multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="wbi",
        program="gaussian",
        workflow="wbi_metal_complex",
        request=(
            "Current project stands. Emit just the command for Wiberg bond "
            "index analysis on pd_complex.xyz, charge +2, multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="irc",
        program="orca",
        workflow="orca_irc_maxpoints",
        request=(
            "ORCA project is configured. Command only: IRC from ts_orca_b.xyz, "
            "forward direction only, charge 0 multiplicity 1."
        ),
    ),
    CorpusItem(
        kind="modred",
        program="orca",
        workflow="orca_modred_two_distances",
        request=(
            "Loaded ORCA project, no YAML changes. Optimize precomplex.xyz with "
            "the 1-3 and 2-4 distances frozen. Charge 0, multiplicity 1. "
            "Command only."
        ),
    ),
    CorpusItem(
        kind="scan",
        program="orca",
        workflow="orca_scan_chinese",
        request=(
            "ORCA 设置已经就绪。只输出命令：扫描 scan_cn.xyz 中原子 2 和 6 之间的"
            "距离，从 1.5 到 3.0 埃，分 16 步。电荷 0，多重度 1。"
        ),
    ),
    CorpusItem(
        kind="sp",
        program="orca",
        workflow="orca_sp_radical",
        request=(
            "ORCA project fine as-is. Single point energy on radical_sp.xyz "
            "with ORCA, charge 0, doublet. Command only please."
        ),
    ),
    CorpusItem(
        kind="ts",
        program="orca",
        workflow="orca_ts_from_guess",
        request=(
            "Refine the saddle-point guess in guess_ts_b.xyz with ORCA using "
            "the loaded project. Charge 0, multiplicity 1. Command only."
        ),
    ),
    CorpusItem(
        kind="neb",
        program="orca",
        workflow="neb_with_ts_guess",
        request=(
            "Do not edit YAML. ORCA NEB from r_b.xyz to p_b.xyz with 12 images "
            "and a TS guess in guess_b.xyz, charge 0 multiplicity 1. Command "
            "only."
        ),
    ),
    CorpusItem(
        kind="opt",
        program="orca",
        workflow="orca_opt_dianion",
        request=(
            "Loaded ORCA project. Command only: optimize dianion.xyz, charge "
            "-2, multiplicity 1."
        ),
    ),
]

AGENTIC_STYLE_CORPUS: list[CorpusItem] = [
    CorpusItem(
        kind="td",
        program="gaussian",
        workflow="methods_section_to_command",
        request=(
            "From my methods note: UV/Vis was computed with Gaussian TD-DFT. "
            "Use the already-loaded project settings and make only the CLI "
            "command for dye_methods.xyz, singlets, 8 states, root 1, neutral "
            "singlet. Do not write or revise project YAML."
        ),
    ),
    CorpusItem(
        kind="scan",
        program="gaussian",
        workflow="lab_notebook_shorthand",
        request=(
            "quick setup: gauss relaxed scan h2o_style.xyz O-H bond 1-2, "
            "12 pts, 0.04 step; keep 1-3 fixed; neutral singlet; use current "
            "project."
        ),
    ),
    CorpusItem(
        kind="modred",
        program="orca",
        workflow="korean_typo_constraint",
        request=(
            "orca로 precomplex_style.xyz 최적화해줘. 형성되는 bond 3-7은 "
            "고정하고, 전하 0 spin singlet. yml 새로 만들지 말고 현재 설정 사용."
        ),
    ),
    CorpusItem(
        kind="neb",
        program="orca",
        workflow="command_vs_yaml_trap",
        request=(
            "The project YAML is already correct. Do not interpret this as a "
            "project setup request. I need the ORCA NEB-TS command for "
            "r_style.xyz to p_style.xyz, 6 images, XTB2 preopt, charge 0 mult 1."
        ),
    ),
    CorpusItem(
        kind="ask",
        program="gaussian",
        workflow="missing_charge_spin_clarification",
        request=(
            "Can you run the TS confirmation for unknown_spin_style.xyz with "
            "Gaussian? The paper only says it was optimized and verified by "
            "one imaginary frequency."
        ),
    ),
    CorpusItem(
        kind="ts",
        program="gaussian",
        workflow="typo_rich_user_query",
        request=(
            "gussian TS opt for proton_style.xyz, max step 8 plz, neutral "
            "singlett. use loaded project and just give the chemsmart command."
        ),
    ),
    CorpusItem(
        kind="wbi",
        program="gaussian",
        workflow="chinese_analysis_request",
        request=(
            "请用当前 Gaussian 项目设置，为 intermediate_style.xyz 生成 Wiberg "
            "bond index 命令，电荷 +1，多重度 1。不要新建 yaml。"
        ),
    ),
    CorpusItem(
        kind="resp",
        program="gaussian",
        workflow="force_field_workflow_request",
        request=(
            "I am preparing MD parameters; make the Gaussian RESP charge command "
            "for ligand_style.xyz, charge -1 multiplicity 1, using current "
            "project settings only."
        ),
    ),
    CorpusItem(
        kind="qmmm",
        program="orca",
        workflow="multiscale_boundary_request",
        request=(
            "For enzyme_active_site_style.xyz, build only the ORCA QM/MM command: "
            "high-level atoms 1-12, total charge -1, multiplicity 1. The project "
            "YAML is already loaded."
        ),
    ),
    CorpusItem(
        kind="sp",
        program="orca",
        workflow="post_dft_energy_refinement",
        request=(
            "After geometry optimization, I only need an ORCA single-point "
            "energy for opt_product_style.xyz, charge 0 mult 1, from the active "
            "project settings."
        ),
    ),
    CorpusItem(
        kind="decline",
        program="gaussian",
        workflow="impossible_without_structure",
        request=(
            "Find the best catalyst for CO2 reduction and run the whole "
            "Gaussian workflow, but I do not have a structure file, charge, "
            "or multiplicity yet."
        ),
    ),
    CorpusItem(
        kind="irc",
        program="gaussian",
        workflow="ts_connectivity_after_failure",
        request=(
            "The TS optimized yesterday; now just create the Gaussian IRC command "
            "for ts_style.xyz in both directions, charge 0 multiplicity 1, using "
            "the loaded project."
        ),
    ),
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", help="DashScope model id.")
    parser.add_argument("--index", type=int, default=0)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--limit", type=int, default=4)
    parser.add_argument(
        "--start-position",
        type=int,
        default=1,
        help="1-based corpus position to start from before index/stride selection.",
    )
    parser.add_argument("--batch-id", default="v17_reasoning")
    parser.add_argument(
        "--task-offset",
        type=int,
        default=0,
        help="Offset added to 1-based corpus indices for globally unique task ids.",
    )
    parser.add_argument(
        "--corpus",
        choices=[
            "batch_8",
            "batch_7",
            "batch_6",
            "batch_5",
            "batch_4",
            "batch_3",
            "batch_2",
            "batch_1",
            "workflow-hard",
            "command-hard",
            "agentic-style",
            "agent-domain",
            "gap-fill",
            "gap-fill-wide",
            "qmmm-deep",
            "rare-deep",
            "repair-fodder",
        ],
        default="agent-domain",
    )
    parser.add_argument(
        "--provider",
        choices=sorted(PROVIDERS),
        default="dashscope",
        help="Which api.env credential and endpoint to use.",
    )
    parser.add_argument(
        "--budget-usd",
        type=float,
        default=0.0,
        help="Stop once metered spend reaches this cap (0 disables the guard).",
    )
    args = parser.parse_args(argv)

    os.environ["PATH"] = "/opt/anaconda3/envs/chemsmart/bin:" + os.environ.get(
        "PATH", ""
    )
    sys.path.insert(0, str(REPO))
    os.chdir(REPO)
    logging.disable(logging.INFO)

    key_name, base_url = PROVIDERS[args.provider]
    key = dotenv_values(REPO / "api.env").get(key_name, "")
    if not key:
        raise SystemExit(f"api.env is missing {key_name}")
    meter = CostMeter(args.model, args.budget_usd or float("inf"))

    run_dir = REPO / "var" / "agent-training" / "runs" / args.model
    os.environ["CHEMSMART_AGENT_TRAINING_DIR"] = str(run_dir)

    import chemsmart.agent.providers as providers_mod
    import chemsmart.agent.tools_command as tools_command
    from chemsmart.agent.core import AgentSession
    from chemsmart.agent.permissions import (
        ApprovalDecision,
        PermissionMode,
        PermissionPolicy,
    )
    from chemsmart.agent.providers import OpenAIProvider
    from chemsmart.agent.registry import ToolRegistry

    class ReasoningProvider(OpenAIProvider):
        """Long-timeout provider; do not disable thinking for thinking models."""

        def chat(
            self,
            messages: list[Any],
            tools: list[Any] | None = None,
            timeout_s: float = 30,
        ) -> dict[str, Any]:
            kwargs: dict[str, Any] = {
                "model": self.default_model,
                "messages": messages,
                "timeout": 300.0,
            }
            if tools:
                kwargs["tools"] = tools
            meter.check()
            payload = self._client.chat.completions.create(
                **kwargs
            ).model_dump()
            meter.record(payload.get("usage"))
            return payload

    provider = ReasoningProvider(key, model=args.model, base_url=base_url)
    providers_mod.get_provider = lambda *a, **k: provider

    corpus = {
        "workflow-hard": WORKFLOW_HARD_CORPUS,
        "command-hard": COMMAND_HARD_CORPUS,
        "agentic-style": AGENTIC_STYLE_CORPUS,
        "agent-domain": AGENT_DOMAIN_CORPUS,
        "qmmm-deep": AGENT_QMMM_CORPUS,
        "rare-deep": AGENT_RARE_CORPUS,
        "repair-fodder": AGENT_REPAIR_FODDER,
        "batch_1": GOAL_BATCH_CORPORA["batch_1"],
        "batch_2": GOAL_BATCH_CORPORA["batch_2"],
        "batch_3": GOAL_BATCH_CORPORA["batch_3"],
        "batch_4": GOAL_BATCH_CORPORA["batch_4"],
        "batch_5": GOAL_BATCH_CORPORA["batch_5"],
        "batch_6": GOAL_BATCH_CORPORA["batch_6"],
        "batch_7": GOAL_BATCH_CORPORA["batch_7"],
        "batch_8": GOAL_BATCH_CORPORA["batch_8"],
        "gap-fill": GAP_FILL_CORPUS,
        "gap-fill-wide": GAP_FILL_WIDE_CORPUS,
    }[args.corpus]
    requests = [
        (position, item)
        for position, item in enumerate(corpus, start=1)
        if position >= args.start_position
        and (position - 1) % args.stride == args.index
    ]
    if args.limit:
        requests = requests[: args.limit]
    plan_path = _write_collection_plan(
        run_dir=run_dir,
        model=args.model,
        batch_id=args.batch_id,
        corpus_name=args.corpus,
        selected=requests,
        task_offset=args.task_offset,
    )
    print(f"collection_plan={plan_path}", flush=True)

    results: list[dict[str, Any]] = []
    for position, item in requests:
        task_index = args.task_offset + position
        task_id = f"{args.corpus}-{task_index:03d}"
        record = _run_one(
            task_index=task_index,
            task_id=task_id,
            corpus_name=args.corpus,
            item=item,
            model=args.model,
            batch_id=args.batch_id,
            provider=provider,
            tools_command=tools_command,
            AgentSession=AgentSession,
            ToolRegistry=ToolRegistry,
            PermissionPolicy=PermissionPolicy,
            PermissionMode=PermissionMode,
            ApprovalDecision=ApprovalDecision,
        )
        record["cost_usd_cumulative"] = round(meter.cost_usd, 4)
        results.append(record)
        print(
            f"{args.model:32s} {task_id:18s} {item.kind:7s} {record['grade']:7s} "
            f"rsn={record.get('reasoning_len', 0):5d} "
            f"{record.get('latency_s', 0):6.1f}s "
            f"${meter.cost_usd:7.4f} "
            f"{str(record.get('command') or record.get('error') or '')[:60]}",
            flush=True,
        )
        if record["grade"] in {"QUOTA", "BUDGET"}:
            break
        time.sleep(1)

    _append_results(run_dir, results)
    tally: dict[str, int] = {}
    for row in results:
        tally[row["grade"]] = tally.get(row["grade"], 0) + 1
    with_reasoning = sum(1 for row in results if row.get("reasoning_len", 0))
    print(
        f"== {args.model}: {tally} | with_reasoning={with_reasoning}/{len(results)} ==",
        flush=True,
    )
    print(
        f"== spend: {json.dumps(meter.summary(), sort_keys=True)} ==",
        flush=True,
    )
    return 0


def _run_one(
    *,
    task_index: int,
    task_id: str,
    corpus_name: str,
    item: CorpusItem,
    model: str,
    batch_id: str,
    provider: Any,
    tools_command: Any,
    AgentSession: Any,
    ToolRegistry: Any,
    PermissionPolicy: Any,
    PermissionMode: Any,
    ApprovalDecision: Any,
) -> dict[str, Any]:
    work = Path(tempfile.mkdtemp(prefix=f"rz-{task_id}-{model[:8]}-"))
    _prepare_workspace(work, item.request, item.program)
    original_cwd = Path.cwd()
    os.chdir(work)
    tools_command.reset_command_tools_state()
    record: dict[str, Any] = {
        "batch_id": batch_id,
        "corpus": corpus_name,
        "task_id": task_id,
        "task_index": task_index,
        "model": model,
        "kind": item.kind,
        "program": item.program,
        "workflow": item.workflow,
        "difficulty": item.difficulty,
        "request": item.request[:160],
    }
    try:
        session = AgentSession(
            provider=provider,
            registry=ToolRegistry.default(
                groups=["synthesis", "project_yaml"]
            ),
            session_root=work / "sessions" / task_id,
            stage_prompt="unified_agent.md",
        )
        policy = PermissionPolicy(
            mode=PermissionMode.DRIVING,
            prompt_risky=True,
        )
        started = time.time()
        out = session.run_loop(
            item.request,
            policy=policy,
            approver=lambda req: ApprovalDecision.ALLOW_ONCE,
        )
        record["latency_s"] = round(time.time() - started, 1)
        synth = _last_synthesis(out)
        if synth:
            command = str(synth.get("command") or "").strip()
            semantic = synth.get("semantic") or {}
            verdict = (
                semantic.get("verdict") if isinstance(semantic, dict) else None
            )
            record["command"] = command
            record["gate"] = verdict
            record["reasoning_len"] = len(synth.get("reasoning") or "")
            status = synth.get("status")
            if (
                status == "ready"
                and verdict in {"ok", "warn"}
                and _kind_is_present(item.kind, command)
            ):
                record["grade"] = "PASS"
            elif status == "needs_clarification":
                record["grade"] = "ASK"
            else:
                record["grade"] = "WRONG"
        else:
            record["grade"] = "NOSYNTH"
    except Exception as exc:  # pragma: no cover - live API path
        message = f"{type(exc).__name__}: {exc}"[:240]
        record["error"] = message
        lowered = message.lower()
        if "reached cap" in lowered:
            record["grade"] = "BUDGET"
        elif "free quota" in lowered:
            record["grade"] = "QUOTA"
        else:
            record["grade"] = "ERR"
    finally:
        os.chdir(original_cwd)
    return record


def _write_collection_plan(
    *,
    run_dir: Path,
    model: str,
    batch_id: str,
    corpus_name: str,
    selected: list[tuple[int, CorpusItem]],
    task_offset: int,
) -> Path:
    """Persist task indices before any API request is sent."""

    plan_dir = run_dir / "collection_plans"
    plan_dir.mkdir(parents=True, exist_ok=True)
    suffix = int(time.time())
    path = plan_dir / f"{batch_id}_{corpus_name}_{suffix}.jsonl"
    with path.open("w", encoding="utf-8") as handle:
        for position, item in selected:
            task_index = task_offset + position
            row = {
                "batch_id": batch_id,
                "corpus": corpus_name,
                "difficulty": item.difficulty,
                "kind": item.kind,
                "model": model,
                "program": item.program,
                "request": item.request,
                "task_id": f"{corpus_name}-{task_index:03d}",
                "task_index": task_index,
                "workflow": item.workflow,
                "workspace_yaml_policy": "one project YAML per task workspace",
            }
            handle.write(
                json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n"
            )
    return path


def _prepare_workspace(work: Path, request: str, program: str) -> None:
    # Keep exactly one workspace project visible. Multiple project YAML files
    # intentionally trigger yaml_check clarification and would poison command
    # synthesis collection with project-selection failures.
    folder = work / ".chemsmart" / program
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "demo.yaml").write_text(_project_yaml(program), encoding="utf-8")
    xyz_text = _placeholder_xyz(_required_atom_count(request))
    for filename in sorted(set(_files_in(request))):
        _write_input_fixture(work / filename, xyz_text, program)


def _write_input_fixture(path: Path, xyz_text: str, program: str) -> None:
    """Create source-valid QRC and trajectory inputs for collection runs."""

    suffix = path.suffix.lower()
    if suffix in {".log", ".out"}:
        source = (
            REPO / "tests/data/GaussianTests/outputs/Pd_insertion_ts_r.log"
            if program == "gaussian"
            else REPO / "tests/data/ORCATests/outputs/sn2_ts.out"
        )
        shutil.copyfile(source, path)
        return
    if "traj" in path.stem.lower():
        source = REPO / "tests/data/ORCATests/outputs/water_opt_trj.xyz"
        shutil.copyfile(source, path)
        return
    path.write_text(xyz_text, encoding="utf-8")


def _project_yaml(program: str) -> str:
    if program == "orca":
        return (
            "gas:\n"
            "  functional: b3lyp\n"
            "  basis: def2svp\n"
            "  freq: true\n"
            "solv:\n"
            "  functional: b3lyp\n"
            "  basis: def2svp\n"
            "  freq: false\n"
            "  solvent_model: smd\n"
            "  solvent_id: water\n"
            "qmmm:\n"
            "  qm_functional: b3lyp\n"
            "  qm_basis: def2svp\n"
            "  qm2_functional: pbe\n"
            "  qm2_basis: 6-31g(d,p)/auto\n"
            "  mm_force_field: AMBER=HardFirst\n"
        )
    return (
        "gas:\n"
        "  functional: b3lyp\n"
        "  basis: def2svp\n"
        "  freq: true\n"
        "solv:\n"
        "  functional: b3lyp\n"
        "  basis: def2svp\n"
        "  freq: false\n"
        "  solvent_model: smd\n"
        "  solvent_id: water\n"
        "qmmm:\n"
        "  high_level_functional: b3lyp\n"
        "  high_level_basis: def2svp\n"
        "  medium_level_functional: pbe\n"
        "  medium_level_basis: 6-31g(d,p)/auto\n"
        "  low_level_force_field: AMBER=HardFirst\n"
    )


def _files_in(text: str) -> list[str]:
    return re.findall(r"[\w가-힣+-]+\.(?:xyz|log|out)", text)


def _required_atom_count(text: str) -> int:
    """Return a fixture atom count compatible with indexed requests.

    Indexed constraints appear in many surface forms — "atoms 3 and 11",
    "atoms 4, 7 and 9", "atoms 1 through 6", "dihedral 5 6 7 8". The fixture
    must contain the highest referenced index or the runtime raises an
    IndexError on a too-small placeholder. Over-provisioning is harmless, so
    every integer inside an atom-index clause is collected and the max taken.
    """

    maxima = [3]
    for left, right in re.findall(r"\b(\d+)\s*-\s*(\d+)\b", text):
        maxima.append(int(left))
        maxima.append(int(right))
    # An index clause opens with an index keyword and runs to sentence end;
    # numbers that clearly belong to non-index quantities are trimmed off so
    # "12 steps" or "charge 0" cannot inflate or distort the count.
    clause_re = re.compile(
        r"\b(?:atoms?|bond|distance|angle|dihedral|torsion|fragment|"
        r"indices|index|high[- ]?level|qm)\b[^.:;\n]*",
        re.I,
    )
    stop_re = re.compile(
        r"\b(?:charge|multiplicit\w*|mult|spin|step|steps|state|states|"
        r"image|images|degree|degrees|root|angstrom\w*|conformer\w*|"
        r"point|points|amplitude|mode)\b",
        re.I,
    )
    for match in clause_re.finditer(text):
        clause = stop_re.split(match.group(0))[0]
        for value in re.findall(r"\b(\d+)\b", clause):
            maxima.append(int(value))
    return max(max(maxima), 3)


def _placeholder_xyz(atom_count: int) -> str:
    if atom_count <= 3:
        return H2O_XYZ
    lines = [str(atom_count), f"{atom_count}-atom indexed placeholder"]
    for index in range(atom_count):
        x = 1.35 * (index % 5)
        y = 1.20 * ((index // 5) % 3)
        z = 1.10 * (index // 15)
        element = "C" if index % 3 else "N"
        lines.append(f"{element} {x:.6f} {y:.6f} {z:.6f}")
    return "\n".join(lines) + "\n"


def _last_synthesis(loop_output: dict[str, Any]) -> dict[str, Any] | None:
    outcomes = loop_output.get("tool_outcomes") or []
    for outcome in reversed(outcomes):
        if getattr(outcome, "name", "") != "synthesize_command":
            continue
        raw = getattr(outcome, "raw_result", None)
        if isinstance(raw, dict):
            return raw
    return None


def _kind_is_present(kind: str, command: str) -> bool:
    if kind == "ask":
        return False
    return kind in command.split()


def _append_results(run_dir: Path, results: list[dict[str, Any]]) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / "reasoning_accum_results.jsonl"
    with path.open("a", encoding="utf-8") as handle:
        for row in results:
            handle.write(
                json.dumps(row, sort_keys=True, ensure_ascii=False) + "\n"
            )


if __name__ == "__main__":
    raise SystemExit(main())
