"""Accumulate hard multi-turn and repair ChemSmart agent episodes.

This live-data helper complements ``reasoning_accum.py``.  It targets the
trajectory types that are underrepresented in the raw store: ask/resume,
command modification, repair_command, project-YAML-to-command chains, and
safe refusal/clarification.  Each selected scenario is indexed and written to a
collection plan before any provider call.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import signal
import sys
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import reasoning_accum
from dotenv import dotenv_values

REPO = Path(__file__).resolve().parents[2]
BASE_URL = reasoning_accum.BASE_URL

# Bumped when this collector's scenario schema or fixture conventions change;
# recorded as dataset_provenance.collector_version on every episode it drives.
COLLECTOR_VERSION = "agentic_workflow_accum@1"


def _cli_schema_version() -> str:
    try:
        from chemsmart import __version__

        return str(__version__)
    except Exception:
        return "unknown"


def _scenario_provenance(
    scenario: "Scenario", *, batch_id: str, model: str
) -> dict[str, str]:
    """Sanitized-by-construction provenance for the episodes this run drives.

    ``scenario_family`` groups every multi-turn variant of one workflow so the
    split builder keeps them on the same side of train/eval; ``fixture_id``
    fingerprints the prepared workspace (derived from the joined turns).
    """

    fixture_id = hashlib.sha256(
        " ".join(scenario.turns).encode("utf-8")
    ).hexdigest()[:12]
    family = (
        scenario.workflow or f"{scenario.program}.{scenario.expected_kind}"
    )
    return {
        "scenario_id": scenario.name,
        "scenario_family": family,
        "batch_id": batch_id,
        "teacher_id": model,
        "fixture_id": fixture_id,
        "schema_version": _cli_schema_version(),
        "collector_version": COLLECTOR_VERSION,
    }


class ScenarioTimeout(RuntimeError):
    """Raised when a live teacher scenario exceeds its collection budget."""


@contextmanager
def _scenario_timeout(timeout_s: float):
    """Bound one scenario so repair loops cannot exhaust free quota."""

    if timeout_s <= 0:
        yield
        return

    previous = signal.getsignal(signal.SIGALRM)

    def _raise_timeout(signum, frame):
        del signum, frame
        raise ScenarioTimeout(f"scenario exceeded {timeout_s:g}s")

    signal.signal(signal.SIGALRM, _raise_timeout)
    signal.setitimer(signal.ITIMER_REAL, timeout_s)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


@dataclass(frozen=True)
class Scenario:
    name: str
    program: str
    expected_kind: str
    workflow: str
    turns: tuple[str, ...]
    difficulty: str = "hard"
    source_reference: str = "ChemSmart CLI and runtime semantic contract"


SCENARIOS: list[Scenario] = [
    Scenario(
        name="ask_resume_ts_charge_spin",
        program="gaussian",
        expected_kind="ts",
        workflow="ask_resume_missing_charge_spin",
        turns=(
            "I need Gaussian TS confirmation for unknown_spin_mt.xyz; the paper only says it has one imaginary frequency.",
            "Use charge 0 and multiplicity 1. Keep the loaded project and make the command only.",
        ),
    ),
    Scenario(
        name="repair_ts_maxstep_not_maxcycle",
        program="gaussian",
        expected_kind="ts",
        workflow="repair_command_semantic_slot",
        turns=(
            "Repair this failed command. Intent: Gaussian TS optimization for proton_repair.xyz, neutral singlet, maxstep 8. Bad command: chemsmart run gaussian -p demo -f proton_repair.xyz -c 0 -m 1 -r maxcycle=8 ts",
        ),
    ),
    Scenario(
        name="modify_scan_add_freeze",
        program="gaussian",
        expected_kind="scan",
        workflow="multi_turn_command_modification",
        turns=(
            "Using the current Gaussian project, make a relaxed scan command for scan_mt.xyz: bond 1-2, 8 steps, 0.06 step, neutral singlet.",
            "Change it: also freeze bond 1-3 and keep the same scan settings.",
        ),
    ),
    Scenario(
        name="repair_orca_neb_route",
        program="orca",
        expected_kind="neb",
        workflow="repair_orca_neb_command",
        turns=(
            "Fix this ORCA command. Intent: NEB-TS from neb_r_mt.xyz to neb_p_mt.xyz with 7 images and XTB2 pre-optimization, charge 0 multiplicity 1. Bad command: chemsmart run orca -p demo -f neb_r_mt.xyz -c 0 -m 1 neb -e neb_p_mt.xyz --nimages 7 -j NEB",
        ),
    ),
    Scenario(
        name="project_yaml_then_opt_command",
        program="gaussian",
        expected_kind="opt",
        workflow="project_yaml_to_command_chain",
        turns=(
            "Create a project YAML named co2chain.yaml for Gaussian gas-phase B3LYP-D3BJ/def2-SVP geometry optimization and harmonic frequency confirmation. CREST/GFN2-xTB conformer search was done externally; include it only as protocol context.",
            "Now use that project to make a Gaussian optimization command for co2_adduct_mt.xyz, charge 0 multiplicity 1.",
        ),
    ),
    Scenario(
        name="basis_correction_then_command",
        program="gaussian",
        expected_kind="opt",
        workflow="project_yaml_correction_chain",
        turns=(
            "Build a Gaussian project YAML named brchain.yaml for B3LYP-D3BJ with def2-SVP, gas phase, opt plus freq.",
            "Correction: Br should use def2-SVPD while light atoms stay def2-SVP. Fix the project YAML and then make the opt command for br_complex_mt.xyz, charge -1 multiplicity 1.",
        ),
    ),
    Scenario(
        name="impossible_decline_no_structure",
        program="gaussian",
        expected_kind="decline",
        workflow="safe_refusal_missing_structure",
        turns=(
            "Find the best catalyst for CO2 reduction and run the full Gaussian workflow, but I do not have a structure file, charge, or multiplicity yet.",
        ),
    ),
    Scenario(
        name="orca_qmmm_boundary_repair",
        program="orca",
        expected_kind="qmmm",
        workflow="hard_multiscale_boundary",
        turns=(
            "Use ORCA QM/MM for enzyme_boundary_mt.xyz. High-level atoms are 1-12, total charge -1, multiplicity 1. If the command cannot be formed, ask for exactly the missing runtime slot.",
        ),
    ),
    Scenario(
        name="tddft_eqsolv_state_correction",
        program="gaussian",
        expected_kind="td",
        workflow="multi_turn_spectroscopy_correction",
        turns=(
            "Set up Gaussian TD-DFT absorption for dye_td_mt.xyz, neutral singlet, 12 singlet states, root 1, and use equilibrium solvation.",
            "Correction: make it 8 triplet states instead, but keep eqsolv and the same project.",
        ),
    ),
    Scenario(
        name="dias_fragment_semantic_repair",
        program="gaussian",
        expected_kind="dias",
        workflow="fragment_index_semantic_repair",
        turns=(
            "Repair this DIAS command. Intent: Gaussian DIAS for diels_alder_mt.xyz with fragment 1 atoms 1-6 and fragment 2 atoms 7-12, charge 0 multiplicity 1. Bad command: chemsmart run gaussian -p demo -f diels_alder_mt.xyz -c 0 -m 1 dias --fragment-indices '[[1,2,3,4,5,6],[7,8,9,10,11,12]]'",
        ),
    ),
    Scenario(
        name="wbi_no_dias_fragments",
        program="gaussian",
        expected_kind="wbi",
        workflow="bond_order_not_dias",
        turns=(
            "I only need Wiberg bond indices for ligand_wbi_mt.xyz with Gaussian, charge 0 multiplicity 1. Do not make this a DIAS job.",
        ),
    ),
    Scenario(
        name="modred_internal_freeze_not_atom_freeze",
        program="gaussian",
        expected_kind="modred",
        workflow="internal_coordinate_disambiguation",
        turns=(
            "Constrain the internal bond distance between atoms 2 and 5 during a Gaussian optimization of bridge_mt.xyz, neutral singlet. This is a modredundant bond constraint, not freezing atoms in Cartesian space.",
        ),
    ),
    Scenario(
        name="orca_ts_trust_radius_recalc",
        program="orca",
        expected_kind="ts",
        workflow="orca_ts_option_preservation",
        turns=(
            "Prepare an ORCA transition-state optimization for organocat_ts_mt.xyz, charge 0 multiplicity 1, using recalc_hess 5 and trust radius 0.12.",
        ),
    ),
    Scenario(
        name="gaussian_irc_after_ts",
        program="gaussian",
        expected_kind="irc",
        workflow="post_ts_reaction_path",
        turns=(
            "From ts_guess_mt.xyz, set up a Gaussian IRC calculation in both directions for a neutral singlet transition state.",
        ),
    ),
    Scenario(
        name="resp_charge_fit_from_optimized_structure",
        program="gaussian",
        expected_kind="resp",
        workflow="charge_derivation_workflow",
        turns=(
            "Prepare a Gaussian RESP charge fitting setup for optimized_ligand_mt.xyz, charge -1 multiplicity 1. This should be a charge derivation job, not a geometry optimization.",
        ),
    ),
    Scenario(
        name="decline_missing_product_for_neb",
        program="orca",
        expected_kind="decline",
        workflow="safe_refusal_missing_endpoint",
        turns=(
            "Run an ORCA NEB-TS calculation for reactant_only_mt.xyz but I do not have a product endpoint structure yet.",
        ),
    ),
    Scenario(
        name="neb_missing_product_then_resume",
        program="orca",
        expected_kind="neb",
        workflow="ask_resume_missing_neb_endpoint",
        turns=(
            "I want an ORCA NEB-TS setup from reactant_resume_mt.xyz, 9 images, charge 0 multiplicity 1, but I forgot to give the product endpoint.",
            "Use product_resume_mt.xyz as the product endpoint. Keep 9 images and make the command.",
        ),
    ),
    Scenario(
        name="korean_ts_maxstep_repair",
        program="gaussian",
        expected_kind="ts",
        workflow="korean_route_repair",
        turns=(
            "Gaussian TS 최적화 명령을 고쳐줘. 의도는 kr_ts_mt.xyz, neutral singlet, maxstep 8이야. 잘못된 명령: chemsmart run gaussian -p demo -f kr_ts_mt.xyz -c 0 -m 1 -r maxcycle=8 ts",
            "maxcycle 말고 TS opt route 안에 maxstep=8이 들어가야 해. 다시 repair해.",
        ),
    ),
    Scenario(
        name="scan_change_coordinate_type",
        program="gaussian",
        expected_kind="scan",
        workflow="multi_turn_coordinate_correction",
        turns=(
            "Set up a Gaussian relaxed scan for ring_open_mt.xyz, neutral singlet, scan bond 3-4 for 6 steps of 0.04.",
            "Actually scan the angle 2-3-4 instead, still 6 steps and 0.04, and freeze bond 1-5.",
        ),
    ),
    Scenario(
        name="yaml_solvation_then_tddft",
        program="gaussian",
        expected_kind="td",
        workflow="project_yaml_update_to_spectroscopy",
        turns=(
            "Create a Gaussian project YAML named solvtd.yaml using CAM-B3LYP/def2-TZVP with SMD acetonitrile for excited-state calculations.",
            "Now make a TD-DFT command for chromophore_solv_mt.xyz, charge 0 multiplicity 1, 15 singlet states and root 2 using that project.",
        ),
    ),
    Scenario(
        name="basis_ambiguous_karlsruhe_then_command",
        program="orca",
        expected_kind="sp",
        workflow="basis_search_then_command",
        turns=(
            "For ORCA, I need a single point on nickel_cat_mt.xyz with the Karlsruhe triple-zeta basis and a matching RI/J auxiliary basis, charge 0 multiplicity 1.",
            "Use the best matching concrete basis names and make the command if the project settings are available; otherwise ask only for the missing project slot.",
        ),
    ),
    Scenario(
        name="dias_then_wbi_correction",
        program="gaussian",
        expected_kind="wbi",
        workflow="kind_correction_dias_to_wbi",
        turns=(
            "Set up Gaussian DIAS-like bond analysis for complex_wbi_mt.xyz, charge 0 multiplicity 1.",
            "Correction: I do not need DIAS fragments. I only need Wiberg bond indices, so switch the command to WBI.",
        ),
    ),
    Scenario(
        name="decline_no_charge_spin_then_resume",
        program="gaussian",
        expected_kind="opt",
        workflow="ask_resume_missing_charge_mult",
        turns=(
            "Optimize radical_candidate_mt.xyz with Gaussian using the current project, but I am not sure about charge and multiplicity.",
            "Use charge 1 and multiplicity 2. Now make the Gaussian optimization command.",
        ),
    ),
    Scenario(
        name="orca_aux_basis_flag_repair",
        program="orca",
        expected_kind="sp",
        workflow="upstream_aux_basis_flag_repair",
        turns=(
            "Repair this ORCA command. Intent: DLPNO single point for aux_repair_mt.xyz, charge 0 multiplicity 1, with def2/J auxiliary basis. Bad command: chemsmart run orca -p demo -f aux_repair_mt.xyz -c 0 -m 1 -a def2/J sp",
        ),
    ),
    Scenario(
        name="tddft_missing_states_then_resume",
        program="gaussian",
        expected_kind="td",
        workflow="ask_resume_missing_excited_state_slots",
        turns=(
            "Prepare a Gaussian TD-DFT job for push_pull_td2_mt.xyz, neutral singlet, but I forgot how many excited states and which spin manifold.",
            "Use 20 singlet states, root 3, and keep equilibrium solvation if the active project supports solvation.",
        ),
    ),
    Scenario(
        name="orca_neb_change_images_and_product",
        program="orca",
        expected_kind="neb",
        workflow="multi_turn_neb_endpoint_and_nimages_change",
        turns=(
            "Set up ORCA NEB-TS from neb_a_mt.xyz to neb_b_mt.xyz, neutral singlet, 5 images.",
            "Correction: use neb_c_mt.xyz as the product endpoint instead and increase to 11 images.",
        ),
    ),
    Scenario(
        name="gaussian_irc_direction_correction",
        program="gaussian",
        expected_kind="irc",
        workflow="multi_turn_irc_direction_correction",
        turns=(
            "Set up Gaussian IRC from saddle_irc_mt.xyz, neutral singlet, only forward direction.",
            "Actually make it reverse direction instead, keeping the same structure and project.",
        ),
    ),
    Scenario(
        name="resp_wrong_opt_to_resp_correction",
        program="gaussian",
        expected_kind="resp",
        workflow="kind_correction_opt_to_resp",
        turns=(
            "Optimize resp_candidate_mt.xyz with Gaussian, charge 0 multiplicity 1, then derive RESP charges.",
            "Correction: do not optimize here. I only want the RESP charge fitting command from the already optimized structure.",
        ),
    ),
    Scenario(
        name="korean_neb_endpoint_resume",
        program="orca",
        expected_kind="neb",
        workflow="korean_ask_resume_missing_endpoint",
        turns=(
            "ORCA NEB-TS 작업을 만들고 싶어. reactant는 kr_neb_r_mt.xyz, 이미지 7개, 전하 0, 다중도 1인데 product 파일을 아직 안 줬어.",
            "product는 kr_neb_p_mt.xyz야. 같은 조건으로 chemsmart 명령을 만들어줘.",
        ),
    ),
    Scenario(
        name="project_yaml_dispersion_correction_then_ts",
        program="gaussian",
        expected_kind="ts",
        workflow="project_yaml_correction_then_ts_command",
        turns=(
            "Create a Gaussian project YAML named tsproto.yaml for B3LYP-D3/def2-SVP gas-phase transition-state searches.",
            "Correction: use D3BJ, not plain D3. Then prepare a TS optimization command for ts_proto_mt.xyz, charge 0 multiplicity 1, maxstep 6.",
        ),
    ),
    Scenario(
        name="wbi_then_add_dias_fragments_correction",
        program="gaussian",
        expected_kind="dias",
        workflow="kind_correction_wbi_to_dias",
        turns=(
            "Make a Gaussian WBI command for frag_rxn_mt.xyz, charge 0 multiplicity 1.",
            "Correction: I actually need DIAS with fragment 1 atoms 1-4 and fragment 2 atoms 5-9. Ask only if the fragment syntax cannot be represented safely.",
        ),
    ),
    Scenario(
        name="orca_ts_opt_not_plain_opt",
        program="orca",
        expected_kind="ts",
        workflow="kind_correction_opt_to_orca_ts",
        turns=(
            "Prepare an ORCA optimization for orca_saddle_mt.xyz, charge -1 multiplicity 2.",
            "Correction: this is a first-order saddle search, so use ORCA transition-state optimization with recalc_hess 10.",
        ),
    ),
    Scenario(
        name="gaussian_scan_missing_step_then_resume",
        program="gaussian",
        expected_kind="scan",
        workflow="ask_resume_missing_scan_step",
        turns=(
            "Set up a relaxed Gaussian scan for scan_missing_step_mt.xyz over bond 1-2, neutral singlet. I know it should be 12 steps but forgot the increment.",
            "Use increment 0.08 and also freeze bond 1-3 during the scan.",
        ),
    ),
    Scenario(
        name="orca_aux_basis_project_update_then_sp",
        program="orca",
        expected_kind="sp",
        workflow="project_yaml_aux_basis_then_command",
        turns=(
            "Create an ORCA project YAML named auxfit.yaml for PBE0/def2-TZVP single points with a def2/J RI auxiliary basis.",
            "Now make the ORCA single-point command for auxfit_substrate_mt.xyz, charge 0 multiplicity 1 using that project.",
        ),
    ),
    Scenario(
        name="decline_impossible_literature_only_then_file",
        program="gaussian",
        expected_kind="opt",
        workflow="decline_then_resume_with_file",
        turns=(
            "Use the molecule described in the paper text and run a Gaussian optimization; I do not have coordinates or a structure file.",
            "I found the structure file: paper_molecule_mt.xyz. It is neutral singlet. Now prepare the optimization command.",
        ),
    ),
    Scenario(
        name="gaussian_ts_calcall_variant_repair",
        program="gaussian",
        expected_kind="ts",
        workflow="route_variant_repair_calcall",
        turns=(
            "Repair this TS command. Intent: Gaussian TS optimization of calcall_ts_mt.xyz, neutral singlet, using calcall instead of calcfc. Bad command: chemsmart run gaussian -p demo -f calcall_ts_mt.xyz -c 0 -m 1 -r 'calcfc calcall' ts",
            "Keep calcall, remove calcfc conflict, and preserve noeigentest in the generated TS route.",
        ),
    ),
    Scenario(
        name="tddft_eqsolv_to_none_correction",
        program="gaussian",
        expected_kind="td",
        workflow="multi_turn_tddft_solvation_correction",
        turns=(
            "Prepare Gaussian TD-DFT for push_pull_solvent_mt.xyz, charge 0 multiplicity 1, 10 singlet states, root 2, with equilibrium solvation.",
            "Correction: remove equilibrium solvation and make it 14 triplet states, root 1. Keep the same input file and project.",
        ),
    ),
    Scenario(
        name="modred_not_cartesian_freeze_correction",
        program="gaussian",
        expected_kind="modred",
        workflow="kind_correction_atom_freeze_to_modred",
        turns=(
            "Optimize clamp_bridge_mt.xyz with Gaussian, neutral singlet, and keep atoms 2 and 5 fixed.",
            "Correction: I meant constrain the internal bond distance between atoms 2 and 5 during optimization, not freeze Cartesian atoms.",
        ),
    ),
    Scenario(
        name="scan_angle_then_dihedral_correction",
        program="gaussian",
        expected_kind="scan",
        workflow="multi_turn_scan_coordinate_rank_change",
        turns=(
            "Set up a relaxed Gaussian scan for torsion_probe_mt.xyz, neutral singlet, scanning angle 2-3-4 for 9 steps of 0.03.",
            "Actually scan the dihedral 1-2-3-4 for 12 steps of 10 degrees, and freeze bond 5-6.",
        ),
    ),
    Scenario(
        name="orca_neb_missing_images_resume",
        program="orca",
        expected_kind="neb",
        workflow="ask_resume_missing_neb_nimages",
        turns=(
            "Create an ORCA NEB-TS command from epoxide_r_mt.xyz to epoxide_p_mt.xyz, neutral singlet, but I have not decided the number of images.",
            "Use 13 images and keep it as NEB-TS. Do not switch this to a plain optimization.",
        ),
    ),
    Scenario(
        name="orca_qmmm_region_correction",
        program="orca",
        expected_kind="qmmm",
        workflow="multi_turn_qmmm_region_correction",
        turns=(
            "Prepare ORCA QM/MM for enzyme_site_mt.xyz, charge -1 multiplicity 1, high-level atoms 1-8.",
            "Correction: high-level atoms should be 1-10 and 27-31. Keep the total charge and multiplicity.",
        ),
    ),
    Scenario(
        name="database_selector_ambiguous_then_record_id",
        program="gaussian",
        expected_kind="opt",
        workflow="database_selector_clarification",
        turns=(
            "Use my ChemSmart database file co2_library.db and optimize the CO2 adduct with Gaussian as neutral singlet.",
            "Use record_id co2-adduct-17 from co2_library.db. Now prepare the Gaussian optimization command.",
        ),
    ),
    Scenario(
        name="fake_runner_local_smoke_request",
        program="orca",
        expected_kind="sp",
        workflow="fake_runner_command_synthesis",
        turns=(
            "Make a local fake-runner ORCA single point for smoke_ligand_mt.xyz, charge 0 multiplicity 1. It should be safe to test without ORCA installed.",
            "Keep it as fake execution and do not submit to a server.",
        ),
    ),
    Scenario(
        name="scheduler_submit_then_fake_correction",
        program="gaussian",
        expected_kind="opt",
        workflow="execution_mode_correction_submit_to_fake",
        turns=(
            "Submit a Gaussian optimization of hpc_probe_mt.xyz, neutral singlet, to the slurm_gpu server.",
            "Correction: do not submit. I only want a local fake-run command that exercises the runner path.",
        ),
    ),
    Scenario(
        name="project_yaml_heavy_atom_then_tddft",
        program="gaussian",
        expected_kind="td",
        workflow="project_yaml_mixed_basis_to_tddft",
        turns=(
            "Create a Gaussian project YAML named bromotd.yaml for CAM-B3LYP-D3BJ with def2-SVP on light atoms and def2-SVPD for Br, acetonitrile SMD, for TD-DFT studies.",
            "Now prepare a TD-DFT command for bromo_dye_mt.xyz, charge 0 multiplicity 1, 18 singlet states, root 4, using that project.",
        ),
    ),
    Scenario(
        name="project_yaml_solvent_switch_then_sp",
        program="orca",
        expected_kind="sp",
        workflow="project_yaml_solvent_switch_to_sp",
        turns=(
            "Build an ORCA project YAML named redoxsp.yaml for PBE0/def2-TZVP single points with CPCM water.",
            "Correction: use acetonitrile instead of water, then make the ORCA single-point command for redox_pair_mt.xyz, charge -1 multiplicity 2.",
        ),
    ),
    Scenario(
        name="orca_dlpno_aux_basis_repair",
        program="orca",
        expected_kind="sp",
        workflow="stale_aux_basis_flag_repair",
        turns=(
            "Repair this ORCA DLPNO single-point command. Intent: dlpno_probe_mt.xyz, neutral singlet, def2/J auxiliary basis. Bad command: chemsmart run orca -p demo -f dlpno_probe_mt.xyz -c 0 -m 1 -a def2/J sp",
            "Use the current ChemSmart ORCA aux-basis flag and keep the job as a single point.",
        ),
    ),
    Scenario(
        name="gaussian_ts_noeigentest_route_repair",
        program="gaussian",
        expected_kind="ts",
        workflow="route_repair_missing_noeigentest",
        turns=(
            "Repair this Gaussian TS command. Intent: aldehyde_shift_mt.xyz, neutral singlet, TS optimization with calcfc and noeigentest. Bad command: chemsmart run gaussian -p demo -f aldehyde_shift_mt.xyz -c 0 -m 1 -o calcfc ts",
            "The generated TS route must include noeigentest and must not leak a Python list into the route.",
        ),
    ),
    Scenario(
        name="decline_method_only_then_structure",
        program="gaussian",
        expected_kind="opt",
        workflow="method_text_then_structure_resume",
        turns=(
            "The paper says B3LYP-D3BJ/def2-SVP was used for minima; run that calculation for my catalyst, but I have not provided coordinates.",
            "Use catalyst_minimum_mt.xyz, charge 0 multiplicity 1. Prepare the Gaussian optimization command only.",
        ),
    ),
    Scenario(
        name="korean_modred_not_scan_correction",
        program="gaussian",
        expected_kind="modred",
        workflow="korean_kind_correction_scan_to_modred",
        turns=(
            "Gaussian에서 kr_bridge_modred_mt.xyz의 결합 4-8을 고정하면서 최적화해줘. 처음에는 scan처럼 생각했는데, 실제로는 값을 변화시키지 않을 거야.",
            "그러니까 relaxed scan이 아니라 modredundant constraint 명령으로 만들어줘. 전하 0, 다중도 1.",
        ),
    ),
    Scenario(
        name="chinese_tddft_missing_nstates_resume",
        program="gaussian",
        expected_kind="td",
        workflow="chinese_ask_resume_tddft_slots",
        turns=(
            "请为 cn_dye_td_mt.xyz 准备 Gaussian TD-DFT 计算，电荷 0，多重度 1，但我还没有给激发态数目。",
            "使用 16 个 singlet states，root 2，并且启用 eqsolv。",
        ),
    ),
    Scenario(
        name="resp_not_two_command_chain",
        program="gaussian",
        expected_kind="resp",
        workflow="separate_workflow_not_shell_chain",
        turns=(
            "Optimize ligand_resp_chain_mt.xyz and then make RESP charges in one command.",
            "Correction: do not chain shell commands. Prepare only the RESP charge fitting command for the already optimized ligand_resp_chain_mt.xyz, charge -1 multiplicity 1.",
        ),
    ),
    Scenario(
        name="gaussian_qmmm_partition_resume",
        program="gaussian",
        expected_kind="qmmm",
        workflow="literature_style_active_region_partition",
        source_reference=(
            "QM/MM active-region refinement: structural atom partitions are "
            "user inputs while method layers remain project-YAML owned."
        ),
        turns=(
            "For an enzyme active-site model, prepare only the Gaussian QM/MM "
            "optimization command for co2_enzyme_qmmm_mt.xyz. The QM region is "
            "atoms 1-8, total charge 0, multiplicity 1; keep the loaded project "
            "settings and do not add method or basis flags.",
            "The MM low layer is atoms 9-16. Keep the same structure and make "
            "the grounded ChemSmart command now.",
        ),
    ),
    Scenario(
        name="orca_qmmm_explicit_qm_region",
        program="orca",
        expected_kind="qmmm",
        workflow="literature_style_multiscale_active_site_refinement",
        source_reference=(
            "QM/MM enzyme refinement: ORCA uses the explicit high-level atom "
            "region and treats the remainder as the MM environment."
        ),
        turns=(
            "A catalytic intermediate from an enzyme mechanism study needs an "
            "ORCA QM/MM optimization. Use enzyme_co2_orca_qmmm_mt.xyz, atoms "
            "1-12 as the high-level region, total charge -1, multiplicity 1, "
            "and the loaded project. Generate only the ChemSmart command.",
        ),
    ),
    Scenario(
        name="gaussian_qrc_from_ts_frequency_output",
        program="gaussian",
        expected_kind="qrc",
        workflow="normal_mode_qrc_from_verified_transition_state",
        source_reference=(
            "Reaction-path workflow: a frequency-characterized transition-state "
            "output supplies a normal mode for displaced QRC optimizations."
        ),
        turns=(
            "I have a Gaussian transition-state frequency output named "
            "co2_insertion_ts_qrc_mt.log. Starting from normal mode 2 with "
            "amplitude 1.2, set up the Gaussian QRC optimization workflow for "
            "a neutral singlet using the loaded project. Return only the "
            "ChemSmart command.",
        ),
    ),
    Scenario(
        name="gaussian_trajectory_tail_refinement",
        program="gaussian",
        expected_kind="traj",
        workflow="trajectory_frame_selection_then_dft_refinement",
        source_reference=(
            "Multistage computational workflow: select representative late "
            "trajectory frames before higher-level Gaussian refinement."
        ),
        turns=(
            "From the molecular-dynamics trajectory md_tail_traj_mt.xyz, use "
            "the loaded Gaussian project to prepare a trajectory workflow that "
            "refines five structures from the end of the trajectory. The system "
            "is a neutral singlet. Give only the grounded ChemSmart command.",
        ),
    ),
]


# Deep (>=4 turn) scenarios that specifically feed the rare kinds measured
# under ten trusted rows on 2026-07-12: orca.modred, orca.qrc, orca.irc,
# gaussian.crest, gaussian.dias. Every scenario is hand-authored, ends in the
# target kind, and is phrased distinctly (chitchat -> partial -> correction ->
# constraint -> command) so the rows do not collapse onto one query skeleton.
RARE4TURN_SCENARIOS: list[Scenario] = [
    Scenario(
        name="rk_orca_modred_distance_freeze_deep",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_constraint_multiturn",
        turns=(
            "I'm studying a pre-reaction complex and want to hold one forming "
            "bond fixed while everything else relaxes in ORCA. Can this be done "
            "with the loaded project?",
            "The structure is rk_precomplex_mt.xyz. I first thought a relaxed "
            "surface scan, but I do NOT want to vary the distance at all.",
            "Right, it's a single frozen internal coordinate: the distance "
            "between atoms 3 and 11 should stay constrained during optimization "
            "— a modredundant constraint, not a Cartesian atom freeze.",
            "Charge is 0 and multiplicity 1. Give me only the ChemSmart command "
            "as a modred constraint job, no project YAML edits.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_angle_correction_deep",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_angle_multiturn",
        turns=(
            "Set up a constrained ORCA optimization for rk_bridge_angle_mt.xyz "
            "using the current project.",
            "The constraint is an angle, not a bond. I keep mixing them up.",
            "Freeze the angle defined by atoms 4, 7 and 9; do not scan it.",
            "Neutral singlet. Command only, please, and keep it as a modred "
            "constraint job rather than a plain opt.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_from_saddle_deep",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_multiturn",
        turns=(
            "A reviewer asked us to confirm that our ORCA saddle point really "
            "connects the intended minima. What's the cleanest way?",
            "We have the transition-state output rk_saddle_orca_mt.out. I'd like "
            "to displace along the imaginary mode and reoptimize.",
            "Yes, both forward and reverse displacements from that saddle.",
            "The system is an anion: charge -1, multiplicity 1. Return only the "
            "ChemSmart command using the loaded ORCA project.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_neutral_deep",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_neutral_multiturn",
        turns=(
            "I need to validate a transition state from an ORCA frequency job.",
            "The output file is rk_ts_freq_orca_mt.out and the imaginary mode is "
            "the proton transfer we care about.",
            "Please set up the quasi-reaction-coordinate follow from that mode.",
            "It's neutral, multiplicity 1, current project. Command only.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_direction_resume_deep",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_multiturn",
        turns=(
            "I want to trace the reaction path away from an ORCA transition "
            "state, but I haven't decided the direction yet.",
            "The saddle structure is rk_orca_ts_mt.xyz and the project is "
            "already loaded.",
            "Let's do both directions from the saddle.",
            "Charge 0, multiplicity 1. Give the ChemSmart IRC command only, and "
            "keep it as an IRC, not a TS reoptimization.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_chinese_deep",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_chinese_multiturn",
        turns=(
            "我有一个 ORCA 过渡态，想沿反应坐标走一遍确认连接的极小点。",
            "过渡态文件是 rk_irc_cn_mt.xyz，项目已经加载好了，不要改 YAML。",
            "两个方向都要走，forward 和 reverse。",
            "电荷 0，多重度 1。只给我 ChemSmart 命令行。",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_then_opt_deep",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_multiturn",
        turns=(
            "I have a very floppy molecule and want to sample its conformers "
            "before any DFT work, using the loaded Gaussian project.",
            "The structure file is rk_flexible_chain_mt.xyz.",
            "Use a CREST conformer search, then take the lowest-energy "
            "conformers forward to a Gaussian optimization.",
            "Neutral singlet. Return only the ChemSmart command, no YAML.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_korean_deep",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_korean_multiturn",
        turns=(
            "유연한 리간드의 형태 이성질체를 먼저 탐색하고 싶어. 현재 Gaussian "
            "프로젝트를 그대로 쓰면 돼.",
            "구조 파일은 rk_ligand_flex_mt.xyz 이고, YAML은 새로 만들지 마.",
            "CREST 형태 탐색을 돌린 다음, 낮은 에너지 구조들을 opt로 이어줘.",
            "전하 0, 다중도 1. ChemSmart 명령어만 알려줘.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_from_irc_deep",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_irc_multiturn",
        turns=(
            "A reviewer wants an activation-strain (distortion/interaction) "
            "analysis along our cycloaddition path. Is that supported?",
            "I have the IRC output rk_cycloaddition_irc_mt.log from Gaussian.",
            "Fragment one is atoms 1 through 6; sample every second point along "
            "the IRC.",
            "Neutral singlet overall. Give only the ChemSmart DIAS command using "
            "the loaded project.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_ts_mode_correction_deep",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_ts_multiturn",
        turns=(
            "Set up a distortion/interaction analysis for rk_da_adduct_ts_mt.log "
            "with the current Gaussian project.",
            "I initially said IRC mode, but this file is a single transition "
            "state, so use TS mode instead.",
            "Fragment one is atoms 1 to 8.",
            "Charge 0, multiplicity 1. Command only, and keep it a DIAS job, not "
            "a plain optimization.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_dihedral_deep",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_dihedral_multiturn",
        turns=(
            "I need to lock a torsion in place while optimizing a rotamer with "
            "ORCA.",
            "The molecule is rk_rotamer_mt.xyz and the project is loaded.",
            "Constrain the dihedral across atoms 5, 6, 7 and 8; keep its current "
            "value fixed, don't scan through it.",
            "Neutral singlet. Only the ChemSmart command, as a modred constraint "
            "job.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_charged_fragments_deep",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_charged_multiturn",
        turns=(
            "We want a distortion/interaction decomposition on an anionic "
            "Diels-Alder path computed in Gaussian.",
            "The IRC file is rk_anionic_irc_mt.log and the loaded project has "
            "our method already.",
            "Fragment one is atoms 1 to 4; take every third point.",
            "Overall charge is -1, multiplicity 1. Return only the ChemSmart "
            "DIAS command.",
        ),
    ),
    # --- batch 2 (positions 13-22): distinct chemistry + phrasing, weighted to
    # the shortest kinds (orca.modred, orca.qrc, gaussian.crest). ---
    Scenario(
        name="rk_orca_modred_hbond_freeze_b2",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_hbond_b2",
        turns=(
            "I'm probing a hydrogen-bonded dimer and need to hold the H-bond "
            "length fixed while the rest relaxes in ORCA.",
            "The dimer geometry is rk_hbond_dimer_mt.xyz, project already set.",
            "Keep the donor-acceptor distance between atoms 6 and 14 fixed; "
            "it's a constraint, not a scan.",
            "Neutral singlet. Only the ChemSmart command, as a modredundant "
            "constraint optimization.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_two_bonds_b2",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_two_bonds_b2",
        turns=(
            "For a concerted TS mimic I want two bonds constrained at once "
            "during an ORCA optimization.",
            "Structure rk_concerted_mt.xyz, loaded project, no YAML changes.",
            "Freeze both the 2-7 distance and the 4-9 distance simultaneously; "
            "hold them, do not scan.",
            "Charge 0, multiplicity 1. Command only, keep it a modred job.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_korean_b2",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_korean_b2",
        turns=(
            "ORCA로 촉매 중간체를 최적화하는데, 특정 결합 하나만 고정하고 싶어.",
            "구조 파일은 rk_catalyst_int_mt.xyz 이고 프로젝트는 이미 로드돼 있어.",
            "원자 3과 8 사이 거리를 고정해줘. 스캔이 아니라 그냥 값을 유지하는 "
            "제약이야.",
            "전하 0, 다중도 1. modred 제약 명령어만 줘.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_dlpno_context_b2",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_dlpno_b2",
        turns=(
            "We located an ORCA saddle for a radical rearrangement and a "
            "reviewer doubts it connects the right species.",
            "The TS output is rk_radical_ts_mt.out, project loaded.",
            "Do the quasi-reaction-coordinate displacement along the imaginary "
            "mode, both directions.",
            "It's a doublet radical: charge 0, multiplicity 2. Command only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_forward_only_b2",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_forward_b2",
        turns=(
            "I need to check which minimum one side of an ORCA transition "
            "state falls to.",
            "The frequency output is rk_halfpath_ts_mt.out and the project is "
            "already configured.",
            "Just follow the imaginary mode; a QRC-style displacement is what I "
            "want.",
            "Cation, charge +1, multiplicity 1. Return only the ChemSmart "
            "command using the loaded project.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_chinese_b2",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_chinese_b2",
        turns=(
            "我们用 ORCA 找到了一个过渡态，想确认它连接的反应物和产物。",
            "过渡态频率输出是 rk_saddle_check_mt.out，项目已加载，别改 YAML。",
            "沿虚频模式做 quasi-reaction coordinate 位移，正反两个方向都要。",
            "电荷 0，多重度 1。只要 ChemSmart 命令行。",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_solvent_b2",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_solvent_b2",
        turns=(
            "Before any DFT, I want a conformational search on a flexible drug "
            "candidate with Gaussian follow-up.",
            "The molecule is rk_druglike_flex_mt.xyz and the project is loaded.",
            "Run a CREST search and carry the lowest conformers into a Gaussian "
            "geometry optimization.",
            "Charge 0, multiplicity 1. Only the ChemSmart command, no YAML "
            "edits.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_anion_b2",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_anion_b2",
        turns=(
            "This carboxylate has many rotatable bonds and I need its "
            "conformers sampled first.",
            "Structure rk_carboxylate_mt.xyz, current Gaussian project, keep "
            "the YAML as is.",
            "Use CREST conformer sampling followed by Gaussian optimization of "
            "the best structures.",
            "It's an anion: charge -1, multiplicity 1. Command only.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_reverse_only_b2",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_reverse_b2",
        turns=(
            "I want to walk down only one side of an ORCA transition state to "
            "confirm the product.",
            "The saddle structure is rk_prod_side_ts_mt.xyz, project loaded.",
            "Follow the reaction path in the reverse direction only, keep it an "
            "IRC.",
            "Charge 0, multiplicity 1. ChemSmart command only, do not turn it "
            "into a TS reoptimization.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_every_point_b2",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_dense_b2",
        turns=(
            "I need a fine-grained distortion/interaction scan along a "
            "cycloaddition IRC in Gaussian.",
            "The IRC output is rk_fine_irc_mt.log, project already has the "
            "method.",
            "Fragment one is atoms 1 to 5; sample every single point this time.",
            "Neutral singlet. Only the ChemSmart DIAS command.",
        ),
    ),
    # --- batch 3 (positions 23-36): the final turn names the loaded 'demo'
    # project explicitly, which the teacher was intermittently dropping; still
    # distinct chemistry/phrasing per item, weighted to crest/modred/qrc. ---
    Scenario(
        name="rk_gaussian_crest_macrocycle_b3",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_macrocycle_b3",
        turns=(
            "A macrocyclic host has too many conformers to guess by hand; I "
            "want them sampled before Gaussian work.",
            "The structure is rk_macrocycle_mt.xyz.",
            "Run a CREST conformer search, then optimize the lowest conformers "
            "with Gaussian.",
            "Neutral singlet. Use the loaded 'demo' project and return only the "
            "ChemSmart command.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_phosphine_b3",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_phosphine_b3",
        turns=(
            "This bulky phosphine ligand has several accessible rotamers I need "
            "to enumerate first.",
            "File is rk_phosphine_mt.xyz; do not create new YAML.",
            "CREST conformer search first, then a Gaussian optimization of the "
            "best ones.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_zwitterion_b3",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_zwitterion_b3",
        turns=(
            "An amino-acid zwitterion is very floppy and I want a conformer "
            "search before anything else.",
            "The geometry is rk_zwitterion_mt.xyz.",
            "Do a CREST search and hand the low-energy conformers to a Gaussian "
            "optimization.",
            "Overall neutral, multiplicity 1. Use the 'demo' project; ChemSmart "
            "command only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_metal_ligand_b3",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_metal_b3",
        turns=(
            "I'm optimizing a metal complex but need to hold the metal-ligand "
            "bond at its crystal value with ORCA.",
            "Structure rk_metal_complex_mt.xyz.",
            "Constrain the distance between atoms 1 and 12; it must stay fixed, "
            "this is not a scan.",
            "Charge +2, multiplicity 1. Use the loaded 'demo' project; give the "
            "modred constraint command only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_pucker_dihedral_b3",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_pucker_b3",
        turns=(
            "I want to lock a ring-puckering dihedral while relaxing the rest of "
            "a sugar with ORCA.",
            "The molecule is rk_sugar_ring_mt.xyz.",
            "Freeze the dihedral defined by atoms 2, 3, 4 and 5; keep its value, "
            "do not scan it.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart modred "
            "command.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_chinese_b3",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_chinese_b3",
        turns=(
            "用 ORCA 优化一个反应中间体，但要保持一个关键键长不变。",
            "结构文件是 rk_intermediate_cn_mt.xyz，项目已经加载。",
            "固定原子 5 和 9 之间的距离，这是约束不是扫描。",
            "电荷 0，多重度 1，用已加载的 'demo' 项目。只给 modred 约束命令。",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_migration_b3",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_migration_b3",
        turns=(
            "We have an ORCA saddle for a 1,2-migration and want to confirm it "
            "links reactant and product.",
            "The frequency output is rk_migration_ts_mt.out.",
            "Displace along the imaginary mode both ways, a quasi-reaction "
            "coordinate follow.",
            "Neutral singlet, loaded 'demo' project. ChemSmart command only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_ringopen_b3",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_ringopen_b3",
        turns=(
            "A reviewer wants proof our ORCA ring-opening saddle connects the "
            "closed and open forms.",
            "The TS output is rk_ringopen_ts_mt.out; keep the YAML as is.",
            "Follow the imaginary mode in both directions from the saddle.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Return only the "
            "ChemSmart command.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_korean_b3",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_korean_b3",
        turns=(
            "ORCA로 찾은 전이상태가 원하는 반응물과 생성물을 잇는지 확인하고 싶어.",
            "진동수 계산 출력은 rk_ts_korean_mt.out 이고 프로젝트는 로드돼 있어.",
            "허수 진동 모드를 따라 양방향으로 변위시켜줘. QRC 방식으로.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. ChemSmart 명령어만 줘.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_solvent_b3",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_solvent_b3",
        turns=(
            "I need to trace both sides of an ORCA transition state for a "
            "solution-phase proton transfer.",
            "The saddle geometry is rk_proton_ts_mt.xyz.",
            "Run the IRC in both directions from the saddle.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Keep it an IRC; "
            "command only.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_forward_b3",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_forward_b3",
        turns=(
            "I only need the product side of an ORCA saddle traced out.",
            "The transition-state file is rk_fwd_ts_mt.xyz.",
            "Follow the reaction path forward only, as an IRC.",
            "Charge -1, multiplicity 1, loaded 'demo' project. ChemSmart command "
            "only, not a TS reoptimization.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_asymmetric_b3",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_asym_b3",
        turns=(
            "I want an activation-strain analysis on an asymmetric "
            "cycloaddition path from Gaussian.",
            "The IRC output is rk_asym_irc_mt.log.",
            "Fragment one is atoms 1 through 7; sample every second point.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart DIAS "
            "command.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_ts_mode_b3",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_ts_b3",
        turns=(
            "For a single transition state I need a distortion/interaction "
            "decomposition in Gaussian.",
            "The TS output is rk_single_ts_mt.log.",
            "Use TS mode; fragment one is atoms 1 to 9.",
            "Charge 0, multiplicity 1, loaded 'demo' project. DIAS command only, "
            "not a plain optimization.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_angle_b3",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_angle_b3",
        turns=(
            "I need to keep a bite angle fixed while optimizing a chelate with "
            "ORCA.",
            "The structure is rk_chelate_mt.xyz.",
            "Constrain the angle across atoms 3, 1 and 8; hold it fixed, this is "
            "not a scan.",
            "Charge 0, multiplicity 1, loaded 'demo' project. modred constraint "
            "command only.",
        ),
    ),
    # --- batch 4 (positions 37-50): fresh chemistry, weighted to the laggards
    # crest / qrc / modred. Final turn names the 'demo' project. ---
    Scenario(
        name="rk_gaussian_crest_peptide_b4",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_peptide_b4",
        turns=(
            "A short peptide has a huge conformational space and I want it "
            "sampled before any Gaussian optimization.",
            "The structure file is rk_peptide_chain_mt.xyz.",
            "Run a CREST conformer search first, then optimize the lowest few "
            "with Gaussian.",
            "Neutral singlet, loaded 'demo' project. Return only the ChemSmart "
            "command.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_diamine_b4",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_diamine_b4",
        turns=(
            "This flexible diamine linker needs a conformer search before DFT.",
            "The geometry is rk_diamine_linker_mt.xyz; keep the YAML unchanged.",
            "Use CREST to sample conformers, then carry the best into a Gaussian "
            "optimization.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_cation_b4",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_cation_b4",
        turns=(
            "An ammonium cation with long alkyl chains has many rotamers to "
            "enumerate first.",
            "The file is rk_ammonium_cat_mt.xyz.",
            "Do a CREST conformer search followed by Gaussian optimization of "
            "the low-energy structures.",
            "It's a cation: charge +1, multiplicity 1, loaded 'demo' project. "
            "ChemSmart command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_korean_b4",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_korean_b4",
        turns=(
            "유연한 크라운 에터의 형태를 먼저 탐색한 뒤 DFT를 하고 싶어.",
            "구조는 rk_crown_ether_mt.xyz 이고 YAML은 그대로 둬.",
            "CREST 형태 탐색 후 낮은 에너지 구조들을 Gaussian opt로 이어줘.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. ChemSmart 명령어만.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_hydride_b4",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_hydride_b4",
        turns=(
            "We found an ORCA saddle for a hydride transfer and need to confirm "
            "the connected minima.",
            "The frequency output is rk_hydride_ts_mt.out.",
            "Displace along the imaginary mode both directions, QRC style.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart command.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_epoxide_b4",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_epoxide_b4",
        turns=(
            "A reviewer wants proof that our ORCA epoxidation saddle links "
            "alkene and epoxide.",
            "The TS output is rk_epox_ts_mt.out; don't touch the YAML.",
            "Follow the imaginary mode in both directions from the saddle.",
            "Charge 0, multiplicity 1, loaded 'demo' project. ChemSmart command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_anion_b4",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_anion_b4",
        turns=(
            "I need to validate an ORCA transition state for an anionic "
            "elimination.",
            "The frequency file is rk_elim_ts_mt.out and the project is loaded.",
            "Set up the quasi-reaction-coordinate follow along the imaginary "
            "mode, both ways.",
            "Anion: charge -1, multiplicity 1, loaded 'demo' project. Command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_chinese_b4",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_chinese_b4",
        turns=(
            "我们用 ORCA 找到一个质子转移过渡态，想确认它连接的两个极小点。",
            "频率输出文件是 rk_pt_ts_cn_mt.out，项目已加载，不要改 YAML。",
            "沿虚频模式做 quasi-reaction coordinate 位移，正反两方向。",
            "电荷 0，多重度 1，用已加载的 'demo' 项目。只要 ChemSmart 命令。",
        ),
    ),
    Scenario(
        name="rk_orca_modred_cage_b4",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_cage_b4",
        turns=(
            "I'm optimizing a cage compound with ORCA but must keep one "
            "cross-cage distance fixed.",
            "The structure is rk_cage_compound_mt.xyz.",
            "Constrain the distance between atoms 7 and 19; hold it, don't scan.",
            "Neutral singlet, loaded 'demo' project. modred constraint command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_backbone_b4",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_backbone_b4",
        turns=(
            "I want to hold a backbone dihedral while relaxing a foldamer with "
            "ORCA.",
            "The geometry is rk_foldamer_mt.xyz.",
            "Freeze the dihedral across atoms 6, 7, 8 and 9; keep its value.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Only the ChemSmart "
            "modred command.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_korean_b4",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_korean_b4",
        turns=(
            "ORCA로 이량체를 최적화하는데 두 단량체 사이 거리를 고정하고 싶어.",
            "구조 파일은 rk_dimer_kr_mt.xyz, 프로젝트는 이미 로드됨.",
            "원자 4와 15 사이 거리를 고정해줘. 스캔 아니고 제약이야.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. modred 명령어만.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_carbene_b4",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_carbene_b4",
        turns=(
            "I need to trace both sides of an ORCA saddle for a carbene "
            "insertion.",
            "The transition-state geometry is rk_carbene_ts_mt.xyz.",
            "Run the IRC in both directions from the saddle.",
            "Neutral singlet, loaded 'demo' project. Keep it an IRC; command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_hda_b4",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_hda_b4",
        turns=(
            "I want an activation-strain analysis along a hetero-Diels-Alder "
            "path from Gaussian.",
            "The IRC output is rk_hda_irc_mt.log.",
            "Fragment one is atoms 1 through 6; sample every second point.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart DIAS "
            "command.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_13dipolar_b4",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_dipolar_b4",
        turns=(
            "For a 1,3-dipolar cycloaddition I need a distortion/interaction "
            "decomposition in Gaussian.",
            "The IRC file is rk_dipolar_irc_mt.log.",
            "Fragment one is atoms 1 to 5; take every second point.",
            "Charge 0, multiplicity 1, loaded 'demo' project. DIAS command only, "
            "not a plain optimization.",
        ),
    ),
    # --- batch 5 (positions 51-65): qrc turns now say "QRC, not IRC" to fix the
    # qrc<->irc confusion that produced `irc -d both` in batch 4. ---
    Scenario(
        name="rk_orca_qrc_notirc_sn2_b5",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_notirc_b5",
        turns=(
            "I have an ORCA saddle for an SN2 and want a quick sanity check on "
            "which minima it falls to.",
            "The frequency output is rk_sn2_check_ts_mt.out.",
            "Run a QRC: displace along the single imaginary normal mode and "
            "reoptimize. This is a QRC job, not an IRC path integration.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart command.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_notirc_tautomer_b5",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_tautomer_b5",
        turns=(
            "We optimized an ORCA transition state for a tautomerization and "
            "need to know the connected tautomers.",
            "The TS frequency file is rk_tautomer_ts_mt.out; keep the YAML.",
            "Use a quasi-reaction-coordinate displacement from the imaginary "
            "mode, then reoptimize. Keep it a QRC job, do not switch to IRC.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_notirc_cation_b5",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_cation_b5",
        turns=(
            "An ORCA saddle for a cationic rearrangement needs a fast "
            "minima-connectivity check.",
            "The output is rk_cation_rearr_ts_mt.out.",
            "Do a QRC displacement along the imaginary mode and reoptimize — a "
            "QRC job specifically, not an IRC.",
            "Cation: charge +1, multiplicity 1, loaded 'demo' project. Command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_notirc_korean_b5",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_korean_b5",
        turns=(
            "ORCA 전이상태가 어떤 극소점으로 떨어지는지 빠르게 확인하고 싶어.",
            "진동수 출력은 rk_qrc_kr_ts_mt.out 이고 프로젝트는 로드됨.",
            "허수 모드를 따라 변위시키고 재최적화하는 QRC 작업이야. IRC 경로 "
            "적분이 아니라 QRC로 해줘.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. ChemSmart 명령어만.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_notirc_radical_b5",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_radical_b5",
        turns=(
            "A radical ORCA saddle needs its connected minima confirmed quickly.",
            "The frequency output is rk_radical2_ts_mt.out.",
            "Run a QRC: displace along the imaginary normal mode and reoptimize. "
            "This must be a QRC job, not an IRC.",
            "Doublet: charge 0, multiplicity 2, loaded 'demo' project. Command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_glycol_b5",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_glycol_b5",
        turns=(
            "A polyethylene-glycol fragment is very floppy; I want its "
            "conformers sampled before Gaussian work.",
            "The structure is rk_peg_fragment_mt.xyz.",
            "Run a CREST conformer search, then optimize the lowest ones with "
            "Gaussian.",
            "Neutral singlet, loaded 'demo' project. ChemSmart command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_thiol_b5",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_thiol_b5",
        turns=(
            "This flexible dithiol has many rotamers I need enumerated first.",
            "The file is rk_dithiol_mt.xyz; don't change the YAML.",
            "Use CREST to sample conformers, then a Gaussian optimization of the "
            "best structures.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_amide_b5",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_amide_b5",
        turns=(
            "A flexible secondary amide needs conformer sampling before DFT "
            "refinement.",
            "The geometry is rk_amide_flex_mt.xyz.",
            "Do a CREST search and carry the low-energy conformers into a "
            "Gaussian optimization.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart command.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_chinese_b5",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_chinese_b5",
        turns=(
            "一个柔性的糖分子构象很多，我想先做构象搜索再用 Gaussian 优化。",
            "结构文件是 rk_sugar_cn_mt.xyz，不要改 YAML。",
            "先用 CREST 搜索构象，再把能量最低的几个交给 Gaussian 优化。",
            "电荷 0，多重度 1，用已加载的 'demo' 项目。只要 ChemSmart 命令。",
        ),
    ),
    Scenario(
        name="rk_orca_modred_spiro_b5",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_spiro_b5",
        turns=(
            "I'm optimizing a spiro compound with ORCA but must hold one "
            "ring-junction distance fixed.",
            "The structure is rk_spiro_mt.xyz.",
            "Constrain the distance between atoms 5 and 13; hold it, this is not "
            "a scan.",
            "Neutral singlet, loaded 'demo' project. modred constraint command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_planar_b5",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_planar_b5",
        turns=(
            "I want to keep a biaryl dihedral fixed while relaxing the rest with "
            "ORCA.",
            "The geometry is rk_biaryl_mt.xyz.",
            "Freeze the dihedral across atoms 3, 4, 5 and 6; keep its current "
            "value, do not scan.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Only the ChemSmart "
            "modred command.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_isomerization_b5",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_isomerization_b5",
        turns=(
            "I need to trace both sides of an ORCA saddle for a cis-trans "
            "isomerization.",
            "The transition-state geometry is rk_isomer_ts_mt.xyz.",
            "Run the IRC in both directions from the saddle.",
            "Neutral singlet, loaded 'demo' project. Keep it an IRC; command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_korean_b5",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_korean_b5",
        turns=(
            "ORCA 전이상태의 양쪽 반응 경로를 모두 추적하고 싶어.",
            "안장점 구조는 rk_irc_kr2_ts_mt.xyz 이고 프로젝트는 로드됨.",
            "안장점에서 양방향으로 IRC를 돌려줘.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. IRC로 유지하고 명령어만.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_simple_irc_b5",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_simple_b5",
        turns=(
            "I want a distortion/interaction analysis along a Diels-Alder IRC "
            "in Gaussian.",
            "The IRC output is rk_da_simple_irc_mt.log.",
            "Fragment one is atoms 1 to 4.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart DIAS "
            "command; use IRC mode.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_ene_b5",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_ene_b5",
        turns=(
            "For an ene reaction I need an activation-strain decomposition in "
            "Gaussian.",
            "The IRC file is rk_ene_irc_mt.log.",
            "Fragment one is atoms 1 to 6.",
            "Charge 0, multiplicity 1, loaded 'demo' project. DIAS command only, "
            "not a plain optimization.",
        ),
    ),
    # --- batch 6 (positions 66-83): crest turns now give an explicit conformer
    # count ("the 5 lowest") so the teacher stops asking how many. ---
    Scenario(
        name="rk_gaussian_crest_urea_b6",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_urea_b6",
        turns=(
            "A flexible bis-urea receptor has many conformers; I want them "
            "sampled before Gaussian DFT.",
            "The structure is rk_bisurea_mt.xyz.",
            "Run a CREST search, then optimize the 5 lowest-energy conformers "
            "with Gaussian.",
            "Neutral singlet, loaded 'demo' project. ChemSmart command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_alkylchain_b6",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_alkyl_b6",
        turns=(
            "A long alkyl thiolate is conformationally floppy and needs sampling "
            "first.",
            "The file is rk_alkyl_thiolate_mt.xyz; keep the YAML.",
            "CREST conformer search, then optimize the 8 lowest conformers with "
            "Gaussian.",
            "Anion: charge -1, multiplicity 1, loaded 'demo' project. Command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_ester_b6",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_ester_b6",
        turns=(
            "This flexible diester has several rotamers I need enumerated "
            "before DFT.",
            "The geometry is rk_diester_mt.xyz.",
            "Do a CREST search and take the 10 lowest-energy conformers into a "
            "Gaussian optimization.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart command.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_amine_b6",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_amine_b6",
        turns=(
            "A flexible triamine needs a conformer search before any "
            "optimization.",
            "The file is rk_triamine_mt.xyz; don't change the YAML.",
            "Use CREST, then optimize the 6 lowest conformers with Gaussian.",
            "Neutral singlet, loaded 'demo' project. ChemSmart command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_diol_b6",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_diol_b6",
        turns=(
            "A flexible 1,4-diol has many H-bonded conformers to sample first.",
            "The structure is rk_14diol_mt.xyz.",
            "CREST conformer search, then optimize the 5 lowest with Gaussian.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_korean_b6",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_korean_b6",
        turns=(
            "유연한 폴리올 분자의 형태를 먼저 탐색하고 싶어.",
            "구조 파일은 rk_polyol_kr_mt.xyz 이고 YAML은 그대로 둬.",
            "CREST 탐색 후 에너지가 가장 낮은 5개 구조를 Gaussian opt로 이어줘.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. ChemSmart 명령어만.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_chinese_b6",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_chinese_b6",
        turns=(
            "一个柔性的多肽片段构象很多，先做构象搜索再优化。",
            "结构文件是 rk_peptide_cn_mt.xyz，别改 YAML。",
            "用 CREST 搜索，然后把能量最低的 8 个构象交给 Gaussian 优化。",
            "电荷 0，多重度 1，用已加载的 'demo' 项目。只要命令。",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_michael_b6",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_michael_b6",
        turns=(
            "I want an activation-strain analysis along a Michael-addition IRC "
            "from Gaussian.",
            "The IRC output is rk_michael_irc_mt.log.",
            "Fragment one is atoms 1 to 5; sample every second point.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart DIAS "
            "command, IRC mode.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_aldol_b6",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_aldol_b6",
        turns=(
            "For an aldol TS I need a distortion/interaction decomposition in "
            "Gaussian.",
            "The transition-state output is rk_aldol_ts_mt.log.",
            "Use TS mode; fragment one is atoms 1 to 8.",
            "Charge 0, multiplicity 1, loaded 'demo' project. DIAS command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_click_b6",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_click_b6",
        turns=(
            "An azide-alkyne click IRC needs an activation-strain analysis in "
            "Gaussian.",
            "The IRC file is rk_click_irc_mt.log.",
            "Fragment one is atoms 1 to 3; sample every third point.",
            "Neutral singlet, loaded 'demo' project. DIAS command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_korean_b6",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_korean_b6",
        turns=(
            "Gaussian IRC를 따라 활성화-변형 분석(distortion/interaction)을 하고 싶어.",
            "IRC 출력은 rk_dias_kr_irc_mt.log 이고 프로젝트는 로드됨.",
            "프래그먼트 1은 원자 1-6, 매 두 번째 점을 사용해.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. DIAS 명령어만.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_out_of_plane_b6",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_oop_b6",
        turns=(
            "I want to hold an out-of-plane bend fixed while optimizing a "
            "porphyrin with ORCA.",
            "The structure is rk_porphyrin_mt.xyz.",
            "Freeze the angle across atoms 2, 9 and 15; keep it, do not scan.",
            "Neutral singlet, loaded 'demo' project. modred constraint command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_disulfide_b6",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_disulfide_b6",
        turns=(
            "Optimizing a peptide with ORCA, I need to hold a disulfide bond "
            "length fixed.",
            "The geometry is rk_disulfide_mt.xyz.",
            "Constrain the distance between atoms 11 and 24; it stays fixed, not "
            "a scan.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Only the ChemSmart "
            "modred command.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_chinese_b6",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_chinese_b6",
        turns=(
            "用 ORCA 优化一个大环，但要固定一个跨环的距离。",
            "结构文件是 rk_macro_cn_mt.xyz，项目已加载。",
            "固定原子 6 和 18 之间的距离，是约束不是扫描。",
            "电荷 0，多重度 1，用已加载的 'demo' 项目。只给 modred 约束命令。",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_grignard_b6",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_grignard_b6",
        turns=(
            "An ORCA saddle for an addition step needs a fast minima check.",
            "The frequency output is rk_addition_ts_mt.out.",
            "Do a QRC: displace along the imaginary mode and reoptimize — a QRC "
            "job, not an IRC.",
            "Neutral singlet, loaded 'demo' project. ChemSmart command only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_decarbox_b6",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_decarbox_b6",
        turns=(
            "We need to confirm the minima connected by an ORCA decarboxylation "
            "saddle.",
            "The TS frequency file is rk_decarb_ts_mt.out; keep the YAML.",
            "Run a quasi-reaction-coordinate displacement and reoptimize, a QRC "
            "job rather than IRC.",
            "Anion: charge -1, multiplicity 1, loaded 'demo' project. Command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_ringclose_b6",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_ringclose_b6",
        turns=(
            "I need to trace both sides of an ORCA saddle for a ring closure.",
            "The transition-state geometry is rk_ringclose_ts_mt.xyz.",
            "Run the IRC in both directions from the saddle.",
            "Neutral singlet, loaded 'demo' project. Keep it an IRC; command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_two_angle_b6",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_twoangle_b6",
        turns=(
            "Optimizing a bidentate complex with ORCA, I must hold two bite "
            "angles fixed at once.",
            "The structure is rk_bidentate_mt.xyz.",
            "Freeze the angle across atoms 2, 1, 3 and the angle across atoms 4, "
            "1, 5; hold both, do not scan.",
            "Charge 0, multiplicity 1, loaded 'demo' project. modred command "
            "only.",
        ),
    ),
    # --- batch 7 (positions 84-105): final push to clear all five kinds. ---
    Scenario(
        name="rk_orca_modred_bond_b7",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_bond_b7",
        turns=(
            "Optimizing an ion pair with ORCA, I need one contact distance held.",
            "The structure is rk_ionpair_mt.xyz.",
            "Constrain the distance between atoms 4 and 10; hold it, not a scan.",
            "Neutral singlet, loaded 'demo' project. modred command only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_torsion_b7",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_torsion_b7",
        turns=(
            "I want a peptide bond torsion locked while relaxing with ORCA.",
            "The geometry is rk_peptidebond_mt.xyz.",
            "Freeze the dihedral across atoms 3, 4, 5 and 6; keep its value.",
            "Charge 0, multiplicity 1, loaded 'demo' project. modred command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_metal_dist_b7",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_metaldist_b7",
        turns=(
            "A metal-carbonyl needs the M-CO distance frozen during ORCA "
            "optimization.",
            "The file is rk_metalco_mt.xyz.",
            "Constrain the distance between atoms 1 and 7; it must stay fixed.",
            "Charge 0, multiplicity 1, loaded 'demo' project. modred command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_korean_b7",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_korean_b7",
        turns=(
            "ORCA로 착물을 최적화하는데 금속-리간드 각도 하나를 고정하고 싶어.",
            "구조는 rk_complex_kr_mt.xyz 이고 프로젝트는 로드됨.",
            "원자 3, 1, 9로 정의되는 각도를 고정해줘. 스캔이 아니야.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. modred 명령어만.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_hbond2_b7",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_hbond2_b7",
        turns=(
            "A catalyst-substrate pair needs its key H-bond distance held in "
            "ORCA.",
            "The structure is rk_cat_sub_mt.xyz.",
            "Constrain the distance between atoms 8 and 21; hold it, not a scan.",
            "Neutral singlet, loaded 'demo' project. modred command only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_ring_b7",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_ring_b7",
        turns=(
            "Optimizing a strained ring with ORCA, I need one ring bond frozen.",
            "The geometry is rk_strained_ring_mt.xyz.",
            "Freeze the distance between atoms 2 and 6; keep it fixed.",
            "Charge 0, multiplicity 1, loaded 'demo' project. modred command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_ester_b7",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_ester_b7",
        turns=(
            "An ORCA saddle for an ester hydrolysis needs a minima check.",
            "The frequency output is rk_hydrolysis_ts_mt.out.",
            "Run a QRC: displace along the imaginary mode and reoptimize — a QRC "
            "job, not an IRC.",
            "Neutral singlet, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_amide_b7",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_amide_b7",
        turns=(
            "We want the minima connected by an ORCA amide-rotation saddle.",
            "The TS frequency file is rk_amiderot_ts_mt.out; keep the YAML.",
            "Displace along the imaginary mode and reoptimize, a QRC job not an "
            "IRC.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_protonation_b7",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_protonation_b7",
        turns=(
            "An ORCA saddle for a protonation step needs its minima confirmed "
            "quickly.",
            "The output is rk_proton_step_ts_mt.out.",
            "Do a QRC displacement along the imaginary mode and reoptimize; keep "
            "it a QRC job, not an IRC.",
            "Cation: charge +1, multiplicity 1, loaded 'demo' project. Command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_qrc_chinese_b7",
        program="orca",
        expected_kind="qrc",
        workflow="rare_orca_qrc_chinese_b7",
        turns=(
            "ORCA 找到一个消除反应过渡态，想确认它连接的极小点。",
            "频率输出是 rk_elim_cn_ts_mt.out，项目已加载，别改 YAML。",
            "沿虚频模式做位移再优化，是 QRC 作业，不是 IRC。",
            "电荷 0，多重度 1，用 'demo' 项目。只要命令。",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_ketone_b7",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_ketone_b7",
        turns=(
            "A flexible diketone needs conformer sampling before Gaussian DFT.",
            "The structure is rk_diketone_mt.xyz.",
            "CREST search, then optimize the 5 lowest conformers with Gaussian.",
            "Neutral singlet, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_sulfon_b7",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_sulfon_b7",
        turns=(
            "A flexible sulfonamide has many rotamers to sample first.",
            "The file is rk_sulfonamide_mt.xyz; keep the YAML.",
            "Use CREST, then optimize the 6 lowest conformers with Gaussian.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_phenol_b7",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_phenol_b7",
        turns=(
            "A polyphenol with rotatable OH groups needs conformer sampling "
            "first.",
            "The geometry is rk_polyphenol_mt.xyz.",
            "Run CREST, then optimize the 8 lowest conformers with Gaussian.",
            "Neutral singlet, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_ether_b7",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_ether_b7",
        turns=(
            "A flexible polyether needs its conformers enumerated before DFT.",
            "The file is rk_polyether_mt.xyz.",
            "CREST conformer search, then optimize the 5 lowest with Gaussian.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_chinese_b7",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_chinese_b7",
        turns=(
            "一个柔性的二胺分子构象很多，先搜索构象再优化。",
            "结构文件是 rk_diamine_cn_mt.xyz，别改 YAML。",
            "用 CREST 搜索，然后把最低的 5 个构象交给 Gaussian 优化。",
            "电荷 0，多重度 1，用 'demo' 项目。只要命令。",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_da_explicit_b7",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_explicit_b7",
        turns=(
            "I want an activation-strain analysis along a Diels-Alder IRC in "
            "Gaussian, IRC mode.",
            "The IRC output is rk_da7_irc_mt.log.",
            "Fragment one is exactly atoms 1-4; the rest is fragment two.",
            "Neutral singlet, loaded 'demo' project. Only the ChemSmart DIAS "
            "command with fragment 1 = atoms 1-4.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_retro_b7",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_retro_b7",
        turns=(
            "For a retro-cycloaddition IRC I need a distortion/interaction "
            "decomposition in Gaussian, IRC mode.",
            "The IRC file is rk_retro_irc_mt.log.",
            "Fragment one is atoms 1-6.",
            "Neutral singlet, loaded 'demo' project. DIAS command only with "
            "fragment 1 = atoms 1-6.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_carbonyl_b7",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_carbonyl_b7",
        turns=(
            "A carbonyl-addition IRC needs an activation-strain analysis in "
            "Gaussian, IRC mode.",
            "The IRC output is rk_carbonyl_irc_mt.log.",
            "Fragment one is atoms 1-3.",
            "Charge 0, multiplicity 1, loaded 'demo' project. DIAS command only "
            "with fragment 1 = atoms 1-3.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_sigmatropic_b7",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_sigmatropic_b7",
        turns=(
            "A sigmatropic-shift IRC needs a distortion/interaction "
            "decomposition in Gaussian, IRC mode.",
            "The IRC file is rk_sigma_irc_mt.log.",
            "Fragment one is atoms 1-5.",
            "Neutral singlet, loaded 'demo' project. DIAS command only with "
            "fragment 1 = atoms 1-5.",
        ),
    ),
    Scenario(
        name="rk_gaussian_dias_korean_b7",
        program="gaussian",
        expected_kind="dias",
        workflow="rare_gaussian_dias_korean_b7",
        turns=(
            "Gaussian IRC를 따라 distortion/interaction 분석을 IRC 모드로 하고 싶어.",
            "IRC 출력은 rk_dias7_kr_irc_mt.log 이고 프로젝트는 로드됨.",
            "프래그먼트 1은 정확히 원자 1-4 야.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. 프래그먼트 1 = 원자 1-4 로 "
            "DIAS 명령어만.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_migration_b7",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_migration_b7",
        turns=(
            "I need to trace both sides of an ORCA saddle for a 1,2-shift.",
            "The transition-state geometry is rk_shift_ts_mt.xyz.",
            "Run the IRC in both directions from the saddle.",
            "Neutral singlet, loaded 'demo' project. Keep it an IRC; command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_irc_addition_b7",
        program="orca",
        expected_kind="irc",
        workflow="rare_orca_irc_addition_b7",
        turns=(
            "Trace both sides of an ORCA saddle for a radical addition.",
            "The saddle geometry is rk_radadd_ts_mt.xyz.",
            "Run the IRC in both directions.",
            "Doublet: charge 0, multiplicity 2, loaded 'demo' project. Keep it "
            "an IRC; command only.",
        ),
    ),
    # --- batch 8 (positions 106-112): mop-up for the last two kinds
    # (modred +3, crest +1). ---
    Scenario(
        name="rk_orca_modred_dist_b8",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_dist_b8",
        turns=(
            "Optimizing a transition-metal amine with ORCA, I need the M-N "
            "distance held.",
            "The structure is rk_mn_amine_mt.xyz.",
            "Constrain the distance between atoms 1 and 5; hold it, not a scan.",
            "Neutral singlet, loaded 'demo' project. modred command only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_dihed_b8",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_dihed_b8",
        turns=(
            "Relaxing a biphenyl with ORCA, I want the inter-ring twist fixed.",
            "The geometry is rk_biphenyl_mt.xyz.",
            "Freeze the dihedral across atoms 5, 6, 7 and 8; keep its value.",
            "Charge 0, multiplicity 1, loaded 'demo' project. modred command "
            "only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_angle_b8",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_angle_b8",
        turns=(
            "A chelating ligand needs its bite angle held during ORCA "
            "optimization.",
            "The file is rk_chelate2_mt.xyz.",
            "Constrain the angle across atoms 4, 1 and 9; hold it, not a scan.",
            "Neutral singlet, loaded 'demo' project. modred command only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_bond2_b8",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_bond2_b8",
        turns=(
            "Optimizing an enolate with ORCA, I need one C-O distance frozen.",
            "The structure is rk_enolate_mt.xyz.",
            "Constrain the distance between atoms 3 and 7; hold it, not a scan.",
            "Anion: charge -1, multiplicity 1, loaded 'demo' project. modred "
            "command only.",
        ),
    ),
    Scenario(
        name="rk_orca_modred_korean_b8",
        program="orca",
        expected_kind="modred",
        workflow="rare_orca_modred_korean_b8",
        turns=(
            "ORCA로 다이머를 최적화하는데 두 조각 사이 거리를 고정하고 싶어.",
            "구조 파일은 rk_dimer2_kr_mt.xyz, 프로젝트는 로드됨.",
            "원자 5와 12 사이 거리를 고정해. 스캔이 아니라 제약이야.",
            "전하 0, 다중도 1, 로드된 'demo' 프로젝트. modred 명령어만.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_final_b8",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_final_b8",
        turns=(
            "A flexible bis-phosphine needs conformer sampling before Gaussian "
            "DFT.",
            "The structure is rk_bisphosphine_mt.xyz.",
            "CREST search, then optimize the 5 lowest conformers with Gaussian.",
            "Neutral singlet, loaded 'demo' project. Command only.",
        ),
    ),
    Scenario(
        name="rk_gaussian_crest_final2_b8",
        program="gaussian",
        expected_kind="crest",
        workflow="rare_gaussian_crest_final2_b8",
        turns=(
            "A flexible glycoside has many conformers to sample first.",
            "The file is rk_glycoside_mt.xyz; keep the YAML.",
            "Use CREST, then optimize the 6 lowest conformers with Gaussian.",
            "Charge 0, multiplicity 1, loaded 'demo' project. Command only.",
        ),
    ),
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--models",
        default=(
            "qwen3-235b-a22b-thinking-2507,"
            "glm-5.2,"
            "deepseek-v4-pro,"
            "qwen3-coder-flash,"
            "kimi-k2.7-code"
        ),
        help="Comma-separated fallback model order.",
    )
    parser.add_argument("--batch-id", default="v17_agentic_workflow")
    parser.add_argument("--start-position", type=int, default=1)
    parser.add_argument("--limit", type=int, default=4)
    parser.add_argument(
        "--scenario-timeout-s",
        type=float,
        default=120.0,
        help="Maximum wall time for one live scenario; 0 disables the guard.",
    )
    parser.add_argument("--task-offset", type=int, default=500)
    parser.add_argument(
        "--groups",
        default="synthesis,project_yaml",
        help="Comma-separated ToolRegistry groups.",
    )
    parser.add_argument(
        "--provider",
        choices=sorted(reasoning_accum.PROVIDERS),
        default="dashscope",
        help="Which api.env credential + endpoint to use for all models.",
    )
    parser.add_argument(
        "--run-label",
        default="",
        help=(
            "Override the run directory name (default: the model id). Use a "
            "unique label to avoid appending into a directory another "
            "collector is actively writing."
        ),
    )
    parser.add_argument(
        "--scenario-set",
        choices=["default", "rare4turn"],
        default="default",
        help="'rare4turn' selects the >=4-turn rare-kind scenarios.",
    )
    args = parser.parse_args(argv)

    os.environ["PATH"] = "/opt/anaconda3/envs/chemsmart/bin:" + os.environ.get(
        "PATH", ""
    )
    sys.path.insert(0, str(REPO))
    os.chdir(REPO)
    logging.disable(logging.INFO)

    key_name, base_url = reasoning_accum.PROVIDERS[args.provider]
    key = dotenv_values(REPO / "api.env").get(key_name, "")
    if not key:
        raise SystemExit(f"api.env is missing {key_name}")

    models = [item.strip() for item in args.models.split(",") if item.strip()]
    if not models:
        raise SystemExit("--models must contain at least one model id")
    if args.run_label and len(models) > 1:
        raise SystemExit("--run-label requires exactly one --models entry")

    def _run_dir_for(model: str) -> Path:
        name = args.run_label or model
        return REPO / "var" / "agent-training" / "runs" / name

    scenario_source = (
        RARE4TURN_SCENARIOS if args.scenario_set == "rare4turn" else SCENARIOS
    )
    selected = [
        (position, scenario)
        for position, scenario in enumerate(scenario_source, start=1)
        if position >= args.start_position
    ][: args.limit]
    if not selected:
        raise SystemExit("no scenarios selected")

    _write_collection_plans(
        models=models,
        batch_id=args.batch_id,
        selected=selected,
        task_offset=args.task_offset,
        run_dir_for=_run_dir_for,
    )

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

    class LongTimeoutProvider(OpenAIProvider):
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
            return self._client.chat.completions.create(**kwargs).model_dump()

    providers: dict[str, LongTimeoutProvider] = {
        model: LongTimeoutProvider(key, model=model, base_url=base_url)
        for model in models
    }

    results: list[dict[str, Any]] = []
    for position, scenario in selected:
        task_index = args.task_offset + position
        task_id = f"agentic-workflow-{task_index:03d}"
        record: dict[str, Any] | None = None
        for model in models:
            run_dir = _run_dir_for(model)
            os.environ["CHEMSMART_AGENT_TRAINING_DIR"] = str(run_dir)
            provider = providers[model]
            providers_mod.get_provider = (
                lambda *a, _provider=provider, **k: _provider
            )
            record = _run_scenario(
                scenario=scenario,
                task_id=task_id,
                task_index=task_index,
                batch_id=args.batch_id,
                model=model,
                provider=provider,
                scenario_timeout_s=args.scenario_timeout_s,
                registry_groups=[
                    item.strip()
                    for item in args.groups.split(",")
                    if item.strip()
                ],
                tools_command=tools_command,
                AgentSession=AgentSession,
                ToolRegistry=ToolRegistry,
                PermissionPolicy=PermissionPolicy,
                PermissionMode=PermissionMode,
                ApprovalDecision=ApprovalDecision,
            )
            _append_result(run_dir, record)
            print(_format_result(record), flush=True)
            results.append(record)
            if record["grade"] != "QUOTA":
                break
            time.sleep(1)
        if record and record["grade"] == "QUOTA":
            print(
                f"{task_id}: all configured models exhausted or unavailable",
                flush=True,
            )

    tally: dict[str, int] = {}
    for row in results:
        tally[row["grade"]] = tally.get(row["grade"], 0) + 1
    print(f"== agentic workflow tally: {tally} ==", flush=True)
    return 0


def _run_scenario(
    *,
    scenario: Scenario,
    task_id: str,
    task_index: int,
    batch_id: str,
    model: str,
    provider: Any,
    scenario_timeout_s: float,
    registry_groups: list[str],
    tools_command: Any,
    AgentSession: Any,
    ToolRegistry: Any,
    PermissionPolicy: Any,
    PermissionMode: Any,
    ApprovalDecision: Any,
) -> dict[str, Any]:
    work = Path(tempfile.mkdtemp(prefix=f"aw-{task_id}-{model[:8]}-"))
    reasoning_accum._prepare_workspace(
        work,
        " ".join(scenario.turns),
        scenario.program,
    )
    original_cwd = Path.cwd()
    os.chdir(work)
    # Attach provenance for every episode AgentSession writes this scenario;
    # the writer reads it from the env var (no core.py change needed).
    from chemsmart.agent.training_log import DATASET_PROVENANCE_ENV

    os.environ[DATASET_PROVENANCE_ENV] = json.dumps(
        _scenario_provenance(scenario, batch_id=batch_id, model=model)
    )
    tools_command.reset_command_tools_state()
    started = time.time()
    turns: list[dict[str, Any]] = []
    record: dict[str, Any] = {
        "batch_id": batch_id,
        "corpus": "agentic-workflow",
        "difficulty": scenario.difficulty,
        "model": model,
        "program": scenario.program,
        "expected_kind": scenario.expected_kind,
        "scenario": scenario.name,
        "source_reference": scenario.source_reference,
        "task_id": task_id,
        "task_index": task_index,
        "turn_count": len(scenario.turns),
        "workflow": scenario.workflow,
    }
    try:
        session = AgentSession(
            provider=provider,
            registry=ToolRegistry.default(groups=registry_groups),
            session_root=work / "sessions" / task_id,
            stage_prompt="unified_agent.md",
        )
        policy = PermissionPolicy(
            mode=PermissionMode.DRIVING, prompt_risky=True
        )
        last_synthesis: dict[str, Any] | None = None
        ask_count = 0
        with _scenario_timeout(scenario_timeout_s):
            for index, turn in enumerate(scenario.turns, start=1):
                out = session.run_loop(
                    turn,
                    policy=policy,
                    approver=lambda req: ApprovalDecision.ALLOW_ONCE,
                )
                synth = reasoning_accum._last_synthesis(out)
                if synth:
                    last_synthesis = synth
                if out.get("ask_user_question"):
                    ask_count += 1
                turns.append(
                    {
                        "turn": index,
                        "ask_user": bool(out.get("ask_user_question")),
                        "assistant_chars": len(
                            out.get("assistant_output") or ""
                        ),
                        "tool_count": len(out.get("tool_outcomes") or []),
                        "synthesis_status": (synth or {}).get("status"),
                        "command": (synth or {}).get("command"),
                    }
                )
        record["latency_s"] = round(time.time() - started, 1)
        record["turns"] = turns
        record["ask_count"] = ask_count
        if last_synthesis:
            command = str(last_synthesis.get("command") or "").strip()
            semantic = last_synthesis.get("semantic") or {}
            verdict = (
                semantic.get("verdict") if isinstance(semantic, dict) else None
            )
            record["command"] = command
            record["gate"] = verdict
            record["reasoning_len"] = len(
                last_synthesis.get("reasoning") or ""
            )
            if (
                scenario.expected_kind != "decline"
                and last_synthesis.get("status") == "ready"
                and verdict in {"ok", "warn"}
                and reasoning_accum._kind_is_present(
                    scenario.expected_kind,
                    command,
                )
            ):
                record["grade"] = "PASS"
            elif last_synthesis.get("status") == "needs_clarification":
                record["grade"] = "ASK"
            else:
                record["grade"] = "WRONG"
        elif scenario.expected_kind == "decline" and ask_count:
            record["grade"] = "PASS_DECLINE"
        elif ask_count:
            record["grade"] = "ASK"
        else:
            record["grade"] = "NOSYNTH"
    except ScenarioTimeout as exc:  # pragma: no cover - live API path
        record["latency_s"] = round(time.time() - started, 1)
        record["turns"] = turns
        record["ask_count"] = ask_count
        record["error"] = str(exc)
        record["grade"] = "TIMEOUT"
    except Exception as exc:  # pragma: no cover - live API path
        message = f"{type(exc).__name__}: {exc}"[:240]
        record["error"] = message
        record["grade"] = "QUOTA" if "free quota" in message.lower() else "ERR"
    finally:
        os.chdir(original_cwd)
        os.environ.pop(DATASET_PROVENANCE_ENV, None)
    return record


def _write_collection_plans(
    *,
    models: list[str],
    batch_id: str,
    selected: list[tuple[int, Scenario]],
    task_offset: int,
    run_dir_for,
) -> None:
    for model in models:
        run_dir = run_dir_for(model)
        plan_dir = run_dir / "collection_plans"
        plan_dir.mkdir(parents=True, exist_ok=True)
        path = (
            plan_dir / f"{batch_id}_agentic-workflow_{int(time.time())}.jsonl"
        )
        with path.open("w", encoding="utf-8") as handle:
            for position, scenario in selected:
                task_index = task_offset + position
                row = {
                    "batch_id": batch_id,
                    "corpus": "agentic-workflow",
                    "difficulty": scenario.difficulty,
                    "expected_kind": scenario.expected_kind,
                    "model": model,
                    "program": scenario.program,
                    "scenario": scenario.name,
                    "source_reference": scenario.source_reference,
                    "task_id": f"agentic-workflow-{task_index:03d}",
                    "task_index": task_index,
                    "turns": list(scenario.turns),
                    "workflow": scenario.workflow,
                    "workspace_yaml_policy": "one project YAML per task workspace",
                }
                handle.write(
                    json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n"
                )
        print(f"collection_plan[{model}]={path}", flush=True)


def _append_result(run_dir: Path, record: dict[str, Any]) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / "agentic_workflow_results.jsonl"
    with path.open("a", encoding="utf-8") as handle:
        handle.write(
            json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n"
        )


def _format_result(record: dict[str, Any]) -> str:
    command = str(record.get("command") or record.get("error") or "")
    return (
        f"{record['model']:32s} {record['task_id']:22s} "
        f"{record['scenario'][:24]:24s} {record['grade']:12s} "
        f"ask={record.get('ask_count', 0)} "
        f"rsn={record.get('reasoning_len', 0):5d} "
        f"{record.get('latency_s', 0):6.1f}s {command[:70]}"
    )


if __name__ == "__main__":
    raise SystemExit(main())
