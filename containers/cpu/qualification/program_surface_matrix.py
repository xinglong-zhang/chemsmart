#!/usr/bin/env python3
"""Provider-free materialization of four licensed-program workflow surfaces.

No Gaussian or ORCA executable is present in the CPU image.  This gate proves
that ChemSmart materializes the requested native inputs. It never claims that
the rendered inputs were executed, and it does not redistribute native output
from either licensed program. The real PySCF qualification separately exercises
the shared typed quasi-harmonic thermochemistry path on a runtime-generated
structured result.
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace

import numpy as np
import yaml

from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    evaluate_quantity_expression,
)
from chemsmart.analysis.result_quantities import ENERGY, make_quantity_value
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianLinkJobSettings,
    GaussianTDDFTJobSettings,
)
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.writer import ORCAInputWriter
from chemsmart.jobs.settings import read_molecular_job_yaml


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _job(settings, molecule, *, label: str):
    return SimpleNamespace(
        settings=settings,
        molecule=molecule,
        label=label,
        jobrunner=SimpleNamespace(num_cores=2, mem_gb=6),
    )


def _gaussian_input(settings, molecule, *, label: str) -> str:
    output = io.StringIO()
    GaussianInputWriter(_job(settings, molecule, label=label))._write_all(
        output
    )
    return output.getvalue()


def _orca_input(settings, molecule, *, label: str) -> str:
    output = io.StringIO()
    ORCAInputWriter(_job(settings, molecule, label=label))._write_all(output)
    return output.getvalue()


def _settings_from_project_yaml(
    *,
    program: str,
    stage: str,
    stage_settings: dict[str, object],
    settings_class,
):
    """Load one stage through the canonical project-YAML routing boundary."""

    document = yaml.safe_dump({stage: stage_settings}, sort_keys=True)
    project_sha256 = hashlib.sha256(document.encode("utf-8")).hexdigest()
    with tempfile.TemporaryDirectory(prefix="chemsmart-surface-") as directory:
        project = Path(directory) / "project.yaml"
        project.write_text(document, encoding="utf-8")
        settings_logger = logging.getLogger("chemsmart.jobs.settings")
        prior_level = settings_logger.level
        try:
            settings_logger.setLevel(logging.ERROR)
            routed = read_molecular_job_yaml(str(project), program=program)
        finally:
            settings_logger.setLevel(prior_level)
    require(stage in routed, f"{program} project did not route stage {stage}")
    return settings_class.from_dict(routed[stage]), project_sha256


def _public_molecules() -> tuple[Molecule, Molecule]:
    nitric_oxide = Molecule(
        symbols=("N", "O"),
        positions=np.asarray(((0.0, 0.0, 0.0), (0.0, 0.0, 1.1508))),
        charge=0,
        multiplicity=2,
    )
    water = Molecule(
        symbols=("O", "H", "H"),
        positions=np.asarray(
            (
                (0.0, 0.0, 0.1174),
                (-0.7570, 0.0, -0.4696),
                (0.7570, 0.0, -0.4696),
            )
        ),
        charge=0,
        multiplicity=1,
    )
    return nitric_oxide, water


def qualify_radical_link(nitric_oxide: Molecule) -> dict[str, object]:
    settings, project_sha256 = _settings_from_project_yaml(
        program="gaussian",
        stage="link",
        settings_class=GaussianLinkJobSettings,
        stage_settings={
            "functional": "UM06-2X",
            "basis": "def2-TZVP",
            "charge": 0,
            "multiplicity": 2,
            "jobtype": "opt",
            "freq": True,
            "stable": "opt",
            "guess": "mix",
            "additional_route_parameters": "NoSymm",
            "title": "Public NO radical stability then opt/freq",
        },
    )
    rendered = _gaussian_input(settings, nitric_oxide, label="no_radical")
    first, linked = (
        section for section in rendered.split("--Link1--") if section.strip()
    )
    first_route = next(
        line for line in first.splitlines() if line.startswith("#")
    )
    linked_route = next(
        line for line in linked.splitlines() if line.startswith("#")
    )
    first_lower, linked_lower = first_route.lower(), linked_route.lower()
    require("um06-2x" in first_lower, "Radical route is not unrestricted")
    require("stable=opt" in first_lower, "Missing radical stability analysis")
    require("guess=mix" in first_lower, "Missing mixed orbital guess")
    require("nosymm" in first_lower, "Missing explicit no-symmetry request")
    require(
        " opt" not in first_lower and "freq" not in first_lower,
        "Stability step incorrectly includes opt/freq",
    )
    for token in (
        "um06-2x",
        "def2tzvp",
        "opt",
        "freq",
        "geom=check",
        "guess=read",
    ):
        require(token in linked_lower, f"Linked radical route omitted {token}")
    require(
        rendered.count("\n0 2\n") == 2, "Doublet identity was not retained"
    )
    return {
        "state": "materialized_not_executed",
        "project_yaml_sha256": project_sha256,
        "first_route": first_route,
        "linked_route": linked_route,
    }


def _typed_focal_cbs() -> dict[str, float]:
    values = {
        "hf_tz": -74.900,
        "hf_qz": -74.940,
        "mp2_corr_tz": -0.200,
        "mp2_corr_qz": -0.225,
        "dlpno_tz": -75.145,
    }
    inputs = tuple(
        make_quantity_value(
            quantity_id=name,
            source_value=value,
            source_unit="hartree",
            value=value,
            unit="hartree",
            dimension=ENERGY,
            evidence_ref=f"artifact:public-{name}#" + "0" * 64,
        )
        for name, value in values.items()
    )
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="public.doublet.ri-mp2-cbs.dlpno-focal",
            inputs=inputs,
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="hf_cbs",
                    operation="scf_exponential_cbs_limit",
                    input_ids=("hf_tz", "hf_qz"),
                    cardinal_numbers=(3, 4),
                    extrapolation_exponent=3.9,
                ),
                QuantityExpressionNodeV1(
                    node_id="mp2_corr_cbs",
                    operation="correlation_inverse_power_cbs_limit",
                    input_ids=("mp2_corr_tz", "mp2_corr_qz"),
                    cardinal_numbers=(3, 4),
                    extrapolation_exponent=3.0,
                ),
                QuantityExpressionNodeV1(
                    node_id="mp2_cbs",
                    operation="add",
                    input_ids=("hf_cbs", "mp2_corr_cbs"),
                ),
                QuantityExpressionNodeV1(
                    node_id="mp2_tz",
                    operation="add",
                    input_ids=("hf_tz", "mp2_corr_tz"),
                ),
                QuantityExpressionNodeV1(
                    node_id="dlpno_focal_correction",
                    operation="subtract",
                    input_ids=("dlpno_tz", "mp2_tz"),
                ),
                QuantityExpressionNodeV1(
                    node_id="focal_total",
                    operation="add",
                    input_ids=("mp2_cbs", "dlpno_focal_correction"),
                ),
            ),
            output_node_ids=(
                "hf_cbs",
                "mp2_corr_cbs",
                "mp2_cbs",
                "dlpno_focal_correction",
                "focal_total",
            ),
        )
    )
    observed = {
        item.quantity_id: float(item.value) for item in receipt.outputs
    }
    require(
        all(math.isfinite(value) for value in observed.values()),
        "Non-finite CBS",
    )
    return observed


def qualify_orca_focal(nitric_oxide: Molecule) -> dict[str, object]:
    settings_and_projects = tuple(
        _settings_from_project_yaml(
            program="orca",
            stage="sp",
            settings_class=ORCAJobSettings,
            stage_settings={
                "ab_initio": "RI-MP2",
                "basis": basis,
                "aux_basis": auxiliary,
                "reference": "rohf",
                "charge": 0,
                "multiplicity": 2,
                "freq": False,
                "title": "Public fixed-geometry NO doublet",
            },
        )
        for basis, auxiliary in (
            ("cc-pVTZ", "cc-pVTZ/C"),
            ("cc-pVQZ", "cc-pVQZ/C"),
        )
    )
    rendered = [
        _orca_input(item, nitric_oxide, label=f"no_rimp2_{index}")
        for index, (item, _project_sha256) in enumerate(
            settings_and_projects, start=3
        )
    ]
    dlpno, dlpno_project_sha256 = _settings_from_project_yaml(
        program="orca",
        stage="sp",
        settings_class=ORCAJobSettings,
        stage_settings={
            "ab_initio": "DLPNO-CCSD(T)",
            "basis": "cc-pVTZ",
            "aux_basis": "cc-pVTZ/C",
            "reference": "rohf",
            "charge": 0,
            "multiplicity": 2,
            "freq": False,
            "title": "Public fixed-geometry NO doublet focal point",
        },
    )
    rendered.append(_orca_input(dlpno, nitric_oxide, label="no_dlpno_tz"))
    for text in rendered:
        lower = text.lower()
        require("* xyz 0 2" in lower, "ORCA focal input lost doublet identity")
        require("hftyp rohf" in lower, "ORCA focal input lost ROHF reference")
        require(
            "freq" not in lower, "Fixed-geometry focal input added frequency"
        )
        require(
            "opt" not in lower, "Fixed-geometry focal input added optimization"
        )
    require(
        "ri-mp2 cc-pvtz cc-pvtz/c" in rendered[0].lower(), "Invalid TZ RI-MP2"
    )
    require(
        "ri-mp2 cc-pvqz cc-pvqz/c" in rendered[1].lower(), "Invalid QZ RI-MP2"
    )
    require(
        "dlpno-ccsd(t) cc-pvtz cc-pvtz/c" in rendered[2].lower(),
        "DLPNO focal input lacks matching AuxC",
    )
    return {
        "state": "materialized_and_typed_arithmetic_not_executed",
        "project_yaml_sha256s": (
            *(
                project_sha256
                for _settings, project_sha256 in settings_and_projects
            ),
            dlpno_project_sha256,
        ),
        "routes": tuple(text.splitlines()[0] for text in rendered),
        "typed_components_hartree": _typed_focal_cbs(),
        "thermochemistry_requested": False,
    }


def qualify_iefpcm(water: Molecule) -> dict[str, object]:
    settings, project_sha256 = _settings_from_project_yaml(
        program="gaussian",
        stage="opt",
        settings_class=GaussianJobSettings,
        stage_settings={
            "functional": "B3LYP",
            "basis": "6-311+G(d,p)",
            "charge": 0,
            "multiplicity": 1,
            "freq": True,
            "solvent_model": "iefpcm",
            "solvent_id": "dichloromethane",
            "title": "Public water IEFPCM dichloromethane opt/freq surface",
        },
    )
    rendered = _gaussian_input(settings, water, label="water_dcm_opt_freq")
    route = next(
        line for line in rendered.splitlines() if line.startswith("#")
    )
    lower = route.lower()
    for token in (
        "opt",
        "freq",
        "b3lyp",
        "6-311+g(d,p)",
        "scrf=(iefpcm,solvent=dichloromethane)",
    ):
        require(token in lower, f"IEFPCM opt/freq route omitted {token}")

    return {
        "state": "materialized_not_executed",
        "project_yaml_sha256": project_sha256,
        "route": route,
        "thermochemistry_boundary": (
            "shared typed engine qualified on runtime-generated PySCF Hessian"
        ),
    }


def qualify_td(water: Molecule) -> dict[str, object]:
    common = {
        "functional": "wB97XD",
        "basis": "def2-SVP",
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "td",
        "freq": False,
        "states": "singlets",
        "root": 1,
        "nstates": 30,
        "title": "Public fixed-geometry water vertical TD surface",
    }
    gas, gas_project_sha256 = _settings_from_project_yaml(
        program="gaussian",
        stage="td",
        settings_class=GaussianTDDFTJobSettings,
        stage_settings=common,
    )
    acetone, acetone_project_sha256 = _settings_from_project_yaml(
        program="gaussian",
        stage="td",
        settings_class=GaussianTDDFTJobSettings,
        stage_settings={
            **common,
            "solvent_model": "smd",
            "solvent_id": "acetone",
        },
    )
    gas_text = _gaussian_input(gas, water, label="water_td_gas")
    acetone_text = _gaussian_input(acetone, water, label="water_td_acetone")
    routes = tuple(
        next(line for line in text.splitlines() if line.startswith("#"))
        for text in (gas_text, acetone_text)
    )
    gas_lower, acetone_lower = (route.lower() for route in routes)
    for route in (gas_lower, acetone_lower):
        for token in ("wb97xd", "def2svp", "td(singlets,nstates=30,root=1)"):
            require(token in route, f"Vertical TD route omitted {token}")
        require("tda" not in route, "Vertical TD request silently enabled TDA")
        require("freq" not in route, "Vertical TD request added frequencies")
        require("opt" not in route, "Vertical TD request changed the geometry")
    require("scrf" not in gas_lower, "Gas-phase TD input contains solvent")
    require(
        "scrf=(smd,solvent=acetone)" in acetone_lower,
        "Acetone TD input lost its SMD environment",
    )
    return {
        "state": "materialized_not_executed",
        "project_yaml_sha256s": (
            gas_project_sha256,
            acetone_project_sha256,
        ),
        "routes": routes,
    }


def qualify() -> dict[str, object]:
    nitric_oxide, water = _public_molecules()
    return {
        "schema_version": "chemsmart.cpu-program-surface-qualification.v1",
        "status": "qualified",
        "engine_executed": False,
        "materialization_boundary": "project_yaml_to_stage_settings_to_native_writer",
        "gaussian_radical_link": qualify_radical_link(nitric_oxide),
        "orca_doublet_focal": qualify_orca_focal(nitric_oxide),
        "gaussian_iefpcm_opt_freq": qualify_iefpcm(water),
        "gaussian_vertical_td": qualify_td(water),
    }


def main() -> int:
    print(
        json.dumps(
            qualify(),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
