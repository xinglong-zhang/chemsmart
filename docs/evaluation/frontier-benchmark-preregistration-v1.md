# Frontier Benchmark Preregistration V1

## Status

`FORJADOR-Q03-ORGANIC-L2-V1` is sealed, held out, and not authorized for
provider or chemistry-engine execution. The machine-readable source of truth is
[`frontier-benchmark-preregistration-v1.json`](frontier-benchmark-preregistration-v1.json).
It contains source locators and byte hashes, but no ground-truth values.

The benchmark source is the CC BY 4.0 preprint [El Agente
Forjador](https://arxiv.org/abs/2604.14609v1) and version 1.0 of its CC0 1.0
[replication dataset](https://doi.org/10.5683/SP3/0YOMKL). The dataset metadata
was retrieved through the Scholars Portal Dataverse API. The registered paper
PDF and ten Q03 artifacts were also retrieved byte-for-byte and independently
SHA-256 hashed on 2026-08-04.

## Corrected Q03 source binding

The earlier draft mixed Organic Level 2 and Organic Level 3 oracle files. The
following binding is authoritative:

| Role | Dataverse file ID | File | SHA-256 |
|---|---:|---|---|
| Task prompt | 1155837 | `question.md` | `d9ebf071cb6e446c661742cf5d9c2eb52d730801afbbc62d873ef78e961f7056` |
| Source scorer | 1155824 | `score_groundtruth_organic_level2.py` | `d634b6bb814c4f5a33498b9b684952184553d181c200d9139b4218292cc88256` |
| Ground truth | 1155886 | `groundtruth_organic_level2.md` | `e202540f1a78515d17e389cd6d2a569dc55001cb2d55986770177e22548a9ab6` |
| Methodology rubric | 1155893 | `organic_molecule_analysis.md` | `814dbe437cbd6aaf847b818f593884c0d8e2e00dec515e22b8f40a74d5ebbc4f` |

Dataverse IDs 1155890 and 1155862 identify the Organic Level 3 scorer and
ground truth. They are excluded from Q03. The correct Level 2 scorer specifies
the HOMO--LUMO gap in Hartree, and the methodology rubric also specifies
Hartree. Therefore `unit_conflict_findings` is empty; the earlier proposed unit
conflict is not supported by the registered evidence.

The scorer, ground truth, and methodology rubric are host-oracle artifacts.
They must never enter a benchmark model context.

## Challenge

The model-visible inputs are the exact prompt and six supplied XYZ files:

| Order | System | Charge | Multiplicity | Dataverse file ID |
|---:|---|---:|---:|---:|
| 1 | caffeine | 0 | 1 | 1155856 |
| 2 | theobromine | 0 | 1 | 1155786 |
| 3 | aspirin | 0 | 1 | 1155785 |
| 4 | methyl salicylate | 0 | 1 | 1155815 |
| 5 | diisopropylamide anion | -1 | 1 | 1155794 |
| 6 | diisopropylammonium cation | +1 | 1 | 1155844 |

The fixed source conditions are gas-phase HF/def2-SVP. Each system requires
optimization and a Hessian-based stationary-point assessment, followed by
per-system coordinates, total energy, point group, dipole, orbital energies,
HOMO--LUMO gap, and Mulliken charges with their declared units.

## ChemSmart adaptation

The source task is adapted without importing its runtime tool-forging model:

- Exact XYZ bytes, charge, multiplicity, HF, def2-SVP, and gas phase are
  adopted.
- The stationary-point workflow becomes a project-YAML-backed `OPT -> HESS`
  command DAG with explicit artifact handoff.
- Reports are evidence-linked views over validated artifacts.
- An imaginary-frequency finding blocks for an explicit replan. It never
  triggers automatic coordinate mutation, retry, or resubmission.
- Model-authored scripts, runtime tool forging, and shell execution are
  rejected.
- Source accuracy and methodology scores remain separate from ChemSmart's
  all-required identity, project, DAG, and preview criteria.
- Forjador tool-reuse efficiency has no direct ChemSmart scoring analogue;
  provider cost remains a separately reported metric.

This preregistration tests future workflow planning and evidence composition.
It does not establish execution, reproduction, generality, or a comparison
with the published Forjador system.
