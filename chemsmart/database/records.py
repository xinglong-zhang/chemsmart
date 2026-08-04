from dataclasses import dataclass


@dataclass(frozen=True)
class AssembledRecord:
    record_id: str
    meta: dict
    results: dict
    molecules: list
    provenance: dict
