import csv
import json
import sqlite3
from typing import Any, Dict, List

import numpy as np

from chemsmart.assembler.records import AssembledRecord


class DataExporter:
    """Utility to export assembled records to JSON/CSV/SQLite."""

    def __init__(self, data, outputfile=None, keys=None):
        self.data = data
        self.outputfile = outputfile
        self.keys = (
            list(keys)
            if isinstance(keys, (list, tuple))
            else (
                [k.strip() for k in keys.split(",")]
                if isinstance(keys, str)
                else None
            )
        )

    def _filter_data(self):
        """Filter keys if provided."""
        filtered = []
        for d in self.data:
            if isinstance(d, AssembledRecord):
                d_dict = {
                    "record_id": d.record_id,
                    "meta": d.meta,
                    "results": d.results,
                    "molecules": d.molecules,
                    "provenance": d.provenance,
                }
            else:
                d_dict = d
            if not self.keys:
                filtered.append(d_dict)
            else:
                entry = {}
                for k in self.keys:
                    if k in d_dict:
                        entry[k] = d_dict[k]
                    elif "meta" in d_dict and k in d_dict["meta"]:
                        entry[k] = d_dict["meta"][k]
                    elif "results" in d_dict and k in d_dict["results"]:
                        entry[k] = d_dict["results"][k]
                filtered.append(entry)
        return filtered

    def _convert(self, obj):
        """Recursively convert NumPy and other non-JSON-native types."""
        # NumPy arrays -> lists
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        # NumPy scalars -> Python scalars
        elif isinstance(obj, np.generic):
            return obj.item()
        # Mapping types
        elif isinstance(obj, dict):
            return {k: self._convert(v) for k, v in obj.items()}
        # Sequence types
        elif isinstance(obj, (list, tuple, set)):
            return [self._convert(x) for x in obj]
        return obj

    def to_json(self):
        filtered = self._filter_data()
        converted = self._convert(filtered)
        with open(self.outputfile + ".json", "w") as f:
            json.dump(converted, f, indent=4)

    def to_csv(self):
        filtered = self._filter_data()
        with open(self.outputfile + ".csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self.keys)
            writer.writeheader()
            writer.writerows(filtered)

    def to_sqlite(self, table_name="records"):
        filtered: List[Dict[str, Any]] = self._filter_data()
        if len(filtered) == 0:
            return
        db_path = self.outputfile + ".sqlite"
        conn = sqlite3.connect(db_path)
        try:
            cols_def = ", ".join(f"{k} TEXT" for k in self.keys)
            conn.execute(
                f"CREATE TABLE IF NOT EXISTS {table_name} (id INTEGER PRIMARY KEY AUTOINCREMENT, {cols_def})"
            )
            placeholders = ",".join("?" for _ in self.keys)
            insert_sql = f"INSERT INTO {table_name} ({', '.join(self.keys)}) VALUES ({placeholders})"
            rows = []
            for row in filtered:
                rows.append([str(row.get(k, "")) for k in self.keys])
            conn.executemany(insert_sql, rows)
            conn.commit()
        finally:
            conn.close()
