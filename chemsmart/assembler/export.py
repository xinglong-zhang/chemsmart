import csv
import json


class DataExporter:
    def __init__(self, data, outputfile="database", keys=None):
        if keys is None:
            keys = [
                "filename",
                "program",
                "version",
                "date",
                "functional",
                "basis_set",
                "chemical_formula",
                "energy",
            ]
        self.data = data
        self.outputfile = outputfile
        self.keys = keys

    def _filter_data(self):
        filtered = []
        for d in self.data:
            entry = {k: d[k] for k in self.keys if k in d}
            filtered.append(entry)
        return filtered

    def to_json(self):
        filtered = self._filter_data()
        db = {str(i): row for i, row in enumerate(filtered, 1)}

        with open(self.outputfile + ".json", "w") as f:
            json.dump(db, f, indent=4)

    def to_csv(self):
        filtered = self._filter_data()
        with open(self.outputfile + ".csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=self.keys)
            writer.writeheader()
            writer.writerows(filtered)
