import logging
from functools import cached_property

from chemsmart.assembler.gaussian import GaussianAssembler
from chemsmart.assembler.orca import ORCAAssembler
from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)


class SingleFileAssembler:
    def __init__(self, filename, index="-1"):
        self.filename = filename
        self.index = index

    @cached_property
    def assemble_data(self):
        assembler = self._get_assembler(self.filename)
        try:
            data = assembler.assemble()
        except Exception as e:
            logger.error(f"Error assembling {self.filename}: {e}")
            return None
        return data

    def _get_assembler(self, file):
        program = get_program_type_from_file(self.filename)
        if program == "gaussian":
            assembler = GaussianAssembler(file, index=self.index)
        elif program == "orca":
            assembler = ORCAAssembler(file, index=self.index)
        else:
            raise ValueError(
                "Unsupported format. Only 'gaussian' and 'orca' are supported."
            )
        return assembler

    def query(self, key, value):
        results = []
        record = self.assemble_data
        d_dict = {
            "record_id": record.record_id,
            "meta": record.meta,
            "results": record.results,
            "molecules": record.molecules,
            "provenance": record.provenance,
        }

        if key in d_dict and d_dict[key] == value:
            results.append(d_dict)
        elif key in d_dict["meta"] and d_dict["meta"][key] == value:
            results.append(d_dict)
        elif key in d_dict["results"] and d_dict["results"][key] == value:
            results.append(d_dict)

        return results
