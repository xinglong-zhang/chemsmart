from chemsmart import __version__ as chemsmart_version
from chemsmart.assembler.utils import file_size, sha256_content, utcnow_iso
from chemsmart.utils.io import get_program_type_from_file


def build_provenance(filename, output):
    return {
        "source_file": filename,
        "source_file_hash": sha256_content(output),
        "source_file_size": file_size(filename),
        "source_file_date": output.date,
        "program": get_program_type_from_file(filename),
        "program_version": output.version,
        "parser": output.__class__.__name__,
        "chemsmart_version": chemsmart_version,
        "assembled_at": utcnow_iso(),
    }
