"""The reader interface the host sensors assume is a contract, not a habit.

``_basin_sensor_inputs``, ``_same_structure_observations`` and the spin
observation read ``output.molecule``, ``final_energy``, ``converged`` and
``vibrational_frequencies`` behind bare ``except`` clauses that return
empty.  ``PySCFOutput`` had ``get_molecule()`` and none of the four, so
every sensor answered nothing for PySCF and called it an absence.  A
reader that cannot answer these names is a reader the sensors cannot
see, and the ladder must be able to say so.
"""

from __future__ import annotations

import importlib

import pytest

from chemsmart.analysis.result_readers import RESULT_READERS

#: What every native-output reader answers.  The xyz reader is data, not a
#: program run, and is held to none of this.
_SENSOR_ATTRIBUTES = (
    "molecule",
    "normal_termination",
    "vibrational_frequencies",
    "charge",
    "multiplicity",
)
#: What every reader of an executable optimising program answers as well.
_OPTIMISER_ATTRIBUTES = ("converged", "final_energy")
_EXECUTABLE = ("orca", "pyscf", "xtb")


def _output_class(reader):
    module, name = reader.parser_id.rsplit(".", 1)
    return getattr(importlib.import_module(module), name)


@pytest.mark.capability("selector:*")
def test_every_native_reader_answers_the_sensor_interface():
    missing = {}
    for program, reader in sorted(RESULT_READERS.items()):
        if program == "xyz":
            continue
        klass = _output_class(reader)
        wanted = _SENSOR_ATTRIBUTES + (
            _OPTIMISER_ATTRIBUTES if program in _EXECUTABLE else ()
        )
        absent = [name for name in wanted if not hasattr(klass, name)]
        if absent:
            missing[program] = absent
    assert not missing, missing
