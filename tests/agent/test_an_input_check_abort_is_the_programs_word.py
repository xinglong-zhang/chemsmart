"""An ORCA input-check abort is typed as the program's own refusal.

REACH-1 po3 (2026-09-06): six saddle-search launches died in ORCA's
input check within a second, on RIJK with an analytical Hessian and on
bare RI with a hybrid functional, and every one was typed
failed_nonconverged_geometry, so the repair menu told the next cycle to
restart from a geometry the run had never reached. The class was
undiagnosed and fell through; the pattern that names it existed.
"""

import pytest

from chemsmart.agent.terminal_states import _classify_failure
from chemsmart.io.native_failure import summarize_orca_native_failure

_BARE_RI = (
    "Error: RI is on but the HF exchange must be handled somehow",
    "Use one of the keywords",
    "     RIJDX  : treat the HF exchange exactly (equals RIJONX)",
    "     RIJCOSX: treat the HF exchange by chain-of-spheres",
    "     RIJK   : treat Coulomb and exchange both by RI",
    "",
    "[file orca_main/main_input_check.cpp, line 3594]: Error (ORCA_MAIN): "
    "... aborting the run",
)
_RIJK_HESSIAN = (
    "WARNING: Analytical Hessian not available with RIJK approximation",
    "  ===> : Skipping actual calculation",
    "[file orca_main/main_input_check.cpp, line 8637]: Error (ORCA_MAIN): "
    "... aborting the run",
)


@pytest.mark.capability("gate:terminal_state_vocabulary")
@pytest.mark.parametrize("tail", (_BARE_RI, _RIJK_HESSIAN))
def test_the_abort_is_classified_and_the_engines_line_is_quoted(tail):
    summary = summarize_orca_native_failure(tail)
    assert summary is not None
    assert summary.error_class == "input_check"
    quoted = "\n".join(summary.engine_lines)
    assert "main_input_check.cpp" in quoted
    assert "RIJ" in quoted
    assert (
        _classify_failure(
            jobtype="ts",
            findings=("orca.result.optimization_not_converged",),
            native_class=summary.error_class,
            converged=False,
            reached=None,
            planned=None,
        )
        == "failed_native"
    )
