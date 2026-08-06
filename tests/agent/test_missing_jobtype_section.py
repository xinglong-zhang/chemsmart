"""A job type whose project section is missing must say which section.

Observed live: an agent session planning a Gaussian excited-state calculation
wrote a project with only a ``gas:`` section.  Gaussian reads ``td`` from its
own section, so the settings were ``None`` and the CLI called ``.merge()`` on
them:

    AttributeError: 'NoneType' object has no attribute 'merge'

The model tried three project variants, hit the same crash each time, and
concluded that the host "is not preview-conformant" for Gaussian TD -- a wrong
conclusion produced entirely by an opaque failure.  An earlier session had
worked around this by adding the section to a probe rather than repairing the
message, which is why it was still here to be hit.
"""

import pytest

from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.orca import ORCAProjectSettings
from chemsmart.settings.project_resolution import require_jobtype_settings


def _project(path, text):
    path.write_text(text)
    return path


@pytest.mark.parametrize(
    "manager,program,body,jobtype",
    [
        (
            GaussianProjectSettings,
            "gaussian",
            "gas:\n  basis: aug-cc-pVDZ\n  functional: B3LYP\n",
            "td",
        ),
        (
            ORCAProjectSettings,
            "orca",
            "gas:\n  basis: def2-SVP\n  functional: B3LYP\n",
            "td",
        ),
    ],
)
def test_a_missing_section_names_itself_and_the_remedy(
    tmp_path, manager, program, body, jobtype
):
    project = manager.from_project(
        str(_project(tmp_path / f"{program}.yaml", body))
    )
    with pytest.raises(ValueError) as failure:
        getattr(project, f"{jobtype}_settings")()
    message = str(failure.value)
    assert program in message
    assert f"'{jobtype}' section" in message
    assert "gas" in message, "the inherited-section remedy must be stated"
    assert "sections present" in message


def test_the_message_lists_what_the_project_did_define(tmp_path):
    project = GaussianProjectSettings.from_project(
        str(
            _project(
                tmp_path / "g.yaml",
                "gas:\n  basis: 6-31G*\n  functional: B3LYP\n",
            )
        )
    )
    message = str(pytest.raises(ValueError, project.td_settings).value)
    for present in ("opt", "sp", "irc"):
        assert present in message, present


def test_a_section_that_exists_is_returned_unchanged(tmp_path):
    project = GaussianProjectSettings.from_project(
        str(
            _project(
                tmp_path / "g.yaml",
                "gas:\n  basis: 6-31G*\n  functional: B3LYP\n",
            )
        )
    )
    assert project.opt_settings() is not None


def test_a_project_with_no_sections_at_all_says_so():
    class _Empty:
        pass

    with pytest.raises(ValueError, match="no job-type sections at all"):
        require_jobtype_settings(
            None, program="orca", jobtype="td", project=_Empty()
        )
