from chemsmart.agent.harness.workflow_state import project_name_from_request


def test_project_name_from_request_keeps_explicit_name() -> None:
    assert project_name_from_request("using the co2 project") == "co2"
    assert project_name_from_request("use project name carbonyl") == "carbonyl"
    assert (
        project_name_from_request("Using ORCA project orca_demo, run a scan.")
        == "orca_demo"
    )


def test_project_name_from_request_ignores_generic_workspace_state_words() -> (
    None
):
    assert project_name_from_request("using the loaded project") == ""
    assert project_name_from_request("with the current project") == ""
    assert (
        project_name_from_request(
            "keep the active Gaussian project and do not edit it"
        )
        == ""
    )
    assert project_name_from_request("keep the loaded project unchanged") == ""
    assert (
        project_name_from_request(
            "Use the loaded Gaussian QM/MM project settings for the active site."
        )
        == ""
    )
