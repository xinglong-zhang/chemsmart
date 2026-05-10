from chemsmart.agent.wizard import (
    ModuleSystem,
    ProgramFinding,
    Topology,
    detect_module_system,
    discover_conda,
    find_program,
    run_software_survey,
)
from chemsmart.agent.wizard.probe import ProbeResult


class StubRunner:
    def __init__(self, local_results=None, ssh_results=None):
        self.local_results = local_results or {}
        self.ssh_results = ssh_results or {}

    def run_local(self, command, timeout_s=15):
        return self.local_results.get(
            tuple(command),
            ProbeResult(
                command=" ".join(command),
                mode="local",
                host=None,
                returncode=1,
                stdout="",
                stderr="",
                duration_s=0.0,
                truncated=False,
            ),
        )

    def run_ssh(self, host, command, timeout_s=15):
        return self.ssh_results.get(
            (host, command),
            ProbeResult(
                command=command,
                mode="ssh",
                host=host,
                returncode=1,
                stdout="",
                stderr="",
                duration_s=0.0,
                truncated=False,
            ),
        )


def _result(command, stdout="", stderr="", returncode=0):
    return ProbeResult(
        command=command,
        mode="local",
        host=None,
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
        duration_s=0.0,
        truncated=False,
    )


def test_detect_module_system_identifies_lmod_from_version_output():
    runner = StubRunner(
        local_results={
            ("type", "module"): _result(
                "type module", "module is a function\n"
            ),
            ("module", "--version"): _result(
                "module --version",
                stderr="Modules based on Lua: Version 8.7.49  Lmod\n",
            ),
        }
    )

    assert detect_module_system(runner, Topology("A", "localhost", [])) == (
        ModuleSystem(
            kind="lmod",
            version="Modules based on Lua: Version 8.7.49  Lmod",
        )
    )


def test_find_program_prefers_executable_on_path():
    runner = StubRunner(
        local_results={
            ("command", "-v", "g16"): _result(
                "command -v g16", "/opt/gaussian/g16\n"
            ),
            ("readlink", "-f", "/opt/gaussian/g16"): _result(
                "readlink -f /opt/gaussian/g16",
                "/apps/gaussian/g16\n",
            ),
        }
    )

    assert find_program(
        runner,
        Topology("A", "localhost", []),
        "gaussian",
        ["g16", "g09"],
        ["gaussian", "g16", "g09"],
    ) == ProgramFinding(
        program="gaussian",
        exefolder="/apps/gaussian",
        source="path",
        module_candidates=[],
        on_path=True,
    )


def test_find_program_returns_single_module_candidate_when_not_on_path():
    runner = StubRunner(
        local_results={
            ("command", "-v", "orca"): _result(
                "command -v orca", returncode=1
            ),
            ("which", "orca"): _result("which orca", returncode=1),
            ("module", "-t", "avail"): _result(
                "module -t avail",
                stderr="orca/5.0.4\n",
            ),
        }
    )

    assert find_program(
        runner,
        Topology("A", "localhost", []),
        "orca",
        ["orca"],
        ["orca"],
    ) == ProgramFinding(
        program="orca",
        exefolder=None,
        source="module",
        module_candidates=["orca/5.0.4"],
        on_path=False,
    )


def test_find_program_sorts_multiple_module_candidates_by_length_then_name():
    runner = StubRunner(
        local_results={
            ("command", "-v", "nciplot"): _result(
                "command -v nciplot", returncode=1
            ),
            ("which", "nciplot"): _result("which nciplot", returncode=1),
            ("module", "-t", "avail"): _result(
                "module -t avail",
                stderr=("nci\n" "nciplot/1.0\n" "chem/nciplot\n"),
            ),
        }
    )

    assert find_program(
        runner,
        Topology("A", "localhost", []),
        "nciplot",
        ["nciplot"],
        ["nciplot", "nci"],
    ) == ProgramFinding(
        program="nciplot",
        exefolder=None,
        source="module",
        module_candidates=["nci", "nciplot/1.0", "chem/nciplot"],
        on_path=False,
    )


def test_discover_conda_returns_env_and_base():
    runner = StubRunner(
        local_results={
            ("printf", "%s\\n", "$CONDA_PREFIX"): _result(
                "printf %s\\n $CONDA_PREFIX",
                "/opt/conda/envs/chemsmart\n",
            ),
            ("conda", "info", "--base"): _result(
                "conda info --base",
                "/opt/conda\n",
            ),
        }
    )

    assert discover_conda(runner, Topology("A", "localhost", [])) == (
        "/opt/conda",
        "/opt/conda/envs/chemsmart",
    )


def test_run_software_survey_collects_all_programs():
    runner = StubRunner(
        local_results={
            ("type", "module"): _result(
                "type module", "module is a function\n"
            ),
            ("module", "--version"): _result(
                "module --version",
                stderr="Lmod 8.7.49\n",
            ),
            ("printf", "%s\\n", "$CONDA_PREFIX"): _result(
                "printf %s\\n $CONDA_PREFIX",
                "/opt/conda/envs/chemsmart\n",
            ),
            ("conda", "info", "--base"): _result(
                "conda info --base",
                "/opt/conda\n",
            ),
            ("command", "-v", "g16"): _result(
                "command -v g16", "/apps/bin/g16\n"
            ),
            ("readlink", "-f", "/apps/bin/g16"): _result(
                "readlink -f /apps/bin/g16", "/apps/bin/g16\n"
            ),
            ("command", "-v", "orca"): _result(
                "command -v orca", returncode=1
            ),
            ("which", "orca"): _result("which orca", returncode=1),
            ("command", "-v", "nciplot"): _result(
                "command -v nciplot", returncode=1
            ),
            ("which", "nciplot"): _result("which nciplot", returncode=1),
            ("module", "-t", "avail"): _result(
                "module -t avail",
                stderr="orca/5.0.4\nnciplot/1.0\n",
            ),
        }
    )

    survey = run_software_survey(runner, Topology("A", "localhost", []))

    assert survey.module_system == ModuleSystem(
        kind="lmod",
        version="Lmod 8.7.49",
    )
    assert survey.conda_base == "/opt/conda"
    assert survey.conda_env == "/opt/conda/envs/chemsmart"
    assert survey.programs["gaussian"].source == "path"
    assert survey.programs["orca"].module_candidates == ["orca/5.0.4"]
    assert survey.programs["nciplot"].module_candidates == ["nciplot/1.0"]
