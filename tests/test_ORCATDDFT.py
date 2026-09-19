"""Tests for the ORCA TDDFT / TDA implementation."""

import os
import re

import pytest
from click.testing import CliRunner

from chemsmart.cli.orca.orca import orca as orca_cli
from chemsmart.jobs.orca import ORCATDDFTJob
from chemsmart.jobs.orca.settings import ORCATDDFTJobSettings
from chemsmart.jobs.orca.tddft import ORCATDDFTJob as _ORCATDDFTJob_alias
from chemsmart.jobs.orca.writer import ORCAInputWriter
from chemsmart.settings.orca import ORCAProjectSettings


def _route_tokens(route_line):
    """Return the route tokens (without the leading ``!``) for exact-token
    checks — using substring matching would confuse ``Freq`` with ``NumFreq``.
    """
    stripped = route_line.lstrip()
    if stripped.startswith("!"):
        stripped = stripped[1:]
    return [tok for tok in stripped.split() if tok]


def _tokens_lower(route_line):
    return {tok.lower() for tok in _route_tokens(route_line)}


def _base_td_settings(project_name):
    project_settings = ORCAProjectSettings.from_project(project_name)
    td_settings = project_settings.td_settings()
    td_settings.charge = 0
    td_settings.multiplicity = 1
    return ORCATDDFTJobSettings(**td_settings.__dict__)


def _write_td_input(tmpdir, xyz_file, jobrunner, settings, label="orca_td"):
    job = ORCATDDFTJob.from_filename(
        filename=xyz_file,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
    )
    ORCAInputWriter(job=job).write(target_directory=tmpdir)
    return open(os.path.join(tmpdir, f"{label}.inp")).read()


class TestORCATDDFTWriterDefaults:
    def test_default_td_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Verify that a bare TD job writes only ``NRoots 3`` and ``TDA false``
        inside ``%tddft``, emits no ``%rel``, omits all unset TD options, and
        does not inherit ``Opt``/``Freq``/``NumFreq`` on the route line."""
        settings = _base_td_settings(orca_td_project_name)
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        assert "%tddft" in content
        assert "NRoots 3" in content
        assert "TDA false" in content
        assert "%rel" not in content
        # unset TD-only options must be omitted
        for keyword in (
            "Triplets",
            "DoSOC",
            "PrintLevel",
            "CPCMEQ",
            "DoNTO",
            "NTOStates",
            "NTOThresh",
            "MaxIter",
            "MaxDim",
            "ETol",
            "RTol",
            "TPrint",
        ):
            assert keyword not in content
        # A bare TD command must not inherit Opt/Freq/NumFreq.
        route_line = content.splitlines()[0]
        assert " Opt" not in route_line
        assert "OptTS" not in route_line
        assert "Freq" not in route_line
        assert "NumFreq" not in route_line

    @pytest.mark.parametrize(
        "field, keyword",
        [
            ("freq", "Freq"),
            ("numfreq", "NumFreq"),
        ],
    )
    def test_settings_preserve_explicit_freq_flag(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
        field,
        keyword,
    ):
        """Verify that ``ORCATDDFTJobSettings`` honours a Python-API caller
        who explicitly passes ``freq=True`` or ``numfreq=True`` (writing
        ``Freq`` / ``NumFreq`` on the route)."""
        project_settings = ORCAProjectSettings.from_project(
            orca_td_project_name
        )
        base = project_settings.td_settings()
        base.charge = 0
        base.multiplicity = 1
        base_dict = base.__dict__.copy()
        # freq and numfreq are mutually exclusive; clear both, then set one.
        base_dict["freq"] = False
        base_dict["numfreq"] = False
        base_dict[field] = True

        settings = ORCATDDFTJobSettings(**base_dict)
        assert getattr(settings, field) is True

        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        route_line = content.splitlines()[0]
        assert keyword in route_line

    def test_td_settings_preserve_inherited_freq_for_cli_to_resolve(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Verify that :class:`ORCATDDFTJobSettings` faithfully preserves an
        inherited ``freq=True`` (from ``ORCAJobSettings`` defaults or a
        parsed ``.log``): the settings class does not silently clear it.

        Clearing inherited Opt/Freq/NumFreq when the caller has not requested
        them via ``additional_route_parameters`` is the CLI layer's job (see
        the CLI-level tests below); this test only pins down the settings
        contract.
        """
        project_settings = ORCAProjectSettings.from_project(
            orca_td_project_name
        )
        base = project_settings.td_settings()
        # Simulate an inherited default that the CLI would otherwise clear.
        base.freq = True
        base.charge = 0
        base.multiplicity = 1
        settings = ORCATDDFTJobSettings(**base.__dict__)
        # Settings class preserves the caller's kwargs verbatim.
        assert settings.freq is True
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        # Documenting the contract: settings class does not scrub Freq.
        assert "Freq" in content.splitlines()[0]

    @pytest.mark.parametrize(
        "opt_excited, freq_flag, numfreq_flag, expected_tokens, forbidden",
        [
            (True, False, False, ["Opt"], ["Freq", "NumFreq"]),
            (True, True, False, ["Opt", "Freq"], ["NumFreq"]),
            (True, False, True, ["Opt", "NumFreq"], []),
            (False, True, False, ["Freq"], ["Opt", "NumFreq"]),
        ],
    )
    def test_route_reflects_opt_and_freq_flags(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
        opt_excited,
        freq_flag,
        numfreq_flag,
        expected_tokens,
        forbidden,
    ):
        """Verify that ``opt_excited`` + ``freq`` / ``numfreq`` project onto
        the ``!`` route line as ``Opt`` / ``Freq`` / ``NumFreq``."""
        settings = _base_td_settings(orca_td_project_name)
        settings.opt_excited = opt_excited
        settings.freq = freq_flag
        settings.numfreq = numfreq_flag
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        route_line = content.splitlines()[0]
        for token in expected_tokens:
            assert token in route_line, route_line
        for token in forbidden:
            # Freq is a substring of NumFreq; use word-ish matching by
            # requiring surrounding whitespace or end-of-line.
            pattern = re.compile(
                rf"(?<![A-Za-z]){re.escape(token)}(?![A-Za-z])"
            )
            assert not pattern.search(route_line), route_line

    def test_additional_route_parameters_appear_on_route_line(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Regression test: ``additional_route_parameters`` on the base
        ``ORCAJobSettings`` must be appended to the ``!`` route line for
        any ORCA job type (previously silently dropped)."""
        from chemsmart.jobs.orca import ORCASinglePointJob
        from chemsmart.settings.orca import ORCAProjectSettings

        project_settings = ORCAProjectSettings.from_project(
            orca_td_project_name
        )
        settings = project_settings.sp_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.additional_route_parameters = "TightSCF"
        job = ORCASinglePointJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_sp_extras",
            jobrunner=orca_jobrunner_no_scratch,
        )
        ORCAInputWriter(job=job).write(target_directory=tmpdir)
        content = open(os.path.join(tmpdir, "orca_sp_extras.inp")).read()
        route_line = content.splitlines()[0]
        assert "TightSCF" in route_line

    def test_additional_route_parameters_are_deduplicated(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Verify that repeated task tokens in ``additional_route_parameters``
        (case-insensitive) do not produce duplicates on the route line."""
        settings = _base_td_settings(orca_td_project_name)
        settings.opt_excited = True
        settings.freq = True
        # user redundantly asked for Opt + FREQ via extras
        settings.additional_route_parameters = "Opt FREQ TightSCF"
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        route_line = content.splitlines()[0]
        assert route_line.lower().count("opt") == 1
        assert route_line.lower().count("freq") == 1
        assert "TightSCF" in route_line


class TestORCATDDFTWriterSOC:
    def test_dosoc_without_soc_type_has_no_rel_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Verify that enabling SOC without ``soc_type`` writes ``DoSOC true``
        and ``Triplets true`` inside ``%tddft`` but omits the ``%rel`` block
        (no ``SOCType`` line either)."""
        settings = _base_td_settings(orca_td_project_name)
        settings.dosoc = True
        settings.triplets = True
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        assert "DoSOC true" in content
        assert "Triplets true" in content
        assert "%rel" not in content
        assert "SOCType" not in content

    def test_soc_type_emits_rel_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Verify that setting ``soc_type`` writes a minimal ``%rel`` block
        containing exactly ``SOCType 3`` alongside the ``%tddft`` block."""
        settings = _base_td_settings(orca_td_project_name)
        settings.dosoc = True
        settings.triplets = True
        settings.soc_type = 3
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        assert "%rel" in content
        assert re.search(r"%rel\s*\n\s*SOCType 3\s*\nend", content)


class TestORCATDDFTWriterNTO:
    def test_nto_written_inside_tddft_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Verify that ``NTOStates`` is written as a plain array parameter
        (no trailing ``end``) and that ``NTOThresh`` remains inside the same
        single ``%tddft`` block."""
        settings = _base_td_settings(orca_td_project_name)
        settings.donto = True
        settings.ntostates = [1, 2, 3]
        settings.ntothresh = 1e-4
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        # NTO lines must live inside the single %tddft block
        tddft_match = re.search(r"%tddft\s*\n(.*?)\nend", content, re.DOTALL)
        assert tddft_match is not None
        block = tddft_match.group(1)
        assert "DoNTO true" in block
        # NTOStates is an array parameter, not a nested block: no trailing ``end``.
        assert "NTOStates 1,2,3" in block
        assert "NTOStates 1,2,3 end" not in block
        # NTOThresh must still be inside the same %tddft block.
        assert "NTOThresh 0.0001" in block
        # only one %tddft/end pair
        assert content.count("%tddft") == 1
        assert content.count("\nend\n") >= 1


class TestORCATDDFTWriterOptionalParams:
    def test_explicit_false_preserved(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Verify that TD options explicitly set to ``False`` (``cpcmeq``,
        ``triplets``) are written as ``false`` rather than being omitted."""
        settings = _base_td_settings(orca_td_project_name)
        settings.cpcmeq = False
        settings.triplets = False
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        assert "CPCMEQ false" in content
        assert "Triplets false" in content

    def test_td_maxiter_does_not_leak_into_scf_block(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Verify that ``td_maxiter`` produces ``MaxIter`` inside ``%tddft``
        only, and does not cause an ``%scf`` block to be emitted."""
        settings = _base_td_settings(orca_td_project_name)
        settings.td_maxiter = 999
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        tddft_match = re.search(r"%tddft\s*\n(.*?)\nend", content, re.DOTALL)
        assert tddft_match is not None
        assert "MaxIter 999" in tddft_match.group(1)
        # no %scf block because scf_maxiter/scf_convergence are unset
        assert "%scf" not in content


class TestORCATDCLI:
    """CLI-level tests exercised through the ``orca td`` command.

    Job execution is patched out; each test only inspects the settings
    object handed to ``ORCATDDFTJob`` or the CLI exit code/output.
    """

    def _invoke(self, args):
        """Invoke the ``orca`` CLI group with an empty context object."""
        runner = CliRunner()
        return runner.invoke(orca_cli, args, obj={}, catch_exceptions=False)

    def test_td_command_is_registered(self):
        """Verify that the ``td`` subcommand is discoverable under ``orca``."""
        assert "td" in orca_cli.commands

    def test_td_default_flow_creates_job(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that ``orca td --nroots 10`` constructs an
        ``ORCATDDFTJobSettings`` with ``nroots=10``, ``tda=False``, ``freq=False``
        and leaves the SOC/NTO/triplets fields unset (``None``)."""
        mock_job = mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        mock_job.return_value = mocker.MagicMock()
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
                "--nroots",
                "10",
            ]
        )
        assert result.exit_code == 0, result.output
        assert mock_job.called
        settings = mock_job.call_args.kwargs["settings"]
        assert isinstance(settings, ORCATDDFTJobSettings)
        assert settings.nroots == 10
        assert settings.tda is False
        assert settings.dosoc is None
        assert settings.triplets is None
        assert settings.soc_type is None
        assert settings.freq is False
        assert settings.donto is None

    def test_td_dosoc_auto_enables_triplets(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that ``--dosoc`` without ``--triplets`` implicitly enables
        triplet excitations and that ``--printlevel`` is propagated, while
        ``soc_type`` stays ``None`` so no ``%rel`` block is requested."""
        mock_job = mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        mock_job.return_value = mocker.MagicMock()
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
                "--dosoc",
                "--printlevel",
                "3",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = mock_job.call_args.kwargs["settings"]
        assert settings.dosoc is True
        assert settings.triplets is True
        assert settings.printlevel == 3
        assert settings.soc_type is None

    def test_td_dosoc_and_no_triplets_conflict(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that ``--dosoc`` combined with an explicit ``--no-triplets``
        aborts the CLI with a non-zero exit code and an error mentioning
        ``--dosoc``."""
        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        runner = CliRunner()
        result = runner.invoke(
            orca_cli,
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
                "--dosoc",
                "--no-triplets",
            ],
            obj={},
        )
        assert result.exit_code != 0
        assert "--dosoc" in result.output

    def test_td_nto_states_implicitly_enables_nto(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that supplying ``--nto-states`` without ``--nto`` implicitly
        sets ``donto=True`` and parses the comma list into ``[1, 2, 3]``."""
        mock_job = mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        mock_job.return_value = mocker.MagicMock()
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
                "--nto-states",
                "1,2,3",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = mock_job.call_args.kwargs["settings"]
        assert settings.donto is True
        assert settings.ntostates == [1, 2, 3]

    def test_td_no_nto_with_states_conflicts(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that combining ``--no-nto`` with ``--nto-states`` aborts the
        CLI with a non-zero exit code and an error mentioning ``--no-nto``."""
        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        runner = CliRunner()
        result = runner.invoke(
            orca_cli,
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
                "--no-nto",
                "--nto-states",
                "1,2",
            ],
            obj={},
        )
        assert result.exit_code != 0
        assert "--no-nto" in result.output

    def test_td_requires_project_td_section(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
    ):
        """Verify that running ``orca td`` against a project YAML that has no
        ``td:`` section aborts with a non-zero exit code and an error that
        references ``td:``."""
        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        runner = CliRunner()
        result = runner.invoke(
            orca_cli,
            [
                "-p",
                os.path.basename(orca_yaml_settings_gas_solv_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
            ],
            obj={},
        )
        assert result.exit_code != 0
        assert "td:" in result.output

    def test_td_default_clears_inherited_freq_and_opt(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_inherited_project_name,
        orca_jobrunner_no_scratch,
        tmpdir,
    ):
        """Verify that a bare ``orca td`` clears Opt/Freq/NumFreq that would
        otherwise be inherited: the ``orca_td_inherited`` YAML omits ``freq``
        so ``ORCAJobSettings`` defaults leak ``freq=True`` into the merged
        TD settings, yet the CLI must produce a vertical TD.  The generated
        ``!`` route line is inspected via token matching (so ``Freq`` isn't
        mistaken for ``NumFreq``) and must still contain the ``%tddft`` block.
        """
        project_settings = ORCAProjectSettings.from_project(
            os.path.basename(orca_td_inherited_project_name)
        )
        # Sanity: the YAML really does hand us a settings object with
        # jobtype='opt' + freq=True inherited from defaults, and no extras.
        pre_cli = project_settings.td_settings()
        pre_cli.jobtype = "opt"  # simulate inherited task from a .log
        pre_cli.freq = True
        pre_cli.additional_route_parameters = None
        mocker.patch.object(
            project_settings, "td_settings", return_value=pre_cli
        )
        mocker.patch(
            "chemsmart.settings.orca.ORCAProjectSettings.from_project",
            return_value=project_settings,
        )

        captured = {}

        def _capture(**kw):
            captured["settings"] = kw["settings"]
            captured["label"] = kw["label"]
            captured["molecule"] = kw["molecule"]
            return mocker.MagicMock()

        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            side_effect=_capture,
        )

        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_inherited_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
            ]
        )
        assert result.exit_code == 0, result.output

        settings = captured["settings"]
        assert settings.opt_excited is False
        assert settings.freq is False
        assert settings.numfreq is False
        assert settings.jobtype == "td"
        assert settings.iroot is None
        assert settings.iroot_mult is None
        assert settings.follow_iroot is None

        # Feed the real settings through the writer to check the actual
        # input file: no Opt/Freq/NumFreq on the route line, %tddft present.
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        route_tokens = _tokens_lower(content.splitlines()[0])
        assert "opt" not in route_tokens
        assert "freq" not in route_tokens
        assert "numfreq" not in route_tokens
        assert "%tddft" in content

    def test_td_default_clears_inherited_numfreq(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_inherited_project_name,
        orca_jobrunner_no_scratch,
        tmpdir,
    ):
        """Same guarantee as the previous test, but with an inherited
        ``numfreq=True`` (e.g. from parsing an ORCA job that used numerical
        frequencies)."""
        project_settings = ORCAProjectSettings.from_project(
            os.path.basename(orca_td_inherited_project_name)
        )
        pre_cli = project_settings.td_settings()
        pre_cli.freq = False
        pre_cli.numfreq = True
        pre_cli.jobtype = "opt"
        pre_cli.additional_route_parameters = None
        mocker.patch.object(
            project_settings, "td_settings", return_value=pre_cli
        )
        mocker.patch(
            "chemsmart.settings.orca.ORCAProjectSettings.from_project",
            return_value=project_settings,
        )

        captured = {}

        def _capture(**kw):
            captured["settings"] = kw["settings"]
            return mocker.MagicMock()

        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            side_effect=_capture,
        )

        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_inherited_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
            ]
        )
        assert result.exit_code == 0, result.output

        settings = captured["settings"]
        assert settings.opt_excited is False
        assert settings.freq is False
        assert settings.numfreq is False

        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        route_tokens = _tokens_lower(content.splitlines()[0])
        assert "numfreq" not in route_tokens
        assert "freq" not in route_tokens
        assert "opt" not in route_tokens

    @pytest.mark.parametrize(
        "extras, expected_opt, expected_freq, expected_numfreq",
        [
            ("Opt", True, False, False),
            ("Freq", False, True, False),
            ("Opt Freq", True, True, False),
            ("opt numfreq", True, False, True),
        ],
    )
    def test_td_extras_activate_excited_state_task(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
        extras,
        expected_opt,
        expected_freq,
        expected_numfreq,
    ):
        """Verify that ``-r Opt|Freq|Opt Freq|opt numfreq`` on the ``orca``
        group flips the TD job into an excited-state task, populates a
        default target root/multiplicity, and consumes the tokens from
        ``additional_route_parameters``."""
        mock_job = mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        mock_job.return_value = mocker.MagicMock()
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                extras,
                "td",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = mock_job.call_args.kwargs["settings"]
        assert settings.opt_excited is expected_opt
        assert settings.freq is expected_freq
        assert settings.numfreq is expected_numfreq
        # additional_route_parameters must not still hold the consumed tokens
        assert settings.additional_route_parameters in (None, "")
        # default target: first singlet
        assert settings.iroot == 1
        assert settings.iroot_mult == "singlet"

    def test_td_extras_conflict_freq_and_numfreq(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that requesting both ``Freq`` and ``NumFreq`` in
        ``-r`` aborts the CLI."""
        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "Freq NumFreq",
                "td",
            ]
        )
        assert result.exit_code != 0
        assert "Freq" in result.output

    def test_td_extras_opt_with_dosoc_is_refused(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that ``-r Opt`` combined with ``--dosoc`` is refused with
        a clear message pointing to ``SOCGrad`` / ``orca inp``."""
        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "Opt",
                "td",
                "--dosoc",
            ]
        )
        assert result.exit_code != 0
        assert "SOCGrad" in result.output or "orca inp" in result.output

    def test_td_root_out_of_range(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that ``--root`` exceeding ``--nroots`` aborts the CLI."""
        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "Opt",
                "td",
                "--nroots",
                "3",
                "--root",
                "5",
            ]
        )
        assert result.exit_code != 0
        assert "--root" in result.output

    def test_td_root_mult_triplet_implies_triplets(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that ``--root-mult triplet`` sets ``iroot_mult`` and
        implicitly enables triplet excitations."""
        mock_job = mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        mock_job.return_value = mocker.MagicMock()
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "Opt",
                "td",
                "--root-mult",
                "triplet",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = mock_job.call_args.kwargs["settings"]
        assert settings.iroot_mult == "triplet"
        assert settings.triplets is True

    def test_td_root_mult_triplet_conflicts_with_no_triplets(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that ``--root-mult triplet`` + ``--no-triplets`` aborts."""
        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "Opt",
                "td",
                "--root-mult",
                "triplet",
                "--no-triplets",
            ]
        )
        assert result.exit_code != 0
        assert "triplet" in result.output.lower()

    def test_td_follow_root_requires_excited_task(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_project_name,
    ):
        """Verify that ``--follow-root`` without ``-r Opt``/``-r Freq``
        aborts (it only makes sense during excited-state optimizations)."""
        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
                "--follow-root",
            ]
        )
        assert result.exit_code != 0
        assert "--follow-root" in result.output

    def test_td_writer_emits_iroot_when_set(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_project_name,
        orca_jobrunner_no_scratch,
    ):
        """Verify that ``iroot``/``iroot_mult``/``follow_iroot`` land inside
        the single ``%tddft`` block when set."""
        settings = _base_td_settings(orca_td_project_name)
        settings.opt_excited = True
        settings.iroot = 2
        settings.iroot_mult = "singlet"
        settings.follow_iroot = True
        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        tddft_match = re.search(r"%tddft\s*\n(.*?)\nend", content, re.DOTALL)
        assert tddft_match is not None
        block = tddft_match.group(1)
        assert "IRoot 2" in block
        assert "IRootMult singlet" in block
        assert "FollowIRoot true" in block
        # Route line contains Opt because opt_excited is True.
        assert "Opt" in content.splitlines()[0]

    # --- YAML additional_route_parameters × CLI -r override matrix --------

    def _capture_td_settings(self, mocker):
        """Return a tuple ``(patch_target, captured_dict)`` for tests that
        want to inspect the settings the CLI hands to :class:`ORCATDDFTJob`
        without mocking the settings/writer transformation itself.
        """
        captured = {}

        def _capture(**kw):
            captured["settings"] = kw["settings"]
            captured["molecule"] = kw["molecule"]
            captured["label"] = kw["label"]
            return mocker.MagicMock()

        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            side_effect=_capture,
        )
        return captured

    def test_td_yaml_sets_task_flags_when_no_r(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_extras_project_name,
    ):
        """Verify that YAML ``additional_route_parameters: "Opt Freq
        TightSCF"`` on the ``td:`` section drives ``opt_excited``/``freq``
        and leaves ``TightSCF`` as the remaining extra, without ``-r``."""
        captured = self._capture_td_settings(mocker)
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_extras_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = captured["settings"]
        assert settings.opt_excited is True
        assert settings.freq is True
        assert settings.numfreq is False
        assert settings.additional_route_parameters == "TightSCF"
        # Default excited-state target populated for Opt/Freq tasks.
        assert settings.iroot == 1
        assert settings.iroot_mult == "singlet"

    def test_td_r_replaces_yaml_extras(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_extras_project_name,
    ):
        """Verify that ``-r Opt`` fully replaces the YAML value: no ``Freq``
        or ``TightSCF`` leaks through."""
        captured = self._capture_td_settings(mocker)
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_extras_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "Opt",
                "td",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = captured["settings"]
        assert settings.opt_excited is True
        assert settings.freq is False
        assert settings.numfreq is False
        assert settings.additional_route_parameters in (None, "")

    def test_td_r_keeps_only_non_task_extras(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_extras_project_name,
    ):
        """Verify that ``-r TightSCF`` keeps the calculation vertical (no
        inherited ``Opt``/``Freq`` from the YAML) and only propagates
        ``TightSCF`` as an extra route token."""
        captured = self._capture_td_settings(mocker)
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_extras_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "TightSCF",
                "td",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = captured["settings"]
        assert settings.opt_excited is False
        assert settings.freq is False
        assert settings.numfreq is False
        assert settings.additional_route_parameters == "TightSCF"

    def test_td_r_empty_clears_yaml_extras(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_extras_project_name,
    ):
        """Verify that ``-r ""`` explicitly clears the YAML value; the
        resulting job is a vertical TD with no route extras."""
        captured = self._capture_td_settings(mocker)
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_extras_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "",
                "td",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = captured["settings"]
        assert settings.opt_excited is False
        assert settings.freq is False
        assert settings.numfreq is False
        assert settings.additional_route_parameters in (None, "")

    def test_td_yaml_numfreq_task(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_numfreq_project_name,
    ):
        """Verify that ``NumFreq`` from the YAML activates ``numfreq`` and
        keeps ``freq``/``opt_excited`` off."""
        captured = self._capture_td_settings(mocker)
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_numfreq_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = captured["settings"]
        assert settings.numfreq is True
        assert settings.freq is False
        assert settings.opt_excited is False
        # Excited-state target still populated because a task was requested.
        assert settings.iroot == 1

    def test_td_yaml_freq_numfreq_conflict_is_rejected(
        self,
        mocker,
        single_molecule_xyz_file,
        orca_td_numfreq_project_name,
    ):
        """Verify that if the effective ``additional_route_parameters``
        (whether from YAML or ``-r``) contains both ``Freq`` and ``NumFreq``,
        the CLI aborts with the existing conflict guard."""
        mocker.patch(
            "chemsmart.jobs.orca.tddft.ORCATDDFTJob",
            autospec=True,
        )
        # -r replaces the YAML value; supplying both keywords through -r is
        # sufficient to exercise the same guard used against a bad YAML.
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_numfreq_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-r",
                "Freq NumFreq",
                "td",
            ]
        )
        assert result.exit_code != 0
        assert "Freq" in result.output

    def test_td_yaml_extras_end_to_end_writer(
        self,
        mocker,
        tmpdir,
        single_molecule_xyz_file,
        orca_td_extras_project_name,
        orca_jobrunner_no_scratch,
    ):
        """CLI → writer integration: YAML supplies ``Opt Freq TightSCF``; the
        CLI must produce settings that write a route line with ``Opt``,
        ``Freq`` and ``TightSCF`` (token-checked so ``Freq`` isn't confused
        with ``NumFreq``), plus a single ``%tddft`` block containing
        ``IRoot``/``IRootMult`` matching the default first-singlet target.
        """
        captured = self._capture_td_settings(mocker)
        result = self._invoke(
            [
                "-p",
                os.path.basename(orca_td_extras_project_name),
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
                "--nroots",
                "5",
            ]
        )
        assert result.exit_code == 0, result.output
        settings = captured["settings"]
        assert settings.opt_excited is True
        assert settings.freq is True

        content = _write_td_input(
            tmpdir,
            single_molecule_xyz_file,
            orca_jobrunner_no_scratch,
            settings,
        )
        route_line = content.splitlines()[0]
        tokens_lower = _tokens_lower(route_line)
        assert "opt" in tokens_lower
        assert "freq" in tokens_lower
        # NumFreq must not sneak in; check exact token match.
        assert "numfreq" not in tokens_lower
        # TightSCF preserved as an extra route token.
        assert any(
            tok.lower() == "tightscf" for tok in _route_tokens(route_line)
        )

        # Exactly one %tddft block, with the default first-singlet target.
        assert content.count("%tddft") == 1
        tddft_match = re.search(r"%tddft\s*\n(.*?)\nend", content, re.DOTALL)
        assert tddft_match is not None
        block = tddft_match.group(1)
        assert "NRoots 5" in block
        assert "TDA false" in block
        assert "IRoot 1" in block
        assert "IRootMult singlet" in block


def test_tddft_module_exports_the_same_class():
    """Regression check: the ``ORCATDDFTJob`` re-exported from
    ``chemsmart.jobs.orca`` is the same class object as the one defined in
    ``chemsmart.jobs.orca.tddft``."""
    assert ORCATDDFTJob is _ORCATDDFTJob_alias


def test_orcatd_registered_in_jobrunner():
    """Verify that ``"orcatd"`` is registered in ``ORCAJobRunner.JOBTYPES`` so
    the shared ORCA runner will accept TD jobs."""
    from chemsmart.jobs.orca.runner import ORCAJobRunner

    assert "orcatd" in ORCAJobRunner.JOBTYPES
