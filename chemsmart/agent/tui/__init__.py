"""Interactive terminal surface for the production ChemSmart Agent."""

from __future__ import annotations

from pathlib import Path


def launch_tui(
    *,
    workspace: str | Path,
    secret_file: str | Path,
    provider: str | None = None,
    provider_config_file: str | Path | None = None,
    execution_envelope_file: str | Path | None = None,
    review_file: str | Path | None = None,
    identity_manifest: str | Path | None = None,
    analysis_completion_file: str | Path | None = None,
    plain: bool = False,
) -> None:
    """Launch the optional Textual application.

    Importing :mod:`chemsmart.agent.tui` never requires Textual.  This keeps
    the normal CLI usable when the ``agent-tui`` extra is not installed.
    """

    try:
        from .app import ChemSmartAgentApp
        from .controller import AgentSessionConfigV1, AgentTuiController
    except ImportError as exc:  # pragma: no cover - exercised without extra
        if exc.name and exc.name.split(".", 1)[0] in {"textual", "rich"}:
            raise RuntimeError(
                "ChemSmart Agent TUI is optional; install "
                "'chemsmart[agent-tui]'"
            ) from exc
        raise

    config = AgentSessionConfigV1(
        workspace=Path(workspace),
        secret_file=Path(secret_file),
        provider=provider,
        provider_config_file=(
            Path(provider_config_file)
            if provider_config_file is not None
            else None
        ),
        execution_envelope_file=(
            Path(execution_envelope_file).expanduser().resolve()
            if execution_envelope_file is not None
            else None
        ),
        review_file=(
            Path(review_file).expanduser().resolve()
            if review_file is not None
            else None
        ),
        identity_manifest=(
            Path(identity_manifest) if identity_manifest is not None else None
        ),
        analysis_completion_file=(
            Path(analysis_completion_file)
            if analysis_completion_file is not None
            else None
        ),
    )
    app = ChemSmartAgentApp(AgentTuiController(config), plain=plain)
    kwargs = {}
    if plain:
        app.animation_level = "none"
        kwargs = {"inline": True, "inline_no_clear": True, "mouse": False}
    app.run(**kwargs)
    summary = getattr(app, "exit_summary", None)
    if summary:
        line = (
            f"ChemSmart Agent · {summary['planning_sessions']} planning "
            f"session(s) · {summary['executions']} execution(s)"
        )
        for report in summary["report_paths"]:
            line += f"\n  report: {report}"
        print(line)


__all__ = ["launch_tui"]
