"""Resolve project settings from an exact YAML path.

The project YAML is ChemSmart's canonical control surface: it carries the
computational-chemistry rationale that the CLI then applies.  A project that
exists only as a *name* under ``~/.chemsmart/<program>/`` cannot be authored at
run time, so any caller that renders a project -- an agent, a test harness, a
scripted parameter sweep -- could not afterwards use the file it had just
written.  It had to install that file into the user's home directory first,
which makes an ephemeral, reproducible project impossible to express.

PySCF and xTB already resolved an exact path and keep their own fail-closed
ordering; ORCA and Gaussian did not.  This module holds the implementation they
share, so every program's ``from_project`` now accepts the same two things: a
path to a project file, or the name of an installed one.
"""

import os

__all__ = ["project_settings_from_path"]


def project_settings_from_path(project, manager_class):
    """Return settings loaded from ``project`` when it names an existing file.

    Returns ``None`` when ``project`` is not a usable filesystem path, leaving
    the caller's name-based resolution untouched.  Resolution therefore stays
    strictly additive: a name that resolved before still resolves the same way.

    Args:
        project: A path to a project YAML, a project name, or ``None``.
        manager_class: The program's settings manager, which must accept
            ``filename=`` and expose ``create()``.

    Returns:
        The loaded project settings, or ``None`` if ``project`` is not a file.
    """

    if project is None:
        return None
    try:
        path = os.fspath(project)
    except TypeError:
        # A non-path object is a name-like value; leave it to the caller.
        return None
    if not os.path.isfile(path):
        return None
    return manager_class(filename=path).create()


def require_jobtype_settings(settings, *, program, jobtype, project):
    """Return per-jobtype settings, or say which section is missing.

    A project that omits the section a job type reads used to return ``None``,
    and the CLI then called ``.merge()`` on it.  The caller saw

        AttributeError: 'NoneType' object has no attribute 'merge'

    which names neither the program, the job type, nor the section to add.  A
    live agent session hit this three times on a Gaussian excited-state job and
    concluded from it that the program could not preview that job at all --
    a wrong conclusion caused entirely by an opaque failure.
    """

    if settings is not None:
        return settings
    present = sorted(
        name.lstrip("_").removesuffix("_settings")
        for name, value in vars(project).items()
        if name.endswith("_settings") and value is not None
    )
    detail = (
        f"sections present: {present}"
        if present
        else "the project defines no job-type sections at all"
    )
    raise ValueError(
        f"{program} project settings define no '{jobtype}' section, so a "
        f"'{jobtype}' job has nothing to merge. Add a '{jobtype}:' section, "
        f"or a 'gas:' section that every job type inherits; {detail}."
    )
