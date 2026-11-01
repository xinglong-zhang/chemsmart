"""Shared pytest configuration.

Fixtures are defined in ``tests/fixtures/<area>.py`` and registered here
via ``pytest_plugins``, so every test can request them by name without
importing. Plain stub classes and helper functions live in
``tests/helpers.py`` and are imported directly.
"""

pytest_plugins = [
    "tests.fixtures.io",
    "tests.fixtures.thermochemistry",
    "tests.fixtures.cli",
    "tests.fixtures.gaussian",
    "tests.fixtures.orca",
    "tests.fixtures.common",
    "tests.fixtures.xtb",
    "tests.fixtures.crest",
    "tests.fixtures.structures",
    "tests.fixtures.jobrunners",
    "tests.fixtures.molecules",
    "tests.fixtures.grouper",
    "tests.fixtures.mixins",
    "tests.fixtures.logs",
    "tests.fixtures.iterate",
    "tests.fixtures.database",
    "tests.fixtures.jobs",
    "tests.fixtures.pymol",
    "tests.fixtures.pka",
]
