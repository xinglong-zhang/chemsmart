"""A declared project parameter must survive the trip into generated input.

The preview validator parses the input ChemSmart generated and compares it with
the settings that were requested.  A parameter the writer emits into the route
but no parser can read back therefore reads as a mismatch, and the model is
told a correct project produced invalid input.

Observed live: an n-butane session reproducing a protocol that calls for
conventional four-index integrals wrote ``ri_approximation: none``, was told
the generated input was invalid across both rotamers and both basis spellings,
and cleared the finding by deleting the key -- the one edit that lets ORCA's
default density fitting back in.  The harness argued the model out of the
scientifically correct input.

The single missing property is repaired in ``chemsmart/io/orca/route.py``.
What is pinned here is the contract, over every declared domain value, so the
next parameter added to the route cannot reopen it.
"""

import pytest

from chemsmart.io.orca.route import ORCARoute
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.settings.capabilities import PROGRAM_CAPABILITIES

#: Declared parameters that reach the generated input through a block rather
#: than the ``!`` route line, so the route parser is the wrong place to look.
#: ``reference`` becomes ``%scf HFTyp``; ``frozen_core`` becomes
#: ``%method FrozenCore``; IRC ``direction`` becomes ``%irc Direction``; and
#: the coupled excited-state settings become a ``%tddft`` block.  Those
#: settings cannot be instantiated or parsed as independent route tokens.
_BLOCK_EMITTED = {
    "direction",
    "frozen_core",
    "nstates",
    "reference",
    "response_method",
    "state_manifold",
}

#: Declared values the writer refuses on a non-relativistic basis, by design.
#: Pairing them with a recontracted basis is a separate concern from whether
#: the route round-trips.
_BASIS_CONSTRAINED = {"relativistic"}


def _domains():
    for name, values in PROGRAM_CAPABILITIES["orca"].project_parameter_domains:
        if name in _BLOCK_EMITTED or name in _BASIS_CONSTRAINED:
            continue
        for value in values:
            yield name, value


@pytest.mark.parametrize("name,value", list(_domains()))
def test_a_declared_route_parameter_reads_back_as_itself(name, value):
    project = {"basis": "6-31G*", name: value}
    if name != "ab_initio":
        project["functional"] = "B3LYP"
    elif str(value).casefold().startswith("dlpno-cc"):
        # The parameter-domain probe still has to construct a scientifically
        # admissible method.  ORCA local coupled cluster needs a correlation
        # fitting (AuxC) space; AutoAux is explicit here so this round-trip
        # test does not depend on any silent writer default.
        project["aux_basis"] = "AutoAux"
    settings = ORCAJobSettings.from_dict(project)
    route = ORCARoute(settings.route_string)
    recovered = getattr(route, name, None)
    assert recovered is not None, (
        f"the route carries {name}={value!r} as {settings.route_string!r} but "
        f"ORCARoute cannot read it back, so the preview validator will report "
        "a correct project as invalid"
    )
    assert str(recovered).lower() == str(value).lower()


def test_conventional_four_index_integrals_survive_the_round_trip():
    """The exact case that made a live session delete a correct key."""

    settings = ORCAJobSettings.from_dict(
        {
            "functional": "B3LYP",
            "basis": "6-31G*",
            "freq": True,
            "ri_approximation": "none",
        }
    )
    assert "NoRI" in settings.route_string
    assert ORCARoute(settings.route_string).ri_approximation == "none"


@pytest.mark.parametrize(
    "choice,keyword",
    [("none", "NoRI"), ("ri", "RI"), ("rijcosx", "RIJCOSX"), ("rijk", "RIJK")],
)
def test_settings_reader_preserves_ri_choice(tmp_path, choice, keyword):
    """The preview parser must retain the RI choice, not only see the token."""

    native_input = tmp_path / "round-trip.inp"
    native_input.write_text(
        f"! PBE def2-SVP {keyword}\n"
        "* xyz 0 1\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 0.74\n"
        "*\n",
        encoding="utf-8",
    )

    observed = ORCAJobSettings.from_filepath(str(native_input))
    assert observed.ri_approximation == choice


def test_an_absent_choice_reads_back_as_absent_not_as_a_default():
    """Reporting a default the project never set is the mirror-image bug."""

    settings = ORCAJobSettings.from_dict(
        {"functional": "B3LYP", "basis": "6-31G*"}
    )
    assert ORCARoute(settings.route_string).ri_approximation is None


def test_every_route_parameter_is_covered_by_this_contract():
    """A parameter added to the route later must be classified, not ignored."""

    declared = {
        name
        for name, _values in PROGRAM_CAPABILITIES[
            "orca"
        ].project_parameter_domains
    }
    covered = (
        {name for name, _ in _domains()} | _BLOCK_EMITTED | _BASIS_CONSTRAINED
    )
    assert declared <= covered, sorted(declared - covered)
