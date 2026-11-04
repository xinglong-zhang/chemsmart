"""Fixtures for chemsmart.utils.mixins tests."""

import os
import tempfile

import pytest
import yaml

from tests.helpers import make_dummy_file, make_mixin_instance


## fixtures for mixins
@pytest.fixture()
def temp_text_file():
    with tempfile.NamedTemporaryFile("w+", delete=False) as tmp:
        tmp.write("Line1\nLine2\n")
        tmp_name = tmp.name
    yield tmp_name
    os.remove(tmp_name)


@pytest.fixture()
def dummy_yaml_file():
    from chemsmart.utils.mixins import YAMLFileMixin

    class DummyYAMLFile(YAMLFileMixin):
        def __init__(self):
            self.filename = "dummy.yaml"
            self._content_lines_string = yaml.dump(
                {"key1": "value1", "key2": "value2"}
            )

        @property
        def content_lines_string(self):
            return self._content_lines_string

    return DummyYAMLFile()


@pytest.fixture()
def temp_folder_with_files():
    with tempfile.TemporaryDirectory() as tmpdir:
        file1 = os.path.join(tmpdir, "test1.txt")
        file2 = os.path.join(tmpdir, "test2.log")
        with open(file1, "w") as f:
            f.write("Test file 1")
        with open(file2, "w") as f:
            f.write("Test file 2")
        yield tmpdir, file1, file2


@pytest.fixture()
def make_frontier_orbital_file(temp_text_file):
    """Factory for a :class:`FileMixin` with the given ``multiplicity``
    (default 1) and alpha/beta occupied/virtual eigenvalue lists (default
    ``None``), for the frontier-orbital properties."""

    def _make(**overrides):
        attrs = dict(
            multiplicity=1,
            alpha_occ_eigenvalues=None,
            beta_occ_eigenvalues=None,
            alpha_virtual_eigenvalues=None,
            beta_virtual_eigenvalues=None,
        )
        attrs.update(overrides)
        return make_dummy_file(temp_text_file, **attrs)

    return _make


@pytest.fixture()
def dummy_orca_file():
    """:class:`ORCAFileMixin` with ``%mdci`` (cutoff/density) and SMD water
    ``%cpcm`` blocks, route ``! B3LYP def2-SVP``."""
    from chemsmart.utils.mixins import ORCAFileMixin

    return make_mixin_instance(
        ORCAFileMixin,
        contents=[
            "%mdci",
            "  cutoff 1e-5",
            "  density 1e-6",
            "%cpcm",
            "  smd true",
            '  solvent "water"',
        ],
        route_string="! B3LYP def2-SVP",
        filename="test.inp",
    )


@pytest.fixture()
def dummy_xtb_file():
    """:class:`XTBFileMixin` for a loose GFN2 ALPB(toluene) gradient opt."""
    from chemsmart.utils.mixins import XTBFileMixin

    return make_mixin_instance(
        XTBFileMixin,
        route_string=(
            "xtb p_benzyne.xyz --opt loose --gfn 2 --alpb toluene "
            "--chrg 0 --uhf 0 --grad"
        ),
        filename="test.out",
    )


@pytest.fixture()
def dummy_crest_file():
    """:class:`CRESTFileMixin` for a constrained GFN2 tight conformer
    search."""
    from chemsmart.utils.mixins import CRESTFileMixin

    return make_mixin_instance(
        CRESTFileMixin,
        route_string=(
            "crest 1a.xyz --cinp constraints.inp --gfn2 --chrg 0 --uhf 0 "
            "--optlev tight"
        ),
        filename="test.out",
    )
