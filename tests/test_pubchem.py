"""Tests for chemsmart.io.molecules.pubchem."""

import json

import pytest
from requests.exceptions import ConnectionError as RequestsConnectionError
from requests.exceptions import HTTPError as RequestsHTTPError
from requests.exceptions import Timeout

from chemsmart.io.molecules.pubchem import (
    _pubchem_2d_to_3d,
    _retryable_pubchem,
    pubchem_search,
    search_pubchem_raw,
)


def _make_response(text=None, status_code=200, raise_exc=None):
    response = MockResponse(text=text, status_code=status_code)
    if raise_exc is not None:
        response.raise_for_status = lambda: (_ for _ in ()).throw(raise_exc)
    return response


class MockResponse:
    def __init__(self, text=None, status_code=200):
        self.text = text
        self.status_code = status_code

    def raise_for_status(self):
        pass


def _http_error(status_code):
    error = RequestsHTTPError(f"HTTP {status_code}")
    error.response = MockResponse(status_code=status_code)
    return error


@pytest.fixture(autouse=True)
def clear_pubchem_cache():
    """pubchem_search is lru_cache'd; avoid cross-test cache pollution."""
    pubchem_search.cache_clear()
    yield
    pubchem_search.cache_clear()


class TestRetryablePubchem:
    def test_timeout_is_retryable(self):
        assert _retryable_pubchem(Timeout()) is True

    def test_http_error_with_retryable_status_is_retryable(self):
        error = _http_error(503)
        assert _retryable_pubchem(error) is True

    def test_http_error_with_non_retryable_status_is_not_retryable(self):
        error = _http_error(404)
        assert _retryable_pubchem(error) is False

    def test_http_error_without_response_is_not_retryable(self):
        error = RequestsHTTPError("boom")
        error.response = None
        assert _retryable_pubchem(error) is False

    def test_generic_request_exception_is_retryable(self):
        assert _retryable_pubchem(RequestsConnectionError()) is True

    def test_non_request_exception_is_not_retryable(self):
        assert _retryable_pubchem(ValueError("boom")) is False


class TestSearchPubchemRawConformersField:
    def test_conformers_field_returns_text_directly(self, mocker):
        mock_get = mocker.patch(
            "chemsmart.io.molecules.pubchem.requests.get",
            return_value=MockResponse(text="CONFORMER_SDF"),
        )
        result = search_pubchem_raw("abc123", "conformers")
        assert result == "CONFORMER_SDF"
        mock_get.assert_called_once()
        assert "conformers/abc123" in mock_get.call_args[0][0]


class TestSearchPubchemRawStandardField:
    def test_raises_on_http_error(self, mocker):
        mocker.patch(
            "chemsmart.io.molecules.pubchem.requests.get",
            return_value=_make_response(
                raise_exc=_http_error(404), status_code=404
            ),
        )
        with pytest.raises(RequestsHTTPError):
            search_pubchem_raw("123", "cid")

    def test_returns_text_when_conformer_lookup_fails_silently(self, mocker):
        # First call: main data succeeds.
        # Second call: nested conformer lookup fails and is suppressed.
        mocker.patch(
            "chemsmart.io.molecules.pubchem.requests.get",
            side_effect=[
                MockResponse(text="MAIN_SDF_DATA"),
                _make_response(raise_exc=_http_error(404), status_code=404),
            ],
        )
        result = search_pubchem_raw("123", "cid")
        assert result == "MAIN_SDF_DATA"

    def test_warns_when_multiple_conformers_found(self, mocker):
        conformer_json = json.dumps({"Conformers": ["c1", "c2"]})
        mocker.patch(
            "chemsmart.io.molecules.pubchem.requests.get",
            side_effect=[
                MockResponse(text="MAIN_SDF_DATA"),
                MockResponse(text=conformer_json),
                # Nested conformer-of-conformer lookup fails; suppressed.
                _make_response(raise_exc=_http_error(404), status_code=404),
            ],
        )
        with pytest.warns(UserWarning, match="more than one"):
            result = search_pubchem_raw("123", "cid")
        assert result == "MAIN_SDF_DATA"


class TestPubchemSearch:
    def test_success_via_3d(self, mocker):
        mocker.patch(
            "chemsmart.io.molecules.pubchem.analyze_input",
            return_value=("search-3d", "cid"),
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.search_pubchem_raw",
            return_value="RAW_SDF_3D",
        )
        fake_molecule = mocker.Mock()
        mocker.patch(
            "chemsmart.utils.utils.sdf2molecule",
            return_value=fake_molecule,
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Molecule", type(fake_molecule)
        )

        result = pubchem_search(cid="search-3d")
        assert result is fake_molecule

    def test_falls_back_to_2d_on_404(self, mocker):
        mocker.patch(
            "chemsmart.io.molecules.pubchem.analyze_input",
            return_value=("search-fallback", "cid"),
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.search_pubchem_raw",
            side_effect=[_http_error(404), "RAW_SDF_2D"],
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem._pubchem_2d_to_3d",
            return_value="RAW_SDF_3D_CONVERTED",
        )
        fake_molecule = mocker.Mock()
        mocker.patch(
            "chemsmart.utils.utils.sdf2molecule",
            return_value=fake_molecule,
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Molecule", type(fake_molecule)
        )

        result = pubchem_search(cid="search-fallback")
        assert result is fake_molecule

    def test_unexpected_http_error_retries_then_raises_retry_error(
        self, mocker
    ):
        mocker.patch("time.sleep")
        mocker.patch(
            "chemsmart.io.molecules.pubchem.analyze_input",
            return_value=("search-retry", "cid"),
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.search_pubchem_raw",
            side_effect=_http_error(500),
        )
        import tenacity

        with pytest.raises(tenacity.RetryError):
            pubchem_search(cid="search-retry")

    def test_2d_fallback_also_fails_returns_none_when_fail_silently(
        self, mocker
    ):
        mocker.patch(
            "chemsmart.io.molecules.pubchem.analyze_input",
            return_value=("search-bothfail", "cid"),
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.search_pubchem_raw",
            side_effect=[_http_error(404), _http_error(404)],
        )
        result = pubchem_search(cid="search-bothfail", fail_silently=True)
        assert result is None

    def test_timeout_retries_then_raises_retry_error(self, mocker):
        mocker.patch("time.sleep")
        mocker.patch(
            "chemsmart.io.molecules.pubchem.analyze_input",
            return_value=("search-timeout", "cid"),
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.search_pubchem_raw",
            side_effect=Timeout("timed out"),
        )
        import tenacity

        with pytest.raises(tenacity.RetryError):
            pubchem_search(cid="search-timeout")

    def test_2d_fallback_also_fails_raises_when_not_fail_silently(
        self, mocker
    ):
        mocker.patch(
            "chemsmart.io.molecules.pubchem.analyze_input",
            return_value=("search-bothfail-raise", "cid"),
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.search_pubchem_raw",
            side_effect=[_http_error(404), _http_error(404)],
        )
        with pytest.raises(ValueError, match="Error getting structure"):
            pubchem_search(cid="search-bothfail-raise", fail_silently=False)


class TestPubchem2dTo3d:
    def test_raises_import_error_when_rdkit_unavailable(self, mocker):
        mocker.patch("chemsmart.io.molecules.pubchem.rdkit", None)
        with pytest.raises(ImportError, match="RDKit package is required"):
            _pubchem_2d_to_3d("SOME_SDF_DATA")

    def test_multi_component_raises_runtime_error(self, mocker):
        mock_mol = mocker.Mock()
        mock_mol.GetPropsAsDict.return_value = {"PUBCHEM_COMPONENT_COUNT": 2}
        mock_mol.GetProp.return_value = "multi_component_name"
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.SDMolSupplier",
            return_value=[mock_mol],
        )
        with pytest.raises(RuntimeError, match="multi-component"):
            _pubchem_2d_to_3d("SOME_SDF_DATA")

    def test_embed_succeeds_on_first_attempt(self, mocker):
        mock_mol = mocker.Mock()
        mock_mol.GetPropsAsDict.return_value = {"PUBCHEM_COMPONENT_COUNT": 1}
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.SDMolSupplier",
            return_value=[mock_mol],
        )
        mock_mol2 = mocker.Mock()
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.AddHs",
            return_value=mock_mol2,
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.EmbedMolecule",
            return_value=0,
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.MolToMolBlock",
            return_value="MOLBLOCK_RESULT",
        )
        result = _pubchem_2d_to_3d("SOME_SDF_DATA")
        assert result == "MOLBLOCK_RESULT"

    def test_embed_succeeds_on_second_attempt_with_random_coords(self, mocker):
        mock_mol = mocker.Mock()
        mock_mol.GetPropsAsDict.return_value = {"PUBCHEM_COMPONENT_COUNT": 1}
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.SDMolSupplier",
            return_value=[mock_mol],
        )
        mock_mol2 = mocker.Mock()
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.AddHs",
            return_value=mock_mol2,
        )
        mock_embed = mocker.patch(
            "chemsmart.io.molecules.pubchem.EmbedMolecule",
            side_effect=[-1, 0],
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.MolToMolBlock",
            return_value="MOLBLOCK_RESULT",
        )
        result = _pubchem_2d_to_3d("SOME_SDF_DATA")
        assert result == "MOLBLOCK_RESULT"
        assert mock_embed.call_count == 2

    def test_embed_fails_both_attempts_raises_runtime_error(self, mocker):
        mock_mol = mocker.Mock()
        mock_mol.GetPropsAsDict.return_value = {"PUBCHEM_COMPONENT_COUNT": 1}
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.SDMolSupplier",
            return_value=[mock_mol],
        )
        mock_mol2 = mocker.Mock()
        mock_mol2.GetProp.return_value = "failed_mol"
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.AddHs",
            return_value=mock_mol2,
        )
        mocker.patch(
            "chemsmart.io.molecules.pubchem.EmbedMolecule",
            side_effect=[-1, -1],
        )
        with pytest.raises(RuntimeError, match="Could not embed"):
            _pubchem_2d_to_3d("SOME_SDF_DATA")

    def test_wrong_number_of_molecules_raises_assertion_error(self, mocker):
        mocker.patch(
            "chemsmart.io.molecules.pubchem.Chem.SDMolSupplier",
            return_value=[],
        )
        with pytest.raises(AssertionError, match="exactly one molecule"):
            _pubchem_2d_to_3d("SOME_SDF_DATA")
