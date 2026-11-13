"""Direct unit tests for chemsmart.io.gaussian.GaussianRefs/BSEMetadata."""

import builtins

import pytest

from chemsmart.io.gaussian import BSEMetadata, GaussianRefs


class TestGaussianRefs:
    def test_reference_lists_are_accessible_via_properties(self):
        refs = GaussianRefs()
        assert refs.gaussian_ab_initio == refs.g_ab_initio
        assert refs.gaussian_dft_fuctionals == refs.g_functionals
        assert refs.gaussian_basis_sets == refs.g_bases
        assert refs.gaussian_semiempirical_methods == refs.g_semiempirical
        assert refs.gaussian_solvation_models == refs.g_solvation_models
        assert (
            refs.gaussian_additional_route_parameters
            == refs.g_additional_route_parameters
        )
        assert (
            refs.gaussian_additional_opt_options
            == refs.g_additional_opt_options
        )
        assert refs.gaussian_dieze_tags == refs.g_dieze_tags


class TestBSEMetadata:
    def test_raises_importerror_when_basis_set_exchange_missing(
        self, monkeypatch
    ):
        real_import = builtins.__import__

        def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == "basis_set_exchange":
                raise ImportError("no module named basis_set_exchange")
            return real_import(name, globals, locals, fromlist, level)

        monkeypatch.setattr(builtins, "__import__", fake_import)

        with pytest.raises(ImportError, match="basis_set_exchange module"):
            BSEMetadata()

    def test_all_bases_names_are_lowercased(self):
        meta = BSEMetadata()
        names = meta.all_bases_names()
        assert len(names) > 0
        assert all(name == name.lower() for name in names)

    def test_all_formats_returns_dict_from_bse(self):
        meta = BSEMetadata()
        formats = meta.all_formats()
        assert isinstance(formats, dict)
        assert "gaussian94" in formats

    def test_list_of_available_formats_returns_format_keys(self):
        meta = BSEMetadata()
        formats = meta.list_of_available_formats()
        assert isinstance(formats, list)
        assert "gaussian94" in formats
        assert set(formats) == set(meta.all_formats().keys())
