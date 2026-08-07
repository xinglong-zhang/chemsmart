from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.settings.crest import CRESTProjectSettings


class TestCRESTJobSettings:
    def test_merge_dict(self):
        settings = CRESTJobSettings.default()
        settings.charge = None
        settings.multiplicity = None

        merged_settings = settings.merge({"charge": -2, "multiplicity": 3})
        assert merged_settings.charge == -2
        assert merged_settings.multiplicity == 3

    def test_get_settings_from_yaml(
        self, crest_yaml_settings_gas_project_name
    ):
        project_settings = CRESTProjectSettings.from_project(
            crest_yaml_settings_gas_project_name
        )
        conformers = project_settings.conformer_settings()

        assert isinstance(conformers, CRESTJobSettings)
        assert conformers.jobtype == "conformers"
        assert conformers.gfn_version == "gfn2"
        assert conformers.energy_window == 6.0
        assert conformers.optimization_level == "vtight"
        assert conformers.nci is False
        assert conformers.solvent_model is None
        assert conformers.solvent_id is None

    def test_get_settings_from_yaml_solv(
        self, crest_yaml_settings_solv_project_name
    ):
        project_settings = CRESTProjectSettings.from_project(
            crest_yaml_settings_solv_project_name
        )
        conformers = project_settings.conformer_settings()

        assert isinstance(conformers, CRESTJobSettings)
        assert conformers.jobtype == "conformers"
        assert conformers.gfn_version == "gfnff"
        assert conformers.optimization_level == "tight"
        assert conformers.solvent_model == "gbsa"
        assert conformers.solvent_id == "water"
