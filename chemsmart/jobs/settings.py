import logging
import os
import yaml

from chemsmart.utils.utils import update_dict_with_existing_keys

logger = logging.getLogger(__name__)


class MolecularJobSettings:
    """Common base job settings for molecular systems using Gaussian and ORCA jobs."""

    def __init__(
        self,
        ab_initio=None,
        functional=None,
        dispersion=None,
        basis=None,
        defgrid=None,
        charge=None,
        multiplicity=None,
        freq=True,
        numfreq=False,
        job_type=None,
        solvent_model=None,
        solvent_id=None,
        additional_route_parameters=None,
        route_to_be_written=None,
        modred=None,
        gen_genecp_file=None,
        heavy_elements=None,
        heavy_elements_basis=None,
        light_elements_basis=None,
        custom_solvent=None,
        forces=False,
        **kwargs,
    ):
        self.ab_initio = ab_initio
        self.functional = functional
        self.dispersion = dispersion
        self.basis = basis
        self.defgrid = defgrid
        self.charge = charge
        self.multiplicity = multiplicity
        self.freq = freq
        self.numfreq = numfreq
        self.job_type = job_type
        self.solvent_model = solvent_model
        self.solvent_id = solvent_id
        self.additional_route_parameters = additional_route_parameters
        self.route_to_be_written = route_to_be_written
        self.modred = modred
        self.gen_genecp_file = gen_genecp_file
        self.heavy_elements = heavy_elements
        self.heavy_elements_basis = heavy_elements_basis
        self.light_elements_basis = light_elements_basis
        self.custom_solvent = custom_solvent
        self.forces = forces


def read_molecular_job_yaml(filename):
    # read in defaults, if exists
    file_directory = os.path.dirname(filename)
    default_file = os.path.join(file_directory, "defaults.yaml")
    if os.path.exists(default_file):
        with open(default_file) as f:
            default_config = yaml.safe_load(f)
    else:
        logger.warning("Default file settings does not exist.\n")
        from pyatoms.jobs.orca.settings import ORCAJobSettings

        default_config = ORCAJobSettings.default().__dict__
        logger.warning(
            f"Using the following pre-set defaults: \n{default_config}"
        )

    # job types
    gas_phase_jobs = [
        "opt",
        "modred",
        "ts",
        "irc",
        "scan",
        "nci",
        "crestopt",
        "dias",
        "resp",
        "saopt",
        "wbi",
    ]
    sp_job = ["sp"]
    td_job = ["td"]
    all_jobs = gas_phase_jobs + sp_job + td_job

    # read in project config
    with open(filename) as f:
        project_config = yaml.safe_load(f)

    # populate job settings for different jobs
    all_project_configs = {}  # store all job settings in a dict

    # check if solv settings exist
    solv_config = project_config.get("solv", None)

    # check if separate gas phase settings exist
    gas_config = project_config.get("gas", None)

    if gas_config is None:
        # no settings for gas phase; using implicit solvation model for all jobs
        for job in all_jobs:
            all_project_configs[job] = (
                default_config.copy()
            )  # populate defaults
            all_project_configs[job]["job_type"] = job  # update job_type
            all_project_configs[job] = update_dict_with_existing_keys(
                all_project_configs[job], solv_config
            )
    else:
        # settings for gas phase exist - also solv settings exist
        for job in gas_phase_jobs:  # jobs using gas config
            all_project_configs[job] = (
                default_config.copy()
            )  # populate defaults
            all_project_configs[job]["job_type"] = job  # update job_type
            all_project_configs[job] = update_dict_with_existing_keys(
                all_project_configs[job], gas_config
            )
        for job in sp_job:  # jobs using solv config
            all_project_configs[job] = (
                default_config.copy()
            )  # populate defaults
            all_project_configs[job]["job_type"] = job  # update job_type
            all_project_configs[job] = update_dict_with_existing_keys(
                all_project_configs[job], solv_config
            )

    # check if td settings exist (optional)
    if "td" in project_config:
        td_config = project_config["td"]
        for job in td_job:  # jobs using td config s
            all_project_configs[job] = (
                default_config.copy()
            )  # populate defaults
            all_project_configs[job]["job_type"] = job  # update job_type
            all_project_configs[job] = update_dict_with_existing_keys(
                all_project_configs[job], td_config
            )

    return all_project_configs
