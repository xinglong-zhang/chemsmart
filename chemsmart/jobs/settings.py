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
        semiempirical=None,
        defgrid=None,
        charge=None,
        multiplicity=None,
        freq=True,
        numfreq=False,
        job_type=None,
        title=None,
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
        input_string=None,
        **kwargs,
    ):
        self.ab_initio = ab_initio
        self.functional = functional
        self.dispersion = dispersion
        self.basis = basis
        self.semiempirical = semiempirical
        self.defgrid = defgrid
        self.charge = charge
        self.multiplicity = multiplicity
        self.freq = freq
        self.numfreq = numfreq
        self.job_type = job_type
        self.title = title
        self.solvent_model = solvent_model
        self.solvent_id = solvent_id
        self.additional_route_parameters = additional_route_parameters
        self.route_to_be_written = route_to_be_written
        self.modred = modred
        self.gen_genecp_file = gen_genecp_file
        self.heavy_elements = heavy_elements
        self.heavy_elements_basis = heavy_elements_basis
        self.light_elements_basis = light_elements_basis

        if custom_solvent is not None:
            if not isinstance(custom_solvent, str):
                raise ValueError(
                    "Custom solvent parameters must be a string! It can be either a string"
                    "giving the path of the custom solvent file or the custom solvent parameters"
                    "in free string format."
                )
            if os.path.exists(custom_solvent) and os.path.isfile(
                custom_solvent
            ):
                self.set_custom_solvent_via_file(custom_solvent)
            else:
                self.custom_solvent = custom_solvent
            # check that the last line of custom_solvent is empty, if not, add an empty line
            if self.custom_solvent[-1] != "\n":
                self.custom_solvent += "\n"
            logger.debug(f"Custom solvent parameters: {self.custom_solvent}")
        else:
            self.custom_solvent = None

        self.forces = forces
        self.input_string = input_string

    def remove_solvent(self):
        self.solvent_model = None
        self.solvent_id = None
        self.custom_solvent = None

    def update_solvent(self, solvent_model=None, solvent_id=None):
        """Update solvent model and solvent identity for implicit solvation.

        Solvent models available: ['pcm', 'iefpcm', 'cpcm', 'smd', 'dipole', 'ipcm', 'scipcm'].
        """
        # update only if not None; do not update to default value of None
        if solvent_model is not None:
            self._check_solvent(solvent_model)
            self.solvent_model = solvent_model

        if solvent_id is not None:
            self.solvent_id = solvent_id

    def _check_solvent(self, solvent_model):
        pass

    def modify_solvent(self, remove_solvent=False, **kwargs):
        if not remove_solvent:
            self.update_solvent(**kwargs)
        else:
            self.remove_solvent()

    def set_custom_solvent_via_file(self, filename):
        if not os.path.exists(os.path.expanduser(filename)):
            raise ValueError(f"File {filename} does not exist!")

        # path to the file for custom_solvent parameters
        with open(filename) as f:
            lines = f.readlines()

        lines = [line.strip() for line in lines]

        self.custom_solvent = "\n".join(lines)

    @classmethod
    def from_dict(cls, settings_dict):
        return cls(**settings_dict)


def read_molecular_job_yaml(filename, program="gaussian"):
    # read in defaults, if exists
    file_directory = os.path.dirname(filename)
    default_file = os.path.join(file_directory, "defaults.yaml")
    default_config = {}
    if os.path.exists(default_file):
        with open(default_file) as f:
            default_config = yaml.safe_load(f)
        logger.debug(
            f"Using the following pre-set defaults: \n{default_config}"
        )
    else:
        logger.warning("Default file settings does not exist.\n")
        if program == "gaussian":
            from chemsmart.settings.gaussian import GaussianJobSettings

            default_config = GaussianJobSettings.default().__dict__
        elif program == "orca":
            from chemsmart.settings.orca import ORCAJobSettings

            default_config = ORCAJobSettings.default().__dict__
        else:
            # other programs may be implemented in future
            pass
        logger.debug(
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
        "crest",
        "dias",
        "resp",
        "set",
        "traj",
        "uvvis",
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
            # turn off freq calculation for single point calculations
            all_project_configs[job]["freq"] = False
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
