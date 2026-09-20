"""
ORCA reaction workflow: optional endpoint opt, NEB-TS guess, TS opt and SP.

Guess reuses ``ORCANEBJob`` with project ``neb_settings()`` (default
``NEB-TS``). Child Opt/TS/SP jobs are stock classes sequenced with
``ChainMixin``.
"""

import os

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.neb import ORCANEBJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import (
    ORCAJobSettings,
    ORCANEBJobSettings,
    ORCATSJobSettings,
)
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.orca.ts import ORCATSJob
from chemsmart.jobs.reaction import ReactionChainMixin

_TS_ONLY_KEYS = (
    "inhess",
    "inhess_filename",
    "hybrid_hess",
    "hybrid_hess_atoms",
    "numhess",
    "recalc_hess",
    "trust_radius",
    "tssearch_type",
    "scants_modred",
    "full_scan",
)
_NEB_ONLY_KEYS = (
    "joboption",
    "nimages",
    "ending_xyzfile",
    "intermediate_xyzfile",
    "restarting_xyzfile",
    "preopt_ends",
)


def _settings_as(settings, cls, drop_keys):
    """Copy ``settings`` as ``cls``, dropping attributes from other jobs."""
    if isinstance(settings, cls):
        return settings.copy()
    data = dict(settings.copy().__dict__)
    for key in drop_keys:
        data.pop(key, None)
    return cls(**data)


class ORCAReactionJob(ReactionChainMixin, ORCAJob):
    """ORCA R/TS/P chain: optional NEB-TS guess, then opt+freq and SP."""

    TYPE = "orcareaction"
    uses_neb = True
    _opt_job_class = ORCAOptJob
    _ts_job_class = ORCATSJob
    _sp_job_class = ORCASinglePointJob

    def __init__(self, molecule, settings=None, **kwargs):
        if not isinstance(settings, ORCAJobSettings):
            raise ValueError(
                f"Settings must be instance of ORCAJobSettings for "
                f"{self.__class__.__name__}, but is {settings} instead!"
            )
        super().__init__(molecule=molecule, settings=settings, **kwargs)

    def _ts_settings_for(self, molecule):
        settings = super()._ts_settings_for(molecule)
        if isinstance(settings, ORCATSJobSettings):
            return settings
        converted = _settings_as(settings, ORCATSJobSettings, _NEB_ONLY_KEYS)
        converted.jobtype = "ts"
        converted.freq = True
        return self._apply_charge_mult(converted, molecule)

    def _neb_settings_for(self, molecule):
        if self.neb_settings is not None:
            base = self.neb_settings
        elif self.opt_settings is not None:
            base = self.opt_settings
        else:
            base = self.settings
        settings = _settings_as(base, ORCANEBJobSettings, _TS_ONLY_KEYS)
        settings.jobtype = "neb"
        settings.freq = False
        if settings.joboption is None:
            settings.joboption = "NEB-TS"
        return self._apply_charge_mult(settings, molecule)

    def _make_guess_job(self, reactant, product, ts_guess=None):
        settings = self._neb_settings_for(reactant)
        settings.preopt_ends = False
        settings.ending_xyzfile = os.path.join(
            self.folder, f"{self.label}_neb_P.xyz"
        )
        if ts_guess is not None:
            settings.intermediate_xyzfile = os.path.join(
                self.folder, f"{self.label}_neb_TS.xyz"
            )
        return ORCANEBJob(
            molecule=reactant,
            settings=settings,
            label=f"{self.label}_neb",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _prepare_guess_inputs(self):
        if self.guess_job is None:
            return
        settings = self.guess_job.settings
        product_path = settings.ending_xyzfile
        if product_path:
            self._optimized_product_for_guess().write_xyz(
                product_path, mode="w"
            )
        ts_path = settings.intermediate_xyzfile
        if ts_path and self.path_ts_guess is not None:
            self.path_ts_guess.write_xyz(ts_path, mode="w")
