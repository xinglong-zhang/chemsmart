"""
Gaussian reaction workflow: optional endpoint opt, QST2/QST3 guess, TS opt and SP.

QST input reuses ``GaussianInputWriter`` and ``settings.route_string``. The
only route change is the optimization keyword (``opt=qst2`` or ``opt=qst3``).
The assembled ``.com`` is stored on ``settings.input_string`` and run as
``GaussianComJob``. Child Opt/TS/SP jobs are stock classes sequenced with
``ChainMixin``.
"""

import logging
import re
from io import StringIO

from chemsmart.jobs.gaussian.job import GaussianComJob, GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.ts import GaussianTSJob
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.jobs.reaction import (
    ReactionChainMixin,
    validate_path_search_structures,
)
from chemsmart.utils.repattern import (
    gaussian_opt_keywords_pattern,
    multiple_spaces_pattern,
)

logger = logging.getLogger(__name__)

_QST2 = "qst2"
_QST3 = "qst3"
_OPT_JOB_TOKENS = ("ts", "qst2", "qst3")


def make_qst_com_job(
    reactant,
    product,
    settings,
    ts_guess=None,
    label=None,
    jobrunner=None,
    **kwargs,
):
    """Create a ``GaussianComJob`` whose ``.com`` is the QST input string."""
    if not isinstance(settings, GaussianJobSettings):
        raise ValueError(
            "settings must be a GaussianJobSettings instance, "
            f"got {type(settings).__name__}."
        )
    validate_path_search_structures(reactant, product, ts_guess=ts_guess)

    settings = settings.copy()
    qst_kind = _QST3 if ts_guess is not None else _QST2
    if label is None:
        label = reactant.get_chemical_formula(empirical=True)
    settings.route_to_be_written = _route_with_qst_opt(settings, qst_kind)

    job = GaussianComJob(
        molecule=reactant,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
    writer = _QSTGaussianInputWriter(
        job=job, product=product, ts_guess=ts_guess
    )
    buf = StringIO()
    writer._write_all(buf)
    job.settings.input_string = buf.getvalue()
    logger.debug(
        "Built %s input for GaussianComJob label=%s",
        qst_kind,
        label,
    )
    return job


def _route_with_qst_opt(settings, qst_kind):
    """Reuse ``settings.route_string``, replacing only the opt keyword."""
    if qst_kind not in (_QST2, _QST3):
        raise ValueError(
            f"qst_kind must be '{_QST2}' or '{_QST3}', got {qst_kind!r}."
        )
    route = settings.copy().route_string

    def _replace_opt(match):
        return _qst_opt_keyword(match.group(0), qst_kind) + " "

    new_route, n_replaced = re.subn(
        gaussian_opt_keywords_pattern,
        _replace_opt,
        route,
        count=1,
        flags=re.IGNORECASE,
    )
    if n_replaced == 0:
        new_route = re.sub(
            r"^(#\S*)",
            rf"\1 opt={qst_kind}",
            route,
            count=1,
        )
    return re.sub(multiple_spaces_pattern, " ", new_route).strip()


def _qst_opt_keyword(opt_token, qst_kind):
    token = opt_token.strip()
    eq_index = token.find("=")
    if eq_index == -1:
        return f"opt={qst_kind}"
    value = token[eq_index + 1 :].strip()
    if value.startswith("(") and value.endswith(")"):
        parts = [
            part.strip() for part in value[1:-1].split(",") if part.strip()
        ]
        replaced = False
        new_parts = []
        for part in parts:
            if part.lower() in _OPT_JOB_TOKENS:
                new_parts.append(qst_kind)
                replaced = True
            else:
                new_parts.append(part)
        if not replaced:
            new_parts.insert(0, qst_kind)
        return f"opt=({','.join(new_parts)})"
    return f"opt={qst_kind}"


class _QSTGaussianInputWriter(GaussianInputWriter):
    """Write two or three QST geometry blocks with the standard writer."""

    def __init__(self, job, product, ts_guess=None):
        super().__init__(job=job)
        self.product = product
        self.ts_guess = ts_guess

    def _geometry_blocks(self):
        title = self.settings.title if self.settings.title else "Gaussian job"
        blocks = [
            (self.job.molecule, f"{title}: reactant"),
            (self.product, f"{title}: product"),
        ]
        if self.ts_guess is not None:
            blocks.append((self.ts_guess, f"{title}: ts guess"))
        return blocks

    def _write_all(self, f):
        self._write_gaussian_header(f)
        self._write_route_section(f)

        molecule = self.job.molecule
        title = self.settings.title
        charge = self.settings.charge
        multiplicity = self.settings.multiplicity
        try:
            for mol, block_title in self._geometry_blocks():
                self.job.molecule = mol
                self.settings.title = block_title
                self.settings.charge = (
                    mol.charge if mol.charge is not None else charge
                )
                self.settings.multiplicity = (
                    mol.multiplicity
                    if mol.multiplicity is not None
                    else multiplicity
                )
                self._write_gaussian_title(f)
                self._write_charge_and_multiplicity(f)
                self._write_cartesian_coordinates(f)
        finally:
            self.job.molecule = molecule
            self.settings.title = title
            self.settings.charge = charge
            self.settings.multiplicity = multiplicity

        self._append_modredundant(f)
        self._append_gen_genecp_basis(f)
        self._append_custom_solvent_parameters(f)
        self._append_job_specific_info(f)
        self._append_other_additional_info(f)


class GaussianReactionJob(ReactionChainMixin, GaussianJob):
    """Gaussian R/TS/P chain: optional QST guess, then opt+freq and SP."""

    TYPE = "g16reaction"
    _opt_job_class = GaussianOptJob
    _ts_job_class = GaussianTSJob
    _sp_job_class = GaussianSinglePointJob

    def __init__(self, molecule, settings=None, **kwargs):
        if not isinstance(settings, GaussianJobSettings):
            raise ValueError(
                f"Settings must be instance of GaussianJobSettings for "
                f"{self.__class__.__name__}, but is {settings} instead!"
            )
        super().__init__(molecule=molecule, settings=settings, **kwargs)

    def _make_guess_job(self, reactant, product, ts_guess=None):
        settings = self._ts_settings_for(reactant)
        settings.freq = False
        return make_qst_com_job(
            reactant=reactant,
            product=product,
            settings=settings,
            ts_guess=ts_guess,
            label=f"{self.label}_qst",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
