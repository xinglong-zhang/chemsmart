# If I comment all these out, I get the following error during run:
#   File "/Users/xinglongzhang/bin/chemsmart/chemsmart/jobs/runner.py", line 192, in from_job
#     raise ValueError(
# ValueError: Could not find any runners for job:
# GaussianOptJob<folder=/Users/xinglongzhang/pyatoms_tests, label=final_prd_opt_scan_gas_opt_opt>.
# Runners in registry: [].
#  Fake: True
#
from .crest import GaussianCrestJob
from .crestopt import GaussianCrestOptJob
from .crestts import GaussianCrestTSJob
from .custom import GaussianCustomJob
from .dias import GaussianDIASJob
from .opt import GaussianOptJob
from .irc import GaussianIRCJob
from .job import GaussianComJob, GaussianGeneralJob, GaussianJob
from .link import GaussianLinkJob
from .modred import GaussianModredJob
from .nci import GaussianNCIJob
from .resp import GaussianRESPJob
from .runner import GaussianJobRunner
from .saopt import GaussianSAOptJob
from .scan import GaussianScanJob
from .singlepoint import GaussianSinglePointJob
from .tddft import GaussianTDDFTJob
from .ts import GaussianTSJob
from .uvvis import GaussianUVVISJob
from .wbi import GaussianWBIJob

jobs = GaussianJob.subclasses()
