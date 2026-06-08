from chemsmart.agent.wizard.scheduler_env import (
    _PBS_BIN_CANDIDATES,
    _SGE_ROOT_CANDIDATES,
    _SLURM_BIN_CANDIDATES,
)


def test_sge_candidate_paths_include_bright() -> None:
    assert "/cm/local/apps/uge/current" in _SGE_ROOT_CANDIDATES
    assert "/cm/local/apps/sge/current" in _SGE_ROOT_CANDIDATES
    assert "/cm/shared/apps/uge/current" in _SGE_ROOT_CANDIDATES
    assert "/cm/shared/apps/sge/current" in _SGE_ROOT_CANDIDATES
    assert "/opt/gridengine" in _SGE_ROOT_CANDIDATES


def test_slurm_candidate_paths_include_ohpc() -> None:
    assert "/cm/local/apps/slurm/current/bin" in _SLURM_BIN_CANDIDATES
    assert "/opt/slurm/current/bin" in _SLURM_BIN_CANDIDATES
    assert "/opt/ohpc/pub/apps/slurm/current/bin" in _SLURM_BIN_CANDIDATES


def test_pbs_candidate_paths_include_pbspro_default() -> None:
    assert "/opt/pbs/default/bin" in _PBS_BIN_CANDIDATES
    assert "/var/spool/pbs/bin" in _PBS_BIN_CANDIDATES
    assert "/cm/local/apps/pbspro/current/bin" in _PBS_BIN_CANDIDATES
