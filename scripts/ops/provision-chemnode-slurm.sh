#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_NAME
readonly LOG_FILE=${CHEMSMART_LOG_FILE:-/var/log/chemsmart-chemnode-slurm-provision.log}
readonly CHEMSMART_HOSTNAME=${CHEMSMART_HOSTNAME:-chemnode2}
readonly CHEMSMART_USER=${CHEMSMART_USER:-chemsmart}
readonly CHEMSMART_ENV_NAME=${CHEMSMART_ENV_NAME:-chemsmart}
readonly MINIFORGE_PREFIX=${MINIFORGE_PREFIX:-/opt/miniforge3}
readonly SCRATCH_DIR=${SCRATCH_DIR:-/home/${CHEMSMART_USER}/scratch}
readonly SLURM_PARTITION_NAME=${SLURM_PARTITION_NAME:-workq}
readonly SLURM_REAL_MEMORY_MB=${SLURM_REAL_MEMORY_MB:-3000}
readonly SLURM_CPUS=${SLURM_CPUS:-$(nproc)}
readonly SELINUX_MODE=${SELINUX_MODE:-permissive}
readonly DISABLE_FIREWALL=${DISABLE_FIREWALL:-true}

TMP_MINIFORGE=""

log() {
  printf '[%s] %s\n' "$(date -Iseconds)" "$*"
}

cleanup() {
  if [[ -n "$TMP_MINIFORGE" && -f "$TMP_MINIFORGE" ]]; then
    rm -f "$TMP_MINIFORGE"
  fi
}

on_error() {
  local exit_code=$?
  local line_no=$1

  log "ERROR: ${SCRIPT_NAME} failed at line ${line_no} with exit ${exit_code}"
  cleanup
  exit "$exit_code"
}

trap 'on_error $LINENO' ERR
trap cleanup EXIT

require_root() {
  if [[ $EUID -ne 0 ]]; then
    printf '%s\n' "${SCRIPT_NAME} must run as root" >&2
    exit 1
  fi
}

setup_logging() {
  mkdir -p "$(dirname "$LOG_FILE")"
  touch "$LOG_FILE"
  chmod 0644 "$LOG_FILE"
  exec > >(tee -a "$LOG_FILE") 2>&1
}

select_miniforge_url() {
  case "$(uname -m)" in
    aarch64 | arm64)
      printf '%s\n' \
        'https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh'
      ;;
    *)
      printf '%s\n' \
        'https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh'
      ;;
  esac
}

assert_bootstrap_channel() {
  if command -v cloud-init >/dev/null 2>&1; then
    log 'cloud-init detected; canonical path still uses GCE startup-script metadata'
    return
  fi

  if systemctl list-unit-files google-startup-scripts.service \
    >/dev/null 2>&1; then
    log 'cloud-init absent; use --metadata-from-file startup-script=...'
    return
  fi

  log 'cloud-init absent and google-startup-scripts.service missing; rerun manually as root'
}

ensure_hostname() {
  local current_host
  local primary_ip

  current_host=$(hostnamectl --static)
  if [[ "$current_host" != "$CHEMSMART_HOSTNAME" ]]; then
    hostnamectl set-hostname "$CHEMSMART_HOSTNAME"
  fi

  if ! grep -Eq "(^|[[:space:]])${CHEMSMART_HOSTNAME}($|[[:space:]])" \
    /etc/hosts; then
    primary_ip=$(hostname -I | awk '{print $1}')
    printf '%s %s\n' "$primary_ip" "$CHEMSMART_HOSTNAME" >> /etc/hosts
  fi
}

ensure_system_account() {
  local user=$1
  local group=$2
  local shell=$3

  if ! getent group "$group" >/dev/null 2>&1; then
    groupadd --system "$group"
  fi

  if ! id "$user" >/dev/null 2>&1; then
    useradd --system -g "$group" -M -r -s "$shell" "$user"
  fi
}

ensure_dnf_package() {
  local package=$1

  if ! rpm -q "$package" >/dev/null 2>&1; then
    dnf install -y "$package"
  fi
}

configure_selinux_and_firewall() {
  if [[ "$SELINUX_MODE" == permissive ]]; then
    setenforce 0 || true
    if grep -q '^SELINUX=' /etc/selinux/config; then
      sed -i 's/^SELINUX=.*/SELINUX=permissive/' /etc/selinux/config
    fi
  fi

  if [[ "$DISABLE_FIREWALL" == true ]]; then
    systemctl disable --now firewalld >/dev/null 2>&1 || true
  fi
}

write_slurm_conf() {
  local tmp_conf

  tmp_conf=$(mktemp -p /var/tmp slurm-conf-XXXXXX)
  cat > "$tmp_conf" <<EOF
ClusterName=${CHEMSMART_HOSTNAME}
SlurmctldHost=${CHEMSMART_HOSTNAME}
AuthType=auth/munge
CryptoType=crypto/munge
MpiDefault=none
ProctrackType=proctrack/linuxproc
ReturnToService=2
SlurmctldPidFile=/var/run/slurmctld.pid
SlurmctldPort=6817
SlurmdPidFile=/var/run/slurmd.pid
SlurmdPort=6818
SlurmdSpoolDir=/var/spool/slurmd
SlurmUser=slurm
StateSaveLocation=/var/spool/slurmctld
SwitchType=switch/none
TaskPlugin=task/none
InactiveLimit=0
KillWait=30
MinJobAge=300
SlurmctldTimeout=120
SlurmdTimeout=300
Waittime=0
SchedulerType=sched/backfill
SelectType=select/cons_tres
SelectTypeParameters=CR_Core_Memory
AccountingStorageType=accounting_storage/none
JobAcctGatherType=jobacct_gather/none
SlurmctldDebug=info
SlurmctldLogFile=/var/log/slurmctld.log
SlurmdDebug=info
SlurmdLogFile=/var/log/slurmd.log
NodeName=${CHEMSMART_HOSTNAME} CPUs=${SLURM_CPUS} RealMemory=${SLURM_REAL_MEMORY_MB} State=UNKNOWN
PartitionName=${SLURM_PARTITION_NAME} Nodes=${CHEMSMART_HOSTNAME} Default=YES MaxTime=24:00:00 State=UP
EOF

  if [[ ! -f /etc/slurm/slurm.conf ]] || ! cmp -s "$tmp_conf" /etc/slurm/slurm.conf; then
    install -m 0644 -o slurm -g slurm "$tmp_conf" /etc/slurm/slurm.conf
  fi

  rm -f "$tmp_conf"
}

ensure_munge_key() {
  if [[ ! -f /etc/munge/munge.key ]]; then
    /usr/sbin/create-munge-key -f
  fi

  chown munge:munge /etc/munge/munge.key
  chmod 0400 /etc/munge/munge.key
}

ensure_chemsmart_user() {
  if ! id "$CHEMSMART_USER" >/dev/null 2>&1; then
    useradd -m -s /bin/bash "$CHEMSMART_USER"
  fi

  install -d -m 0755 -o "$CHEMSMART_USER" -g "$CHEMSMART_USER" \
    "$SCRATCH_DIR"
}

ensure_miniforge() {
  local miniforge_url

  miniforge_url=$(select_miniforge_url)
  if [[ ! -x "$MINIFORGE_PREFIX/bin/conda" ]]; then
    TMP_MINIFORGE=$(mktemp -p /var/tmp miniforge-slurm-XXXXXX.sh)
    curl -fsSL "$miniforge_url" -o "$TMP_MINIFORGE"
    bash "$TMP_MINIFORGE" -b -p "$MINIFORGE_PREFIX"
    "$MINIFORGE_PREFIX/bin/conda" config --system --set auto_update_conda false
  fi

  # shellcheck disable=SC1091
  source "$MINIFORGE_PREFIX/etc/profile.d/conda.sh"

  if ! conda env list | awk 'NR > 2 {print $1}' | grep -qx \
    "$CHEMSMART_ENV_NAME"; then
    conda create -y -n "$CHEMSMART_ENV_NAME"
  fi

  chown -R root:root "$MINIFORGE_PREFIX"
}

write_completion_marker() {
  printf 'FIRSTBOOT_DONE\n' > /var/log/chemsmart-firstboot.done
}

verify_runtime() {
  systemctl is-active --quiet munge
  systemctl is-active --quiet slurmctld
  systemctl is-active --quiet slurmd
  sinfo -Nel | grep -q "$CHEMSMART_HOSTNAME"
  scontrol show node "$CHEMSMART_HOSTNAME" | grep -q 'State='
  test -d /var/spool/slurmctld
  test -d /var/spool/slurmd
  test -x "$MINIFORGE_PREFIX/bin/conda"
  test -d "$SCRATCH_DIR"
}

main() {
  require_root
  setup_logging

  log "starting ${SCRIPT_NAME}"
  assert_bootstrap_channel
  ensure_hostname
  configure_selinux_and_firewall
  ensure_system_account slurm slurm /sbin/nologin
  ensure_dnf_package epel-release
  ensure_dnf_package dnf-plugins-core
  dnf config-manager --set-enabled crb
  ensure_dnf_package curl
  ensure_dnf_package munge
  ensure_dnf_package slurm
  ensure_dnf_package slurm-slurmctld
  ensure_dnf_package slurm-slurmd
  ensure_munge_key
  install -d -m 0755 -o slurm -g slurm /var/spool/slurmctld
  install -d -m 0755 -o slurm -g slurm /var/spool/slurmd
  write_slurm_conf
  ensure_chemsmart_user
  ensure_miniforge
  systemctl enable munge slurmctld slurmd
  systemctl restart munge
  systemctl restart slurmctld
  systemctl restart slurmd
  verify_runtime
  write_completion_marker
  log 'Slurm provisioning completed successfully'
}

main "$@"
