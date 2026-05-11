#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_NAME
readonly LOG_FILE=${CHEMSMART_LOG_FILE:-/var/log/chemsmart-chemnode-pbs-provision.log}
readonly CHEMSMART_HOSTNAME=${CHEMSMART_HOSTNAME:-chemnode1}
readonly CHEMSMART_USER=${CHEMSMART_USER:-chemsmart}
readonly CHEMSMART_ENV_NAME=${CHEMSMART_ENV_NAME:-chemsmart}
readonly MINIFORGE_PREFIX=${MINIFORGE_PREFIX:-/opt/miniforge3}
readonly SCRATCH_DIR=${SCRATCH_DIR:-/home/${CHEMSMART_USER}/scratch}
readonly WORK_QUEUE=${WORK_QUEUE:-workq}
readonly OHPC_RELEASE_RPM=${OHPC_RELEASE_RPM:-https://repos.openhpc.community/OpenHPC/2/EL_8/x86_64/ohpc-release-2-1.el8.x86_64.rpm}

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

ensure_dnf_package() {
  local package=$1

  if ! rpm -q "$package" >/dev/null 2>&1; then
    dnf install -y "$package"
  fi
}

ensure_ohpc_repo() {
  if ! rpm -q ohpc-release >/dev/null 2>&1; then
    dnf install -y "$OHPC_RELEASE_RPM"
  fi
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
    TMP_MINIFORGE=$(mktemp -p /var/tmp miniforge-pbs-XXXXXX.sh)
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

set_config_value() {
  local file=$1
  local key=$2
  local value=$3

  if grep -q "^${key}=" "$file"; then
    sed -i "s|^${key}=.*$|${key}=${value}|" "$file"
  else
    printf '%s=%s\n' "$key" "$value" >> "$file"
  fi
}

configure_pbs() {
  local pbs_path

  pbs_path=/opt/pbs/bin:/opt/pbs/sbin:$PATH

  set_config_value /etc/pbs.conf PBS_SERVER "$CHEMSMART_HOSTNAME"
  set_config_value /etc/pbs.conf PBS_START_SERVER 1
  set_config_value /etc/pbs.conf PBS_START_SCHED 1
  set_config_value /etc/pbs.conf PBS_START_COMM 1
  set_config_value /etc/pbs.conf PBS_START_MOM 1

  systemctl enable pbs
  systemctl restart pbs

  PATH=$pbs_path qmgr -c "list queue ${WORK_QUEUE}" >/dev/null 2>&1 || \
    PATH=$pbs_path qmgr -c "create queue ${WORK_QUEUE}"
  PATH=$pbs_path qmgr -c \
    "set queue ${WORK_QUEUE} queue_type = Execution"
  PATH=$pbs_path qmgr -c "set queue ${WORK_QUEUE} enabled = True"
  PATH=$pbs_path qmgr -c "set queue ${WORK_QUEUE} started = True"
  PATH=$pbs_path qmgr -c "set server default_queue = ${WORK_QUEUE}"

  if ! PATH=$pbs_path pbsnodes -av "$CHEMSMART_HOSTNAME" >/dev/null 2>&1; then
    PATH=$pbs_path qmgr -c "create node ${CHEMSMART_HOSTNAME}"
  fi
}

write_completion_marker() {
  printf 'FIRSTBOOT_DONE\n' > /var/log/chemsmart-firstboot.done
}

verify_runtime() {
  local pbs_path

  pbs_path=/opt/pbs/bin:/opt/pbs/sbin:$PATH

  systemctl is-active --quiet pbs
  PATH=$pbs_path qmgr -c 'list server' >/dev/null
  PATH=$pbs_path qmgr -c "list queue ${WORK_QUEUE}" | grep -q \
    "Queue ${WORK_QUEUE}"
  PATH=$pbs_path pbsnodes -av "$CHEMSMART_HOSTNAME" | grep -q 'state = '
  test -d "$SCRATCH_DIR"
  test -x "$MINIFORGE_PREFIX/bin/conda"
}

main() {
  require_root
  setup_logging

  log "starting ${SCRIPT_NAME}"
  assert_bootstrap_channel
  ensure_hostname
  ensure_dnf_package epel-release
  ensure_dnf_package curl
  ensure_ohpc_repo
  ensure_dnf_package openpbs-server-ohpc
  ensure_chemsmart_user
  ensure_miniforge
  configure_pbs
  verify_runtime
  write_completion_marker
  log 'OpenPBS provisioning completed successfully'
}

main "$@"
