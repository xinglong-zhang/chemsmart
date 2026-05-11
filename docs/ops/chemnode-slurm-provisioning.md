# chemnode2 Slurm provisioning runbook

> Scope: single-node Rocky Linux 9.7 GCE instance in `asia-northeast1-a`
> matching the audited `chemnode2` runtime.
>
> Canonical bootstrap path: `startup-script` metadata that points at
> `scripts/ops/provision-chemnode-slurm.sh`.
>
> Audit date: 2026-05-12 (KST)

## Why this runbook exists

`chemnode2` was created with a cloud-init `user-data` payload, but the audited
runtime shows:

- `cloud-init` is not installed on the Rocky 9 image
- `google-startup-scripts.service` logged `No startup scripts to run`
- the live node was recovered manually by running `/tmp/chemnode2-setup.sh`
- the recovery script tried `chown slurm:slurm` before a guaranteed `slurm`
  user existed, so spool ownership needed manual remediation later

This runbook replaces that ad-hoc recovery path with a single idempotent
startup script.

## Drift summary: declared vs actual

| Fact | Declared / leftover artifact | Audited runtime |
| --- | --- | --- |
| Bootstrap channel | Instance metadata `user-data` cloud-config | `cloud-init` absent; use `startup-script` metadata instead |
| Scheduler | SLURM in `user-data`, but no canonical repo script | `slurm-22.05.9`, `slurmctld`, `slurmd`, `munge` active |
| Recovery path | `/usr/local/sbin/firstboot.sh` in metadata | file absent; manual `/tmp/chemnode2-setup.sh` performed the install |
| Ordering | `chown slurm:slurm` before a guaranteed system user | live node repaired with `useradd --system slurm` + spool recreation |
| Runtime config | implicit only | `/etc/slurm/slurm.conf`, `/var/spool/slurmctld`, `/var/spool/slurmd` |
| Host hardening | permissive / disabled only in temp script | SELinux is `Permissive`, `firewalld` is disabled |

## Prerequisites

1. A GCP project with billing enabled and the Compute Engine API enabled.
2. `gcloud` authenticated with permission to create instances, firewall rules,
   and read serial logs.
3. An SSH public key available either through project metadata or OS Login.
4. An ingress firewall rule that allows TCP/22 from the operator IP.
5. This repository checked out locally so
   `scripts/ops/provision-chemnode-slurm.sh` can be passed with
   `--metadata-from-file startup-script=...`.
6. Read the companion wizard contract in fork PR
   `Hongjiseung-ROK/chemsmart#91` (current validation head
   `dbafc5e1b03cc3a067601b23701bd0676288e039`).

## Expected hourly cost

Approximate on-demand base cost in `asia-northeast1` (Tokyo):

- `e2-medium`: about `$0.03350571 / hour`
- 20 GiB standard persistent disk: about
  `20 * $0.000054795 = $0.00109590 / hour`
- expected baseline total: about `$0.03460 / hour`

This excludes network egress, snapshots, premium images, and any future price
changes.

## Canonical provisioning flow

### 1. Set environment variables

```bash
export PROJECT_ID='REPLACE_ME'
export ZONE='asia-northeast1-a'
export INSTANCE_NAME='chemnode2'
export MACHINE_TYPE='e2-medium'
export IMAGE_PROJECT='rocky-linux-cloud'
export IMAGE_FAMILY='rocky-linux-9-optimized-gcp'
export BOOT_DISK_SIZE_GB='20'
export OPERATOR_CIDR='203.0.113.10/32'
```

### 2. Create the SSH firewall rule once per project

```bash
gcloud compute firewall-rules create chemsmart-ssh \
  --project "$PROJECT_ID" \
  --network default \
  --direction INGRESS \
  --action ALLOW \
  --rules tcp:22 \
  --source-ranges "$OPERATOR_CIDR" \
  --target-tags chemsmart
```

If the rule already exists, reuse it instead of creating a duplicate.

### 3. Create the instance with `startup-script` metadata

```bash
gcloud compute instances create "$INSTANCE_NAME" \
  --project "$PROJECT_ID" \
  --zone "$ZONE" \
  --machine-type "$MACHINE_TYPE" \
  --image-project "$IMAGE_PROJECT" \
  --image-family "$IMAGE_FAMILY" \
  --boot-disk-size "${BOOT_DISK_SIZE_GB}GB" \
  --tags chemsmart \
  --metadata-from-file \
  startup-script=scripts/ops/provision-chemnode-slurm.sh
```

Do **not** attach the script as `user-data`. The audited Rocky image does not
ship `cloud-init`, so `user-data` is silently ignored.

### 4. Watch first boot

```bash
gcloud compute instances tail-serial-port-output "$INSTANCE_NAME" \
  --project "$PROJECT_ID" \
  --zone "$ZONE"
```

The script also logs to `/var/log/chemsmart-chemnode-slurm-provision.log`.

## Post-provision health checks

After SSH access works, verify the runtime as root or via `sudo`:

```bash
sudo systemctl is-active munge
sudo systemctl is-active slurmctld
sudo systemctl is-active slurmd
sudo sinfo -Nel
sudo scontrol show node chemnode2
sudo test -d /var/spool/slurmctld
sudo test -d /var/spool/slurmd
sudo test -d /home/chemsmart/scratch
sudo test -x /opt/miniforge3/bin/conda
```

Expected results:

- `munge`, `slurmctld`, and `slurmd` are all `active`
- partition `workq` exists and the node is `idle`
- both Slurm spool directories are owned by `slurm:slurm`
- Miniforge base exists at `/opt/miniforge3`
- scratch exists at `/home/chemsmart/scratch`

## Wizard contract

Mirror the live validation contract from fork PR
`Hongjiseung-ROK/chemsmart#91` (validation head
`dbafc5e1b03cc3a067601b23701bd0676288e039`). The validated excerpt is:

```yaml
SERVER:
  SCHEDULER: SLURM
  QUEUE_NAME: workq
  NUM_HOURS: 24
  MEM_GB: 3
  NUM_CORES: 2
  SUBMIT_COMMAND: sbatch
  HOST: chemsmart@35.189.152.165
  SCRATCH_DIR: ~/scratch
GAUSSIAN:
  EXEFOLDER: null
  CONDA_ENV: 'source /opt/miniforge3/etc/profile.d/conda.sh

    conda activate chemsmart'
ORCA:
  EXEFOLDER: null
  CONDA_ENV: 'source /opt/miniforge3/etc/profile.d/conda.sh

    conda activate chemsmart'
NCIPLOT:
  EXEFOLDER: null
  CONDA_ENV: 'source /opt/miniforge3/etc/profile.d/conda.sh

    conda activate chemsmart'
```

Operational notes that sit beside that YAML contract:

- local runtime scratch should exist at `/home/chemsmart/scratch`
- `sinfo --json` is not required for correctness on this host; the wizard fix
  in PR #91 falls back to `scontrol show node`
- the provisioning script must guarantee `slurm` user creation before any
  ownership-sensitive directory work

## Teardown

Delete the instance when the audit environment is no longer needed:

```bash
gcloud compute instances delete "$INSTANCE_NAME" \
  --project "$PROJECT_ID" \
  --zone "$ZONE"
```

If the firewall rule was created only for this cluster and is no longer needed:

```bash
gcloud compute firewall-rules delete chemsmart-ssh \
  --project "$PROJECT_ID"
```

## Remaining manual decisions

This runbook intentionally leaves these as operator-supplied inputs:

- project ID and billing setup
- SSH key / OS Login posture
- operator ingress CIDR
- whether to keep using the default VPC or a dedicated network
- whether to keep the currently audited SELinux-permissive /
  firewalld-disabled lab baseline, or replace it with site-specific hardening
