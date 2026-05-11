# chemnode1 OpenPBS provisioning runbook

> Scope: single-node Rocky Linux 8.10 GCE instance in `asia-northeast1-a`
> matching the audited `chemnode1` runtime.
>
> Canonical bootstrap path: `startup-script` metadata that points at
> `scripts/ops/provision-chemnode-pbs.sh`.
>
> Audit date: 2026-05-12 (KST)

## Why this runbook exists

`chemnode1` is not reproducible from the repository state that existed before
this audit:

- instance metadata still declares a cloud-init `user-data` payload for SGE
- the live node is actually Rocky 8.10 + OpenHPC OpenPBS 23.06.06
- `cloud-init` is not installed on the audited image, so `user-data` is not a
  valid bootstrap channel
- `google-startup-scripts.service` reported `No startup scripts to run`
- `dnf history` shows the live node was rebuilt manually with
  `ohpc-release-2-1.el8` and `openpbs-server-ohpc`

This runbook makes the provisioning path explicit and re-runnable.

## Drift summary: declared vs actual

| Fact | Declared / leftover artifact | Audited runtime |
| --- | --- | --- |
| Bootstrap channel | Instance metadata `user-data` cloud-config | `cloud-init` absent; use `startup-script` metadata instead |
| Scheduler | Grid Engine / SGE packages + `inst_sge` | `openpbs-server-ohpc-23.06.06-270.ohpc.3.2.x86_64` |
| Repo source | none documented | `ohpc-release-2-1.el8.x86_64.rpm` installed before OpenPBS |
| Services | `sgemaster.default`, `sge_execd.default` in metadata | `pbs.service` enabled and active |
| Runtime config path | `/tmp/sge_autoconfig.conf`, `/usr/local/sbin/firstboot.sh` | neither file exists on the live node |
| Scratch / runtime | implied only | `/home/chemsmart/scratch`, `/opt/miniforge3/envs/chemsmart` |

## Prerequisites

1. A GCP project with billing enabled and the Compute Engine API enabled.
2. `gcloud` authenticated with permission to create instances, firewall rules,
   and read serial logs.
3. An SSH public key available either through project metadata or OS Login.
4. An ingress firewall rule that allows TCP/22 from the operator IP.
5. This repository checked out locally so
   `scripts/ops/provision-chemnode-pbs.sh` can be passed with
   `--metadata-from-file startup-script=...`.
6. Read the companion wizard contract in fork PR
   `Hongjiseung-ROK/chemsmart#90` (merged as `5a2ffd05`).

## Expected hourly cost

Approximate on-demand base cost in `asia-northeast1` (Tokyo):

- `e2-medium`: about `$0.03350571 / hour`
- 30 GiB standard persistent disk: about
  `30 * $0.000054795 = $0.00164385 / hour`
- expected baseline total: about `$0.03515 / hour`

This excludes network egress, snapshots, premium images, and any future price
changes.

## Canonical provisioning flow

### 1. Set environment variables

```bash
export PROJECT_ID='REPLACE_ME'
export ZONE='asia-northeast1-a'
export INSTANCE_NAME='chemnode1'
export MACHINE_TYPE='e2-medium'
export IMAGE_PROJECT='rocky-linux-cloud'
export IMAGE_FAMILY='rocky-linux-8-optimized-gcp'
export BOOT_DISK_SIZE_GB='30'
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
  startup-script=scripts/ops/provision-chemnode-pbs.sh
```

Do **not** attach the script as `user-data`. The audited Rocky image does not
ship `cloud-init`, so `user-data` is silently ignored.

### 4. Watch first boot

```bash
gcloud compute instances tail-serial-port-output "$INSTANCE_NAME" \
  --project "$PROJECT_ID" \
  --zone "$ZONE"
```

The script also logs to `/var/log/chemsmart-chemnode-pbs-provision.log`.

## Post-provision health checks

After SSH access works, verify the runtime as root or via `sudo`:

```bash
sudo systemctl is-active pbs
sudo /opt/pbs/bin/qmgr -c 'list server'
sudo /opt/pbs/bin/qmgr -c 'list queue workq'
sudo /opt/pbs/bin/pbsnodes -av chemnode1
sudo test -d /home/chemsmart/scratch
sudo test -x /opt/miniforge3/bin/conda
```

Expected results:

- `pbs.service` is `active`
- default queue is `workq`
- node `chemnode1` reports `state = free`
- `chemsmart` scratch exists at `/home/chemsmart/scratch`
- Miniforge base exists at `/opt/miniforge3`

## Wizard contract

Mirror the live validation contract from fork PR
`Hongjiseung-ROK/chemsmart#90` (merge commit
`5a2ffd05cb045538fac103a1a039a2cde4f19b85`). The validated excerpt is:

```yaml
SERVER:
  SCHEDULER: PBS
  QUEUE_NAME: workq
  MEM_GB: 3
  NUM_CORES: 2
  SUBMIT_COMMAND: qsub
  HOST: chemsmart@34.146.86.156
GAUSSIAN:
  EXEFOLDER: null
  CONDA_ENV: |
    source /opt/miniforge3/etc/profile.d/conda.sh
    conda activate chemsmart
ORCA:
  EXEFOLDER: null
  CONDA_ENV: |
    source /opt/miniforge3/etc/profile.d/conda.sh
    conda activate chemsmart
NCIPLOT:
  EXEFOLDER: null
  CONDA_ENV: |
    source /opt/miniforge3/etc/profile.d/conda.sh
    conda activate chemsmart
```

Operational notes that sit beside that YAML contract:

- scratch should exist at `/home/chemsmart/scratch`
- the runtime scheduler package is OpenPBS from OpenHPC, not SGE
- the provisioning script must leave `qsub`, `qmgr`, and `pbsnodes` usable
  without manual repair

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
- whether to harden SELinux / firewall policy beyond the currently audited
  single-node lab baseline
