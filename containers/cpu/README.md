# ChemSmart CPU research image

This directory defines the private, portable `linux/amd64` ChemSmart research
runtime used for Ubuntu CPU-server qualification and benchmark work.

The image is built from ChemSmart source commit
`24e4a0e89700c88fb49bfcd673003f59a2d9cf69`. Container recipes and
documentation may advance independently, so OCI metadata records both the
ChemSmart source revision and the recipe revision.

## Included scientific runtime

| Component | Version | Role |
| --- | ---: | --- |
| ChemSmart | 3.1.4 | project YAML, CLI, agent runtime, analysis |
| Python | 3.10 | controller and PySCF compute interpreter |
| PySCF | 2.14.0 | CPU electronic structure, optimization, Hessian |
| xTB | 6.7.1 | CPU GFN calculations |
| NumPy | 1.26.4 | pinned numerical ABI |
| h5py | 3.16.0 | structured PySCF HDF5 results |
| geomeTRIC | 1.1.1 | PySCF geometry optimization |
| pyscf-dispersion | 1.5.0 | PySCF dispersion corrections |
| PyMOL open source | 3.1.0 | ChemSmart molecular runtime dependency |
| Open Babel | locked conda build | molecular format support |
| FFmpeg | 7.1.1 | visualization media dependency |

GPU4PySCF, CUDA, Gaussian, and ORCA are intentionally absent. A GPU request
must fail; it is never substituted with a CPU PySCF calculation. Gaussian and
ORCA will be integrated only in a separately qualified deployment.

## Image identity and private tags

The private package is:

```text
ghcr.io/hongjiseung-rok/chemsmart
```

- `cpu-24e4a0e` is the immutable candidate tag.
- `3.1.4-cpu` and `cpu-main` are promotion tags applied to the already-tested
  candidate digest. Promotion never rebuilds the image.

The package owner must keep the GHCR package private, remove inherited public
repository access if GitHub exposes it, and grant only the intended package
and Actions access.

Authenticate with a token that has `read:packages`, without placing it in a
shell history entry:

```bash
export GHCR_READ_TOKEN='REDACTED'
printf '%s' "$GHCR_READ_TOKEN" | \
  docker login ghcr.io --username Hongjiseung-ROK --password-stdin
unset GHCR_READ_TOKEN

docker pull ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu
```

## Runtime contract

The image runs as non-root UID/GID 1000. `/workspace`, `/scratch`, and `/tmp`
must be writable; the remaining filesystem can be read-only. xTB also needs an
unlimited process stack. A production-shaped invocation is:

```bash
mkdir -p "$PWD/workspace" "$PWD/scratch"
chmod u+rwx "$PWD/workspace" "$PWD/scratch"

docker run --rm -it \
  --read-only \
  --ulimit stack=-1:-1 \
  --cpus 8 \
  --memory 32g \
  --tmpfs /tmp:rw,nosuid,nodev,size=4g \
  --mount "type=bind,src=$PWD/workspace,dst=/workspace" \
  --mount "type=bind,src=$PWD/scratch,dst=/scratch" \
  ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu \
  chemsmart --help
```

Mounted directories must permit UID 1000 to write. On a multi-user server,
use ACLs or administrator-managed ownership instead of making them globally
writable.

The baked `local.yaml` declares two CPU cores, 8 GB RAM, zero GPUs, local
PySCF and xTB executables, and `/scratch`. CLI resource flags override these
conservative defaults. A site-specific server file may be mounted read-only:

```bash
--mount type=bind,src=/secure/config/local.yaml,\
dst=/home/chemsmart/.chemsmart/server/local.yaml,readonly
```

Do not put provider credentials in a server YAML, image layer, build argument,
OCI label, or repository file.

## ChemSmart YAML and CLI examples

Place a molecular geometry and project YAML in the mounted workspace. PySCF
single point:

```bash
docker run --rm --read-only --ulimit stack=-1:-1 \
  --cpus 4 --memory 16g --tmpfs /tmp:rw,size=2g \
  --mount "type=bind,src=$PWD/workspace,dst=/workspace" \
  --mount "type=bind,src=$PWD/scratch,dst=/scratch" \
  ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu \
  chemsmart run --server local --num-cores 4 --num-gpus 0 \
  --mem-gb 16 --no-scratch pyscf \
  --project /workspace/pyscf.yaml \
  --filename /workspace/water.xyz sp
```

xTB optimization:

```bash
docker run --rm --read-only --ulimit stack=-1:-1 \
  --cpus 4 --memory 16g --tmpfs /tmp:rw,size=2g \
  --mount "type=bind,src=$PWD/workspace,dst=/workspace" \
  --mount "type=bind,src=$PWD/scratch,dst=/scratch" \
  ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu \
  chemsmart run --server local --num-cores 4 --num-gpus 0 \
  --mem-gb 16 xtb --project /workspace/xtb.yaml \
  --filename /workspace/water.xyz opt
```

For an interactive shell, replace the final command with `/bin/bash`.

## Preview-only live agent observation

The image contains an Alibaba Token Plan profile for
`deepseek-v4-flash-0731`, maximum reasoning effort, thinking continuation, and
no fallback provider. The credential is supplied only at runtime as a
read-only secret file:

```bash
install -m 600 /dev/null /secure/provider.env
printf 'ALIBABA_TOKEN_PLAN_KEY=%s\n' "$ALIBABA_TOKEN_PLAN_KEY" \
  > /secure/provider.env

docker run --rm --read-only --ulimit stack=-1:-1 \
  --tmpfs /tmp:rw,size=2g \
  --mount "type=bind,src=$PWD/workspace,dst=/workspace" \
  --mount "type=bind,src=$PWD/scratch,dst=/scratch" \
  --mount type=bind,src=/secure/provider.env,\
dst=/run/secrets/provider.env,readonly \
  ghcr.io/hongjiseung-rok/chemsmart:3.1.4-cpu \
  chemsmart agent plan --provider alibaba-token-plan \
  --task-file /workspace/task.md \
  --secret-file /run/secrets/provider.env \
  --workspace /workspace
```

`agent plan` is preview-only. It may write project YAML and render safe
previews, but it does not execute PySCF or xTB.

## Qualification

The manual candidate workflow
`.github/workflows/cpu-container-candidate.yml` runs natively on GitHub's
Ubuntu 22.04 x86-64 runner. Before publication it:

1. checks imports, versions, the non-root/read-only boundary, and excluded
   programs;
2. previews PySCF and xTB project YAML through the live CLI;
3. runs real neutral-singlet water PySCF B3LYP/def2-SVP SP, D3BJ SP,
   geomeTRIC OPT, and optimized-geometry HESS calculations;
4. runs real GFN2-xTB SP, OPT, and optimized-geometry HESS calculations;
5. checks normal termination, state, finite energies, convergence, HDF5,
   Hessian symmetry, frequencies, and optimized-geometry handoff; and
6. runs one preview-only live agent observation.

Only `qualification-summary.json` and the public agent answer/tool trajectory
are uploaded. Raw calculations, HDF5 files, verbose logs, event stores,
credentials, and provider-private reasoning are not published.

The calculation checks qualify the installed engines. The agent answer is
read directly by a computational chemist; it is not accepted or rejected by
an answer-key, exact-DAG, token, or tool-count score. Neither form of
qualification claims that a large scientific benchmark has been reproduced.

## Rebuilding the locks

Regenerate locks only on an intentional dependency update. `uv` 0.11.23 and
`conda-lock` 3.0.4 were used for this candidate:

```bash
uv pip compile pyproject.toml containers/cpu/science-requirements.in \
  --python-version 3.10 \
  --python-platform x86_64-manylinux_2_28 \
  --generate-hashes \
  --only-binary=pyscf \
  --no-emit-package chemsmart \
  --output-file containers/cpu/requirements-linux-amd64.lock

uv pip compile containers/cpu/build-requirements.in \
  --python-version 3.10 \
  --python-platform x86_64-manylinux_2_28 \
  --generate-hashes \
  --output-file containers/cpu/build-requirements.lock

uvx --from conda-lock==3.0.4 conda-lock lock \
  --file containers/cpu/environment.yml \
  --platform linux-64 --kind lock --micromamba --without-cuda \
  --lockfile /tmp/chemsmart-cpu-conda-lock.yml
uvx --from conda-lock==3.0.4 conda-lock render \
  --kind explicit --platform linux-64 \
  --filename-template 'containers/cpu/conda-{platform}.lock' \
  /tmp/chemsmart-cpu-conda-lock.yml
```

Review the diff and rerun the full x86 qualification after any lock change.
