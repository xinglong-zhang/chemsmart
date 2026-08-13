# syntax=docker/dockerfile:1.7

FROM mambaorg/micromamba:2.8.0@sha256:de51d1d07fb8bed487d5d04ab2034ae05a0f0972bb05fab64c6c923a0c1c888b AS builder

# This throw-away stage installs a prefix; the final runtime is non-root.
# hadolint ignore=DL3002,DL3066
USER root
ENV PATH=/opt/conda/bin:${PATH} \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1

COPY containers/cpu/conda-linux-64.lock /tmp/conda-linux-64.lock
RUN micromamba install --yes --name base \
      --file /tmp/conda-linux-64.lock \
    && micromamba clean --all --yes

WORKDIR /src
COPY containers/cpu/build-requirements.lock /tmp/build-requirements.lock
COPY containers/cpu/requirements-linux-amd64.lock /tmp/requirements-linux-amd64.lock
COPY pyproject.toml README.md LICENSE LICENSE.LESSER MANIFEST.in ./
COPY chemsmart ./chemsmart

RUN python -m pip install \
      --no-cache-dir \
      --require-hashes \
      --no-build-isolation \
      --requirement /tmp/build-requirements.lock \
    && python -m pip wheel \
      --no-cache-dir \
      --no-deps \
      --no-build-isolation \
      --wheel-dir /tmp/wheels \
      . \
    && python -m pip install \
      --no-cache-dir \
      --require-hashes \
      --no-build-isolation \
      --only-binary=pyscf \
      --requirement /tmp/requirements-linux-amd64.lock \
    && python -m pip install \
      --no-cache-dir \
      --no-deps \
      /tmp/wheels/chemsmart-3.1.4-py3-none-any.whl \
    && python -m pip check \
    && find /opt/conda -type d -name __pycache__ -prune -exec rm -rf '{}' + \
    && rm -rf /root/.cache /tmp/wheels /tmp/*.lock

FROM ubuntu:22.04@sha256:3b06811b2afd352be909dd088a004166d665dc76d38b13eada33522a9d915c6f AS runtime

ARG CHEMSMART_SOURCE_REVISION=unknown
ARG CONTAINER_RECIPE_REVISION=unknown

LABEL org.opencontainers.image.title="ChemSmart CPU Research Runtime" \
      org.opencontainers.image.description="ChemSmart 3.1.4 with pinned PySCF and xTB CPU backends" \
      org.opencontainers.image.source="https://github.com/Hongjiseung-ROK/chemsmart" \
      org.opencontainers.image.version="3.1.4-cpu" \
      org.opencontainers.image.revision="${CONTAINER_RECIPE_REVISION}" \
      io.chemsmart.source-revision="${CHEMSMART_SOURCE_REVISION}" \
      io.chemsmart.runtime-transaction-regression="safe-preview-preflight-materialization-v1"

ENV PATH=/opt/conda/bin:${PATH} \
    HOME=/home/chemsmart \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    TZ=UTC \
    CHEMSMART_AGENT_SERVER=local \
    PYTHONNOUSERSITE=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    OMP_STACKSIZE=4G \
    OMP_MAX_ACTIVE_LEVELS=1 \
    XDG_CACHE_HOME=/tmp/.cache \
    MPLCONFIGDIR=/tmp/matplotlib \
    TMPDIR=/tmp

# Scientific packages are fully locked. Ubuntu runtime libraries come from
# the digest-pinned base's supported apt repository so security fixes remain
# available without baking brittle archive versions into this recipe.
# hadolint ignore=DL3008
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
      ca-certificates \
      libgl1 \
      libglu1-mesa \
      libgomp1 \
      libice6 \
      libsm6 \
      libxext6 \
      libxrender1 \
      tini \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --gid 1000 chemsmart \
    && useradd --uid 1000 --gid 1000 --create-home --shell /bin/bash chemsmart \
    && install -d -o chemsmart -g chemsmart -m 0755 \
      /workspace /scratch /opt/chemsmart/qualification \
      /home/chemsmart/.chemsmart/server /home/chemsmart/.chemsmart/agent \
    && install -d -m 0755 /usr/share/licenses/chemsmart

COPY --from=builder /opt/conda /opt/conda
COPY --chown=chemsmart:chemsmart containers/cpu/config/server/local.yaml /home/chemsmart/.chemsmart/server/local.yaml
COPY --chown=chemsmart:chemsmart containers/cpu/config/usersettings.yaml /home/chemsmart/.chemsmart/usersettings.yaml
COPY --chown=chemsmart:chemsmart containers/cpu/config/agent/agent.yaml /home/chemsmart/.chemsmart/agent/agent.yaml
COPY --chown=chemsmart:chemsmart containers/cpu/qualification /opt/chemsmart/qualification
COPY --chown=chemsmart:chemsmart containers/cpu/entrypoint.sh /usr/local/bin/chemsmart-container-entrypoint
COPY LICENSE LICENSE.LESSER /usr/share/licenses/chemsmart/

RUN chmod 0755 /usr/local/bin/chemsmart-container-entrypoint \
    && chown -R chemsmart:chemsmart /opt/chemsmart /home/chemsmart \
      /workspace /scratch

USER 1000:1000
WORKDIR /workspace

ENTRYPOINT ["/usr/bin/tini", "--", "/usr/local/bin/chemsmart-container-entrypoint"]
CMD ["chemsmart", "--help"]
