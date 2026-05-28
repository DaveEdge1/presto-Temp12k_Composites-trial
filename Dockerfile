# Temperature 12k multi-method composite container (R + Python).
#
# Base = davidedge/lipd_webapps:lipdbase2, which already ships the R stack the
# composite methods need (lipdR, geoChronR, compositeR) inside an renv project
# rooted at /. We add a small Python venv for the LiPD-pickle adapter, the GAM
# method (pygam), and the consensus / NetCDF / figure steps.
#
# Layers ordered least->most frequently changed so script edits rebuild fast.

# ── Layer 1: base (R + lipdR/geoChronR/compositeR via renv at /) ──────────
FROM davidedge/lipd_webapps:lipdbase2

ENV DEBIAN_FRONTEND=noninteractive TZ=UTC

# ── Layer 2: Python venv (heavy, cached until deps change) ────────────────
# The base image's python3 ships without pip/ensurepip, so install
# python3-venv + python3-pip via apt first, then build an isolated venv.
# numpy/pandas/scipy/xarray/netcdf4/matplotlib + pygam (GAM) + pyyaml.
# netcdf4 manylinux wheels bundle libnetcdf, so no apt system lib needed.
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3-venv python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --no-cache-dir --upgrade pip setuptools wheel && \
    /opt/venv/bin/pip install --no-cache-dir \
        numpy pandas scipy pyyaml xarray netcdf4 "matplotlib>=3.6" pygam && \
    /opt/venv/bin/python -c "import numpy, pandas, scipy, yaml, xarray, netCDF4, matplotlib, pygam; print('venv deps OK')"

ENV PYTHON_BIN=/opt/venv/bin/python

# ── Layer 3: app source (cheap; iterates often) ───────────────────────────
WORKDIR /app
COPY scripts/        /app/scripts/
COPY config/         /app/config/
COPY reference_data/ /app/reference_data/
COPY entrypoint.sh   /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh

# CI mounts at runtime:
#   /proxies/lipd_legacy.pkl     (RO)  proxy data from lipdverse
#   /app/config/user_config.yml  (RO)  overwritten per run by PReSto
#   /results                     (RW)  outputs land here
ENV LIPD_PICKLE=/proxies/lipd_legacy.pkl \
    PRESTO_CONFIG=/app/config/user_config.yml \
    PRESTO_REFDATA=/app/reference_data \
    PRESTO_OUTPUT=/results

ENTRYPOINT ["/app/entrypoint.sh"]
