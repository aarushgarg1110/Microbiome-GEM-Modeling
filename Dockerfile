# Stage 1: uv binary
FROM ghcr.io/astral-sh/uv:0.7.3 AS uv

# Stage 2: Python 3.11 builder
FROM python:3.11-slim AS builder
ENV UV_COMPILE_BYTECODE=1 UV_PYTHON_DOWNLOADS=0 UV_LINK_MODE=copy

WORKDIR /app

# Set environment for CPLEX
ENV CPLEX_HOME=/cplex/ibm/cplex
ENV PATH="${CPLEX_HOME}/bin/x86-64_linux:$PATH"

# Copy installed CPLEX files
COPY cplex /opt/ibm/cplex

# Install everything using uv
RUN --mount=from=uv,source=/uv,target=/bin/uv \
    --mount=type=cache,target=/root/.cache/uv \
    --mount=type=bind,source=uv.lock,target=uv.lock \
    --mount=type=bind,source=pyproject.toml,target=pyproject.toml \
    uv sync --locked --no-install-project --no-dev

# Stage 3: Runtime
FROM python:3.11-slim

# Minimal runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libglib2.0-0 libgomp1 libquadmath0 libexpat1 && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /app /app
COPY --from=builder /opt/ibm/cplex /opt/ibm/cplex
COPY --from=builder /app/.venv /app/.venv
COPY migemox /app/migemox
COPY test_data_input /test_data_input

ENV PATH="/app/.venv/bin:$PATH"

COPY entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh
ENTRYPOINT ["/app/entrypoint.sh"]