#!/bin/sh
docplex config --upgrade /opt/ibm/cplex/ibm/ILOG/CPLEX_Studio2212
PYTHONPATH=/app exec python /app/migemox/pipeline/main.py "$@"