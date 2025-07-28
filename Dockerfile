FROM python:3.10-slim

# Set environment for CPLEX
ENV CPLEX_HOME=/opt/ibm/cplex
ENV PATH="${CPLEX_HOME}/bin/x86-64_linux:$PATH"
ENV LD_LIBRARY_PATH="${CPLEX_HOME}/lib/x86-64_linux:$LD_LIBRARY_PATH"

# Minimal runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libglib2.0-0 libgomp1 libquadmath0 && \
    rm -rf /var/lib/apt/lists/*

# Copy installed CPLEX files (from host)
COPY cplex /opt/ibm/cplex

# Install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install CPLEX Python bindings (must match Python version)
COPY cplex/python/3.10/x86-64_linux /tmp/cplex_py/
RUN cd /tmp/cplex_py && python setup.py install && rm -rf /tmp/cplex_py

# Copy your modeling script
COPY script.py /app/script.py
WORKDIR /app

CMD ["python", "script.py"]