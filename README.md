# Microbiome-GEM-Modeling

Microbiome-GEM-Modeling is a computational framework for constructing, analyzing, and simulating genome-scale metabolic models (GEMs) of microbial communities. This project enables researchers to integrate metagenomic data, reconstruct metabolic networks, and perform in silico experiments to investigate microbiome function and interactions.

## Features

- Automated reconstruction of GEMs from annotated genomes
- Integration of metagenomic abundance data
- Community-level metabolic modeling and simulation
- Support for constraint-based modeling with CPLEX
- Modular pipeline for reproducible analyses

## Getting Started

To get started, clone this repository and follow the instructions in [Summary.md](./Summary.md) for detailed workflow guidance.

The modeling pipeline is containerized for reproducibility. See the `Dockerfile` for environment details.

## Citing & Credits

This project builds on the work of the COBRApy community and leverages IBM CPLEX for optimization. Please cite relevant tools and datasets as appropriate for your research.

- COBRApy: [https://github.com/opencobra/cobrapy](https://github.com/opencobra/cobrapy)
- IBM CPLEX: [https://www.ibm.com/products/ilog-cplex-optimization-studio](https://www.ibm.com/products/ilog-cplex-optimization-studio)

Project lead: [Aarush Garg, Zomorrodi Lab]

For questions or contributions, please open an issue or pull request.

## How to Run This Project with Docker

This project requires IBM CPLEX with an Academic License, which cannot be redistributed directly. You must manually extract and prepare the CPLEX files before building or running the modeling scripts.

### 1. Install Docker

If you do not have Docker installed, download and install it from [https://www.docker.com/get-started/](https://www.docker.com/get-started/). Follow the platform-specific instructions for your operating system.

### 2. Prepare CPLEX Files

1. **Download the CPLEX installer** from IBM SkillsBuild or the IBM Academic Initiative and place it in a directory named `cplex` inside your project root.

2. **Extract CPLEX using a Debian container:**
    ```bash
    docker run -it --rm -v "${PWD}:/mnt" debian:bullseye bash
    apt update && apt install -y openjdk-11-jre-headless
    cd /mnt/cplex/
    ./cplex_studio2212.linux_x86_64.bin
    # Accept all license agreements and follow prompts
    cd ../opt
    tar czf /mnt/cplex.tar.gz ibm
    exit
    ```

3. **Unpack the CPLEX files on your host:**
    ```bash
    tar xzf cplex.tar.gz
    # You should now have Microbiome-GEM-Modeling/cplex/ibm/ILOG/CPLEX_Studio2212
    ```

### 3. Build and Run the Docker Image

1. **Build the Docker image:**
    ```bash
    docker build -t cobra .
    ```

2. **Run the modeling pipeline:**
    ```bash
    docker run --rm -v "${PWD}:/mnt" cobra
    ```

   - To run interactively (for debugging or development):
     ```bash
     docker run -it --rm -v "${PWD}:/mnt" cobra bash
     docplex config --upgrade ./cplex/ibm/ILOG/CPLEX_Studio2212
     python src/pipeline/pipeline.py
     ```

For further details on workflow and usage, see [Summary.md](./Summary.md).



