# ComPy: A Python Pipeline for Personalized Microbiome Community Modeling

## Introduction

This document provides a detailed description of the ComPy pipeline, a Python-based workflow for creating and simulating personalized microbiome community metabolic models. The pipeline is an implementation of the core functionalities of the `mgPipe` workflow from the MATLAB COBRA Toolbox. The conceptual framework is based on the methods described in "The Microbiome Modeling Toolbox" by Heinken et al. (*Bioinformatics*, 2019), while the specific formulation of the coupling constraints is adapted from the methods detailed in Heinken et al. (*Gut Microbes*, 2013).

The primary goal of ComPy is to transform species abundance data from metagenomic samples into individualized, predictive metabolic models. These models can then be used to simulate the metabolic output of a microbial community under specific dietary conditions, providing insights into host-microbiome interactions in health and disease.

## Pipeline Workflow and Structure

The pipeline is divided into two main stages, each executed by a series of core functions. The structure ensures a logical flow from raw data to interpretable results, mirroring the original `mgPipe` process.

#### # Stage 1: Personalized Community Model Construction
*   **Top-level script:** `compy.py`
*   **Core Functions:**
    *   `compy()`: Main function that orchestrates the entire model building process.
    *   `build_global_model()`: Loads all individual species models and merges them into a single, unified community model, applying necessary modifications.
    *   `build_sample_model()`: Takes the global model and customizes it for a specific sample by pruning absent species and adding an abundance-weighted community biomass reaction.

#### # Stage 2: Diet-Constrained Simulation and Analysis
*   **Top-level script:** `diet_adaptation.py`
*   **Core Functions:**
    *   `simulate_microbiota_models()`: The main driver that manages the simulation for all samples in parallel.
    *   `adapt_vmh_diet_to_agora()`: Pre-processes a standard diet file to make it compatible with the community models, adding essential metabolites and adjusting bounds.
    *   `process_single_sample()`: The workhorse function that applies the diet to a single sample model, runs optimizations, and performs Flux Variability Analysis (FVA).

## Stage 1: Personalized Community Model Construction

This stage, orchestrated by the `compy.py` script, constructs a unique community model for each sample by integrating individual species models according to their abundance.

### # Step 1.1: From Single Species to Community Member (`species_to_community`)

The process begins with individual species models (e.g., AGORA `.mat` files). The `species_to_community` function, called within `build_global_model`, transforms each model to allow its integration into a larger community structure.

*   **Reaction & Metabolite Tagging:** To prevent overlap and track species-specific contributions, all intracellular reactions and metabolites are tagged with a species-specific prefix. For example, an intracellular metabolite `ade[c]` in the *B. theta* model becomes `B_theta_ade[c]`, and a reaction `PGI` becomes `B_theta_PGI`.

*   **Compartment Modification:** The extracellular space is redefined to enable inter-species communication. The original extracellular compartment `[e]` is converted into a species-specific lumen compartment `[u]`. For example, a metabolite `ade[e]` that could be exchanged with the environment now becomes `B_theta_ade[u]`, a metabolite specific to that species but residing in the shared lumen.

*   **New Reactions - Inter-Species Exchange (IEX):** To allow species to share and compete for metabolites, the `_create_inter_species_exchange` helper function establishes a common lumen compartment `[u]`. Species-specific metabolites are connected to this shared pool via new IEX reactions of the format `ade[u] <=> B_theta_ade[u]`. These are named e.g., `B_theta_IEX_ade[u]tr` and are fully reversible, allowing a species to either secrete a metabolite into the common lumen or take it up.

### # Step 1.2: Integrating the Host-Microbiome Interface (`add_diet_fecal_compartments`)

After all species models for a community are assembled, `add_diet_fecal_compartments` (called by `build_global_model`) adds compartments representing the interface with the host.

*   **Removed Reactions:** All original exchange reactions from the base models (e.g., `EX_ade(e)`) are removed. This is a critical step to ensure all nutrient uptake and waste secretion occurs through the new, controlled interface.

*   **New Compartments & Reactions:** The function introduces the diet `[d]` and fecal `[fe]` compartments. It then creates a chain of new transport reactions for every metabolite in the common lumen, connecting the diet to the fecal output:
    1.  **Diet Exchange:** `EX_metabolite[d]: metabolite[d] <=> `
    2.  **Diet-to-Lumen Transport (DUt):** `DUt_metabolite: metabolite[d] -> metabolite[u]`
    3.  **Lumen-to-Fecal Transport (UFEt):** `UFEt_metabolite: metabolite[u] -> metabolite[fe]`
    4.  **Fecal Exchange:** `EX_metabolite[fe]: metabolite[fe] <=> `

### # Step 1.3: Formulating Biomass Coupling Constraints (`build_global_coupling_constraints`)

A standard constraint-based model relies on the stoichiometric matrix **S** to enforce mass balance (`S * v = 0`). However, this alone is insufficient for community models. To address this, `build_global_coupling_constraints` introduces additional linear constraints that link the flux of every reaction within a species directly to that species' growth rate.

*   **Mathematical Formulation:** The formulation is based on the principles outlined in Heinken et al. (2013). For each reaction `i` belonging to a specific species, its flux (`v_i`) is coupled to that species' biomass reaction flux (`v_biomass`):
    *   **For irreversible reactions** (`v_i >= 0`): `v_i ≤ c × v_biomass`
    *   **For reversible reactions:** `v_i ≤ c × v_biomass` (forward) and `v_i ≥ -c × v_biomass` (reverse)

*   **Implementation:** These inequalities are rewritten into the standard linear constraint format `A*x <= b`. For example, `v_i - c × v_biomass ≤ 0`. Here, `c` is the **coupling factor**, set to `400` as in the original `mgPipe` implementation. This ensures that if `v_biomass` is zero, `v_i` must also be zero, making the model more biologically realistic. These rules are compiled into a sparse **coupling matrix (C)** and associated vectors.

### # Step 1.4: Building the Sample-Specific Model (`build_sample_model`)

Finally, `build_sample_model` customizes the global model for each individual sample.

*   **Pruning:** The `prune_zero_abundance_species` helper is used to remove all reactions and metabolites associated with species whose abundance is below a defined threshold (`1e-7`). The corresponding rows in the coupling matrix are also removed via `prune_coupling_constraints_by_species`.
*   **Community Biomass Reaction:** The `com_biomass` function creates a single, unified objective function for the community. It consumes the individual biomass metabolites of each present species, weighted by that species' relative abundance, according to the stoichiometry: `Σ (abundance_i * species_i_biomass[c]) -> microbeBiomass[u]`.
*   **Objective Function:** The model's objective is set to maximize the excretion of this community biomass product: `model.objective = "EX_microbeBiomass[fe]"`.

The final sample-specific model, complete with all its constraints, is saved as a `.mat` file using `make_mg_pipe_model_dict`.

## Stage 2: Diet-Constrained Simulation and Analysis

This stage, driven by `diet_adaptation.py`, applies dietary constraints and simulates the metabolic potential of the constructed models.

### # Step 2.1: Dietary Adjustment (`adapt_vmh_diet_to_agora`)

The `adapt_vmh_diet_to_agora` function pre-processes a standard diet file to make it compatible with the community models. It renames diet reaction IDs (e.g., `EX_glc(e)` becomes `Diet_EX_glc_D[d]`) and ensures feasibility by adding essential metabolites (`ESSENTIAL_METS`, `UNMAPPED_METS`) and relaxing bounds for certain `MICRONUTRIENTS`.

### # Step 2.2: The Full FBA Problem and its Components

The `process_single_sample` function sets up and solves a constrained Flux Balance Analysis (FBA) problem. The helper function `build_constraint_matrix` in `diet_adaptation.py` illustrates how these components are assembled from the model file.

*   **Objective Function:** `Maximize c' * v`, where `c` is a vector with `1` at the position of the objective reaction (`EX_microbeBiomass[fe]`) and `0` elsewhere.
*   **Constraints:**
    1.  **Mass Balance:** `S * v = b`, where `S` is the stoichiometric matrix.
    2.  **Flux Bounds:** `lb <= v <= ub`, where the lower bounds of diet reactions are set by the diet file.
    3.  **Coupling Constraints:** `C * v <= d`, where `C` is the coupling matrix built in Stage 1.

### # Step 2.3: Flux Variability Analysis (FVA) - A Deep Dive

The goal of FVA is to determine the metabolic flexibility of the community by calculating the minimum and maximum possible flux for each reaction while the system remains in a near-optimal state. The implementation of this step is crucial and is handled by `_analyze_metabolite_fluxes`. This pipeline contains two parallel implementations for this task, revealing both the underlying mechanics and the practical, high-performance approach.

#### # Part A: The Manual, "From-Scratch" FVA Implementation

This approach, preserved in commented-out code, demonstrates how the FVA problem is constructed from its fundamental components.

1.  **Assembling the Full Problem (`build_constraint_matrix`):** This function loads all mathematical components from the `.mat` file: the stoichiometric matrix (`S`), the coupling matrix (`C`), their respective right-hand-side vectors (`b` and `d`), constraint sense vectors (`csense`, `dsense`), and reaction bounds (`lb`, `ub`). It combines them into a single, large constraint system (`A*v <= rhs`).
2.  **Creating a Solver-Ready Model (`build_optlang_model`):** This function translates the raw matrices into a structured optimization model that a solver like CPLEX can understand. It explicitly defines each flux `v_i` as a variable and builds each constraint equation from the rows of the `A` matrix.
3.  **Executing Sequential FVA (`run_sequential_fva`):** This function performs the core FVA logic:
    *   First, it solves the primary FBA problem to find the maximum community biomass, `Z_opt`.
    *   It then adds a new constraint to the model: `v_biomass >= 0.9999 * Z_opt`. This forces the model to only consider solutions that achieve near-optimal growth.
    *   Finally, it iterates through every exchange reaction of interest. For each one, it temporarily sets it as the model's objective and solves twice: once to minimize its flux and once to maximize it. This process reveals the full range of metabolic activity possible *while the community is thriving*.

#### # Part B: The Abstracted, High-Performance FVA Implementation

This is the active method used by the pipeline for efficiency. It leverages the highly optimized `cobrapy` library to achieve the same result.

1.  **Applying Coupling Constraints (`couple_constraints`):** This function is the critical bridge. It reads the coupling matrices (`C`, `d`, `dsense`) from the `.mat` file and, using the `optlang` API, adds them directly to the `cobra.Model` object's internal solver interface. Once this is done, the COBRA model itself is inherently aware of the coupling rules for all subsequent calculations.
2.  **Using `cobra.flux_variability_analysis`:** With the model now fully constrained, the pipeline can call the built-in `flux_variability_analysis` function. This function executes the same logic as the manual method (constraining the objective and iterating) but does so far more efficiently, using optimized code and parallel processing capabilities.

Both methods solve the same biological problem, but the second approach is used in practice for its significant performance advantage.

### # Step 2.4: Reporting Final Fluxes

The FVA results are used by `_analyze_metabolite_fluxes` to calculate the net metabolic capabilities for each metabolite. For each metabolite `met`, the following are calculated from the FVA minima and maxima:

*   **Net Production:** `NetProd[met] = abs(min_flux_diet[Diet_EX_met[d]] + max_flux_fecal[EX_met[fe]])`
*   **Net Uptake:** `NetUptake[met] = abs(max_flux_diet[Diet_EX_met[d]] + min_flux_fecal[EX_met[fe]])`

Here, `max_flux_fecal` is the maximum secretion into feces, and `min_flux_diet` is the maximum uptake from the diet (a large negative value). The Net Production value represents the community's *de novo* synthesis capability for that metabolite. The final values are collected across all samples and written to `inputDiet_net_secretion_fluxes.csv` and `inputDiet_net_uptake_fluxes.csv`.

---
## References

1.  Heinken, A., Baldini, F., Heirendt, L., Magnusdottir, S., Fleming, R. M. T., & Thiele, I. (2019). The Microbiome Modeling Toolbox: from microbial interactions to personalized microbial communities. *Bioinformatics*, 35(13), 2332-2334.
2.  Heinken, A., Sahoo, S., Fleming, R. M. T., & Thiele, I. (2013). Systems-level characterization of a host-microbe metabolic symbiosis in the mammalian gut. *Gut Microbes*, 4(1), 28-40.