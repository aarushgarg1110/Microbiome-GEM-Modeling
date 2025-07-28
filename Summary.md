# MiGEMox: A Python Toolbox for Microbiome Genome-Scale Modeling

## 1. Introduction

This document provides a detailed description of the Microbiome Genome-Scale Modeling Toolbox (MiGEMox), a Python-based workflow for constructing and simulating personalized GEnome-scale Models (GEMs) of microbiome metabolism. While this initial release of the pipeline implements the core functionalities of the mgPipe workflow from "The Microbiome Modeling Toolbox in MATLAB by Heinken et al. (Bioinformatics, 2019), we plan to substantially expand extend its capabilities—beyond the original MATLAB-based toolbox—to support new modeling features, scalable analyses, and integration with modern computational technologies for microbiome systems biology.

## 2. Overview of the Pipeline Workflow and Structure

The pipeline is divided into two main stages, each executed by a series of core functions. The structure ensures a logical flow from raw data to interpretable results.

#### # Stage 1: Personalized Community Model Construction
*   **Top-level script:** `migemox.py`
*   **Core Functions:**
    *   `migemox()`: Main function that orchestrates the entire model building process.
    *   `build_global_model()`: Loads all individual organism GEMs found in the abundance input file (from model directory provided by the AGORA database) and merges them into a single, global community model subject to necessary modifications (see Section 3 for more details)
    *   `build_sample_model()`: Takes the global model and customizes it for a specific sample by pruning absent organisms and adding an abundance-weighted community biomass reaction (see Section 3 for more details).

#### # Stage 2: Diet-Constrained Simulation and Analysis
*   **Top-level script:** `diet_adaptation.py`
*   **Core Functions:**
    *   `simulate_microbiota_models()`: The main driver that manages the simulation for all samples in parallel.
    *   `adapt_vmh_diet_to_agora()`: Pre-processes a standard diet file to make it compatible with the community models, adding essential metabolites and adjusting bounds.
    *   `process_single_sample()`: The workhorse function that applies the diet to a single sample model, runs optimizations, and performs Flux Variability Analysis (FVA).

## 3. Internal Mechanics of Community Model Construction and Simulation

In this section, we provide a comprehensive overview of the workflow and implementation at each stage. Additional comments have been included in the relevant sections of the code for greater clarity. 

### Stage 1: Personalized Community Model Construction

This stage, orchestrated by the `migemox.py` script, constructs a unique community model for each sample by integrating individual organism models according to their abundance. For the remainder of this document, we use "organism" as a stand-in for any taxonomic level. One can substitute "organism" with the appropriate taxonomic category relevant to their work, such as strain or genus. 

#### # Step 1.1: From Single Species to Community Member (`species_to_community()`)

The process begins with individual GEMs for all organisms present in the abundance input file, and provided in the AGORA database. The `species_to_community()` function, called within `build_global_model()`, reformats each organism GEM by changing reaction and metabolites formatting as outlined below to allow its integration into a larger community structure.

*   **Reaction & Metabolite Tagging:** To prevent overlap and track organism-specific contributions, all intracellular reactions and metabolites are tagged with a organism-specific prefix. For example, an intracellular metabolite `ade[c]` in the *B. theta* model becomes `B_theta_ade[c]`, and a reaction `PGI` becomes `B_theta_PGI`.

*   **Compartment Modification:** The extracellular space is redefined to enable inter-species communication. The original extracellular compartment `[e]` is converted into a shared lumen compartment `[u]`. For example, a metabolite `ade[e]` that could be exchanged with the environment now becomes `B_theta_ade[u]`, a metabolite specific to that organism but residing in the shared lumen.

*   **New Reactions - Internal Exchange (IEX):** To allow organisms to share and compete for metabolites, the `_create_inter_species_exchange()` helper function establishes a common lumen compartment `[u]`. Organism-specific metabolites are connected to this shared pool via new IEX reactions of the format `ade[u] <=> B_theta_ade[u]`. These are named e.g., `B_theta_IEX_ade[u]tr` and are fully reversible, allowing an organism to either secrete a metabolite into the common lumen or take it up.

At this stage, we now have a global community model consisting of all AGORA-native organisms found in the abundance input file.

#### # Step 1.2: Incorporating New Compartments Representing the Host-Microbiome Interface within the Global Community Model (`add_diet_fecal_compartments()`)

After all organism GEMs are assembled into a first-draft global community model, add_diet_fecal_compartments() (called by build_global_model()) processes the draft model further to add compartments representing the interface with the host.

*   **Removed Reactions:** All original exchange reactions from the base models (e.g., `EX_ade(e)`) are removed. This is a critical step to ensure all nutrient uptake and waste secretion occurs through the new lumen compartment [u].

*   **New Compartments & Reactions:** This function then introduces the diet `[d]` and fecal `[fe]` compartments. It then creates a chain of new transport reactions for every distinct metabolite in the sharied lumen, connecting the diet to the fecal output. The list of distinct metabolites in the lumen comprises all unique metabolites for which at least one species has an associated exchange reaction in individual species GEMs. For each of these lumen metabolites, the following reactions are added:
    1.  **Diet Exchange:** `EX_metabolite[d]: metabolite[d] <=> `, LB = -1000, UB = 10000
    2.  **Diet-to-Lumen Transport (DUt):** `DUt_metabolite: metabolite[d] -> metabolite[u]`, LB = 0, UB = 10000
    3.  **Lumen-to-Fecal Transport (UFEt):** `UFEt_metabolite: metabolite[u] -> metabolite[fe]`, LB = 0, UB = 10000
    4.  **Fecal Exchange:** `EX_metabolite[fe]: metabolite[fe] <=> `, LB = -1000, UB = 10000

#### # Step 1.3: Formulating Biomass Coupling Constraints (`build_global_coupling_constraints()`)

A standard constraint-based model relies on the stoichiometric matrix $S$ to enforce steady-state mass balance ($`S * v = 0`$). However, this alone is insufficient for community modeling. For example, a microbial species may be unable to grow in a community (having zero biomass flux) under certain conditions, yet can carry non-zero fluxes for exchange reactions that secrete metabolites into the lumen. This can result in unrealistic solutions particularly for inter-species metabolic interactions. To prevent this, one needs to couple the flux of each reaction within the GEM for each species to the flux of the `biomassXXX` reaction in that species as outlined in (Heinken et. al, But Micorbes 2013). This was implemented in the `build_global_coupling_constraints()` function, which introduces additional linear constraints linking the flux of every reaction within a species directly to that organism's growth rate (biomass reaction flux).

*   **Mathematical Formulation:** The formulation is based on the principles outlined in Heinken et al. Gut Microbes (2013)—with parameter $u$ set to 0 and $c$ to 400. For each reaction `i` belonging to a specific organism, its flux [Equation] is coupled to that organism' biomass reaction flux (`v_biomass`) as follows: 
    *   **For irreversible reactions:** $v_i - (c × v_{biomass}) ≤ u$
    *   **For reversible reactions:**
        * **Forward Direction:** $v_i - (c × v_{biomass}) ≤ u$
        * **Reverse Direction:** $v_i + (c × v_{biomass}) ≥ u$

*   **Implementation:** Here, $c$ is the **coupling factor** set to `400` and $u$ is a parameter that accounts for the required flux needed to maintain cellular function when the cell is not dividing if $u>0$, following Heinken et al. Gut Microbes (2013)and the `mgPipe` implementation from the MMT. Setting $u$ to 0 ensures that if `v_biomass` is 0, then `v_i` is also forced to 0, making the model more biologically realistic. These rules are compiled into a sparse **coupling matrix (C)** and associated vectors.

#### # Step 1.4: Building the Sample-Specific Models (`build_sample_model()`)

Finally, `build_sample_model` customizes the global model (integration of all organism GEMs found in abundance input file) for each individual microbiome sample based on taxonomic profiling results. This function contains the following steps:

*   **Pruning:** The `prune_zero_abundance_species()` function is used to remove all reactions and metabolites associated with species whose abundance is below a defined threshold (`1e-7`). The corresponding rows in the coupling matrix are also removed via `prune_coupling_constraints_by_species()`.
*   **Community Biomass Reaction:** The `com_biomass()` creates a biomass reaction for the community. It consumes the individual biomass metabolites of each present organism in the microbiome sample with their stoichiometric coefficients being  organism's relative abundance from taxanomic profiling results. The only metabolite on the right-hand side of this reaction is the `microbeBiomass[u]` metabolite:

    * $Σ (abundance_i * species_i\_biomass[c]) -> microbeBiomass[u]$
    * The left-hand side of the above reaction are relative abundances of organisms 1..X in the model.
    * Additionally, the `UFEt_microbeBiomass` transport reaction, $microbeBiomass[u] --> microbeBiomass[fe]$, and the `EX_microbeBiomass[fe]` fecal exchange reaction, $microbeBiomass[fe] <=> $, are added to prevent blockage of the community biomass reaction.

*   **Objective Function:** The model's objective is set to maximize the flux of the community biomass fecal exchange reaction: `model.objective = "EX_microbeBiomass[fe]"`.

The final sample-specific model, complete with all its constraints, is saved as a `.mat` file using the `make_mg_pipe_model_dict()` function.

### Stage 2: Diet-Constrained FBA Simulation and Analysis

This stage, driven by diet_adaptation.py, applies dietary constraints and simulates the metabolic potential of the microbiome using FBA simulation of the constructed microbiome GEM. The following core functions are included within this script:

#### # Step 2.1: Dietary Adjustment (`adapt_vmh_diet_to_agora`)

The `adapt_vmh_diet_to_agora()` function pre-processes an in silico diet file to make it compatible with the community models. It renames diet reaction IDs (e.g., `EX_glc(e)` becomes `Diet_EX_glc_D[d]`) and ensures feasibility by adding essential metabolites (`ESSENTIAL_METS`, `UNMAPPED_METS`) and relaxing bounds for certain `MICRONUTRIENTS`.

* ESSENTIAL_METS: This contains a predefined list of essential metabolites that are critical for the growth (biomass production) of individual species or community growth. These metabolites include water (H₂O), oxygen or alternative electron acceptors (depending on environment), Ions (e.g., Na⁺, K⁺, Cl⁻, phosphate, etc.), cofactors and trace elements, protons (H⁺), carbon sources (e.g., glucose, acetate, etc., depending on media). The full list can be found within the code. That metabolites are assumed to be available in the environment (e.g., gut lumen) And must be included in the medium, if not present in the diet to ensure feasibility and realistic behavior of flux simulations. Any missing essential metabolites is included in the diet  with an uptake limit of 0.1 (LB = -0.1).  

* UNMAPPED_METS: The Diet Designer tool takes a diet file or nutritional table (daily consumption of different food items in g or mg per person) and translates it into exchange fluxes in AGORA models and bounds on these fluxes. Each dietary compound has to be manually or programmatically mapped to one or more model exchange reactions. UNMAPPED_METS is a predefined list of exchange reactions corresponding to dietary or host-derived compounds that exist in AGORA models, but are not currently mapped in the Diet Designer (i.e., not included in the diet’s nutrient table or composition data). However, they are known to be potentially available in the gut (via food, host secretions, or microbial metabolism). These metabolties should still be allowed in the lumen, despite not being in the mapped diet, because they are known to be present in the gut as their uptake is still biologically plausible. The full list of these metabolites can be found in the codes. The metabolties are added to the diet with an uptake limit of 50 mmol/gDW/h (LB = -50).

* Choline, `EX_chol(e)`, is added to the diet as well if not already present. It's uptake limit is set to 41.251 mmol/gDW/h (LB = -41.251).

* MICRONUTRIENTS: This is a predefined list of trace compounds (e.g., vitamins, minerals, cofactors, and functional dietary fibers) that are required in very small amounts by microbes for growth or metabolic function and may be present in the diet, but are assigned extremely low uptake values (e.g., < 1e-6 mol/day/person) in the original diet definitions. If not sufficiently available, these meteabolites may artificially limit growth of individual species or communities in simulation. The full list of these metabolites is available in the code. 
    * If the uptake limit in the diet is equal to or below 0.1 mmol/gDW/h, their uptake limit is increased by a factor of 100 (LB is multiplided by 100). 
    * If the micronutrient being considered is `pnto_R` and its uptake limit is still below 0.1 mmol/gDW/h, then its uptake limit is set to 0.1, otherwise it is left alone. 
    * If the micronutrient being considered is one of ['fol', 'arab_L', 'xyl_D', 'amp', 'nh4', 'cobalt2'] and its uptake limit is less than 1 mmol/gDW/h, then its uptake limit is set to 1, otherwise it is left alone
    * These relaxations of constraints are only applied to present micronutrients, none are manually added if not present.

#### # Step 2.2a: The FBA Formulation

The `process_single_sample()` function sets up and solves a constrained Flux Balance Analysis (FBA) problem to predict the set of metabolites that a given microbiome sample can produce and secrete into lumen. The helper function `build_constraint_matrix()` illustrates how these components are assembled from the model file.

*   **Objective Function:** $Maximize c * v$, where $c$ is a vector with `1` at the position of the objective reaction (`EX_microbeBiomass[fe]`) and `0` elsewhere.
*   **Constraints:**
    1.  **Mass Balance:** $S * v = b$, where $S$ is the stoichiometric matrix.
    2.  **Flux Bounds:** $lb <= v <= ub$, where the lower bounds of diet reactions are set by the diet file.
    3.  **Coupling Constraints:** $C * v <= d$, where $C$ is the coupling matrix built in Stage 1, and $d$ is usually 0, meaning that the constraint is just a relationship between reaction fluxes.

#### # Step 2.2b: Flux Variability Analysis (FVA) - A Deep Dive

The goal of FVA is to determine the metabolic flexibility of the community by calculating the minimum and maximum the flux of each reaction while the system remains in a near-optimal state for a given objective function. The implementation of this step is crucial and is handled by` _analyze_metabolite_fluxes()`, which is called within `process_single_sample()`. This pipeline contains two parallel implementations for this task, revealing both the underlying mechanics and the practical, high-performance approach.

##### # Part A: The Manual, "From-Scratch" FVA Implementation

This approach, preserved in commented-out code, demonstrates how the FVA problem is constructed from its fundamental components.

1.  **Assembling the Full Problem (`build_constraint_matrix()`):** This function loads all mathematical components from the `.mat` file: the stoichiometric matrix (`S`), the coupling matrix (`C`), their respective right-hand-side vectors (`b` and `d`), constraint sense vectors (`csense`, `dsense`) that tell you if `constraint_i` is ≤, ≥, or = to `b` and `d`, and reaction bounds (`lb`, `ub`). It combines them into a single, large constraint matrix.
2.  **Creating a Solver-Ready Model (`build_optlang_model()`):** This function translates the raw matrices into a structured optimization model that a solver like CPLEX can understand. It explicitly defines each flux `v_i` as a variable and builds each constraint equation from the rows of the `A` matrix.
3.  **Executing Sequential FVA (`run_sequential_fva()`):** This function performs the core FVA logic:
    *   First, it solves the primary FBA problem to maximize community biomass, `Z_opt`. This is a check to see if the model is feasible
    *   It then adds a new constraint to the model: $$v_biomass >= 0.9999 * Z_opt$. This forces the model to only consider solutions that achieve near-optimal growth.
    *   Finally, it iterates through every exchange reaction of interest (EX_met[fe] reactions). For each one, it temporarily sets it as the model's objective and solves twice: once to minimize its flux and once to maximize it. This process reveals the full range of metabolic activity possible *while the community is thriving*.

##### # Part B: The Abstracted, High-Performance FVA Implementation

This is the active method used by the workflow for efficiency. It leverages the highly optimized `cobrapy` library to achieve the same result.

1.  **Applying Coupling Constraints (`couple_constraints()`):** This function is the critical bridge. It reads the coupling constraints (`C`, `d`, `dsense`) from the `.mat` file and, using the `optlang` API (a sympy-based optimization modeling language), adds them directly to the `cobra.Model` object's internal solver interface. Once this is done, the COBRA model itself is inherently aware of the coupling rules for all subsequent calculations.
2.  **Using `cobra.flux_variability_analysis()`:** With the model now fully constrained, the pipeline can call COBRApy's built-in `flux_variability_analysis()` function. This function executes the same logic as the manual method (constraining the objective and iterating) but does so far more efficiently, using optimized code and parallel processing capabilities. Essentially, this method already accounts for the model's stochiometric matrix, and we only have to construct the additional coupling constraints, rather than reconstructing all constraints as demoed above.

Both methods solve the same biological problem, but the second approach is used in practice for its performance advantage.

#### # Step 2.4: Reporting Final Fluxes

The FVA results are used by `_analyze_metabolite_fluxes()` to report the following results for each metabolite, consistent with those in the Microbiome Modeling Toolbox.

*   **Net Production:** `NetProd[met] = abs(min_flux_diet[Diet_EX_met[d]] + max_flux_fecal[EX_met[fe]])`
    * Net Production is defined as the absolute value of the difference between the maximum fecal secretion flux and the maximum dietary uptake flux for metabolite `met`. Keep in mind, When running FVA on `Diet_EX_met[d]` reactions, both min_flux_diet and max_flux_diet only contain values <= 0. Since these are negative values, min_flux_diet is referred to as maximum dietary uptake and max_flux_diet is referred to as minimum dietary uptake.
    * Note the direcetions of the exchange reactions: $EX_met[d]: met[d] <=> $
*   **Net Uptake:** `NetUptake[met] = abs(max_flux_diet[Diet_EX_met[d]] + min_flux_fecal[EX_met[fe]])`
    * Net uptake is defined as the absolute value of the difference between the minimum fecal secretion flux and the minimum dietary uptake flux for the metabolite met.
    * Again, note that the minimum dietary uptake is determined by maximizing the flux dietary exchange reaction given its reaction direction. 

Additionally, we have decided to report two new outputs:

*   **Min Net Fecal Excretion:** `MinNetFeEx[met] = min_flux_fecal[EX_met[fe]] + min_flux_diet[Diet_EX_met[d]]`
    * Represents the minimum net fecal secretion potential of the community (worst-case seceraio for secretion) and is defined min fecal secretion – max diet uptake for the metabolite `met`.

*   **Raw FVA Results:** This simply holds the raw flux values computed from running FVA on the fecal and dietary exchange reactions for each metbolite. This will report the `min_flux_diet`, `max_flux_diet`, `min_flux_fecal`, `max_flux_fecal` as described above.

The Net Production value represents the community's de novo synthesis capability for that metabolite. The final values are collected for all lumen metabolites across all samples and written to `inputDiet_net_secretion_fluxes.csv` and `inputDiet_net_uptake_fluxes.csv`. Notably, we believe that “Net Uptake” is a rather confusing term and the way it is defined, does not provide any useful biological insights, but we keep it as is for now, so our results can be compared directly with those of the Microbiome Modeling Toolbox. 

---

## 4. Downstream Analysis: Strain-level Contributions to Metabolites of Interest (`predict_microbe_contributions()`)

The `predict_microbe_contribution.py` script, through its core `predict_microbe_contributions()` function, performs targeted analysis to determine the contribution of individual microbial strains to the overall community's metabolic exchange. This is particularly useful for metabolites where the community-wide secretion potential shows significant differences between different conditions (e.g., disease cases and controls). Even if the community-wide potential does not differ, individual strain contributions might still vary significantly, providing crucial insights.

The function operates on the diet-constrained models that are generated and saved in the "Diet" subfolder during Stage 2.

*   **Workflow for Contribution Analysis:**
    *   **Input Models:** The function takes as input the personalized community models that have already been constrained with the simulated dietary regimes, found in a subfolder called "Diet" within the results folder.
    *   **Target Metabolites:** It allows defining a specific list of metabolites to analyze (e.g., acetate and formate). By default, if no list is provided, it analyzes all exchanged metabolites (specifically, Internal Exchange (IEX) reactions) within the models.
    *   **Flux Variability Analysis (FVA) on IEX Reactions:** For each diet-constrained model, FVA is performed on the internal exchange (IEX) reactions. These are reactions that connect organism-specific metabolites to the shared lumen (`species_name_IEX_met[u]tr`). The FVA determines the minimum and maximum possible flux for each of these IEX reactions while the overall community biomass production (`EX_microbeBiomass[fe]`) remains at a near-optimal level (typically 99.99% of the optimal objective) [4].
        *   Crucially, before running FVA, the `couple_constraints()` function (from `diet_adaptation.py`) is called to add the biomass coupling constraints to the model, ensuring that the FVA results reflect biologically realistic flux distributions within each strain [1][3].
    *   **Output Metrics:** The function produces three main outputs:
        *   `minFluxes`: Shows the minimum fluxes in the reverse direction through all analyzed internal exchange reactions. This often corresponds to the potential for secretion by the microbe into the lumen [4].
        *   `maxFluxes`: Shows the maximum fluxes in the forward direction through all analyzed internal exchange reactions. This often corresponds to the potential for uptake by the microbe from the lumen [4].
        *   `fluxSpans`: Represents the distance between the minimal and maximal fluxes for each internal exchange reaction, indicating the metabolic flexibility of that particular exchange [4].

By performing this analysis, users can gain detailed insights into the specific roles individual strains play in the overall metabolic landscape of the microbiome, contributing to the uptake or secretion of key metabolites under various dietary conditions.


## References

1.  Heinken, A., Baldini, F., Heirendt, L., Magnusdottir, S., Fleming, R. M. T., & Thiele, I. (2019). The Microbiome Modeling Toolbox: from microbial interactions to personalized microbial communities. *Bioinformatics*, 35(13), 2332-2334.
2.  Heinken, A., Sahoo, S., Fleming, R. M. T., & Thiele, I. (2013). Systems-level characterization of a host-microbe metabolic symbiosis in the mammalian gut. *Gut Microbes*, 4(1), 28-40.