"""
MiGEMox Pipeline Main Entry Point

This script orchestrates the entire MiGEMox workflow, from model construction
and diet-constrained simulation to downstream analysis of strain contributions.
It serves as the high-level control script, calling functions from
various specialized modules.
"""

import pandas as pd
import os
import argparse
from typing import Optional, List
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

import docplex
import cplex

# Import functions from our new modules
from src.pipeline.io_utils import get_individual_size_name
from src.pipeline.model_builder import build_mgpipe_models
from src.pipeline.diet_adapter import adapt_vmh_diet_to_agora
from src.pipeline.simulation import simulate_microbiota_models
from src.pipeline.analysis import collect_flux_profiles
from src.downstream_analysis.predict_microbe_contribution import predict_microbe_contributions

def run_migemox_pipeline(abun_filepath: str, mod_filepath: str, res_filepath: str,
                         diet_filepath: str, workers: int = 1, solver: str = 'cplex',
                         biomass_bounds: tuple = (0.4, 1.0), analyze_contributions: bool = False):
    """
    Main function to run the MiGEMox pipeline.

    Args:
        abun_filepath: Path to the abundance CSV file.
        mod_filepath: Path to the directory containing organism model files (.mat).
        res_filepath: Base directory for saving output models and results.
        diet_filepath: Path to the VMH diet file.
        workers: Number of parallel workers to use for sample processing.
        solver: Optimization solver to use (e.g., 'cplex', 'gurobi').
        biomass_bounds: Tuple (lower, upper) for community biomass reaction.
        analyze_contributions: Boolean, whether to run strain contribution analysis.
    """
    print("--- MiGEMox Pipeline Started ---")
    
    build_mgpipe_models(
        abun_filepath=abun_filepath,
        mod_filepath=mod_filepath,
        out_filepath=f'{res_filepath}/Personalized_Models',
        workers=workers
    )

    clean_samp_names, organisms, ex_mets = get_individual_size_name(
        abun_file_path=abun_filepath,
        mod_path=mod_filepath
    )
    
    # 2. Adapt Diet
    print("\n--- Stage 2: Adapting Diet and Running Simulations ---")

    # 3. Simulate Microbiota Models
    exchanges, net_production, net_uptake = simulate_microbiota_models(
        sample_names=clean_samp_names,
        ex_mets=ex_mets,
        model_dir=f'{res_filepath}/Personalized_Models',
        diet_file=diet_filepath,
        res_path=res_filepath,
        biomass_bounds=biomass_bounds,
        solver=solver,
        workers=2
    )

    # 4. Collect Flux Profiles and Save
    print("\n--- Collecting and Saving Simulation Results ---")
    net_secretion_df, net_uptake_df = collect_flux_profiles(
        samp_names=clean_samp_names,
        exchanges=exchanges,
        net_production=net_production,
        net_uptake=net_uptake,
        res_path=res_filepath
    )

    print(f"Net secretion and uptake results saved to {res_filepath}.")

    # 5. Run Strain Contribution Analysis (Optional)
    if analyze_contributions:
        print("\n--- Downstream Analysis: Predicting Strain Contributions ---")

        # VMH_ID from met names like EX_ac[fe] -> ac
        mets = [x.split('[')[0] for x in ex_mets]
        diet_mod_dir_for_contributions = os.path.join(res_filepath, 'Diet')
        
        min_fluxes_df, max_fluxes_df, flux_spans_df = predict_microbe_contributions(
            diet_mod_dir=diet_mod_dir_for_contributions,
            mets_list=mets,
            solver=solver,
            workers=workers
        )
        print("Strain contribution analysis completed and results saved.")

    print("\n--- MiGEMox Pipeline Finished ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the MiGEMox (Microbiome Genome-Scale Modeling) pipeline.")
    parser.add_argument("--abun_filepath", type=str, required=True,
                        help="Path to the abundance CSV file (e.g., test_data_input/normCoverageReduced.csv)")
    parser.add_argument("--mod_filepath", type=str, required=True,
                        help="Path to the directory containing organism model files (.mat) (e.g., test_data_input/AGORA103)")
    parser.add_argument("--res_filepath", type=str, default="Results",
                        help="Base directory for saving all output models and results (e.g., Python_Models)")
    parser.add_argument("--diet_filepath", type=str, required=True,
                        help="Path to the VMH diet file (e.g., test_data_input/AverageEU_diet_fluxes.txt)")
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of parallel workers to use for sample processing (default: 1)")
    parser.add_argument("--solver", type=str, default="cplex",
                        help="Optimization solver to use (e.g., 'cplex', 'gurobi', 'glpk') (default: cplex)")
    parser.add_argument("--biomass_lb", type=float, default=0.4,
                        help="Lower bound for community biomass growth (default: 0.4)")
    parser.add_argument("--biomass_ub", type=float, default=1.0,
                        help="Upper bound for community biomass growth (default: 1.0)")
    parser.add_argument("--analyze_contributions", action="store_true",
                        help="Set this flag to run strain contribution analysis.")

    args = parser.parse_args()

    # Combine biomass bounds
    biomass_bounds_tuple = (args.biomass_lb, args.biomass_ub)

    run_migemox_pipeline(
        abun_filepath=args.abun_filepath,
        mod_filepath=args.mod_filepath,
        res_filepath=args.res_filepath,
        diet_filepath=args.diet_filepath,
        workers=args.workers,
        solver=args.solver,
        biomass_bounds=biomass_bounds_tuple,
        analyze_contributions=args.analyze_contributions
    )

# Example:
# run_migemox_pipeline(abun_filepath="test_data_input/normCoverageReduced.csv",
#                      mod_filepath='test_data_input/AGORA103',
#                      res_filepath='Results',
#                      diet_filepath='test_data_input/AverageEU_diet_fluxes.txt',
#                      workers=2,
#                      solver='cplex',
#                      biomass_bounds=(0.4, 1.0),
#                      analyze_contributions=True)