"""
Analysis Module for MiGEMox Pipeline

This module provides functions for collecting and organizing simulation
results, such as net production and net uptake fluxes, into structured
data formats (e.g., Pandas DataFrames) for easier interpretation and
further downstream analysis.
"""

import pandas as pd
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def collect_flux_profiles(samp_names: List[str], exchanges: List[str],
                          net_production: Dict[str, Dict[str, float]],
                          net_uptake: Dict[str, Dict[str, float]],
                          res_path: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Converts net production and net uptake dictionaries, which contain results
    for each sample, into structured Pandas DataFrames.

    Args:
        res_path: path to results directory where flux.csv files should be stored
        samp_names: List of sample identifiers.
        exchanges: List of exchange reaction IDs (metabolites) analyzed.
        net_production: Dictionary where keys are sample names and values are
                        dictionaries of metabolite-to-flux for net production.
        net_uptake: Dictionary where keys are sample names and values are
                    dictionaries of metabolite-to-flux for net uptake.

    Returns:
        Tuple of:
            - net_secretion_df (pd.DataFrame): DataFrame of net production (secretion)
                                               fluxes with metabolites as index and samples as columns.
            - net_uptake_df (pd.DataFrame): DataFrame of net uptake fluxes with
                                            metabolites as index and samples as columns.
    """
    # Ensure exchanges are sorted for consistent DataFrame indexing
    exchanges_sorted = sorted(exchanges)

    res_path = Path.cwd() / 'Results' if not res_path else Path(res_path)
    logger.info(f"Results will be saved to: {res_path}")

    # Prepare data for net production DataFrame
    prod_data, uptk_data = {}, {}
    for samp in samp_names:
        prod_data[samp] = [net_production.get(samp, {}).get(r, 0.0) for r in exchanges_sorted]
        uptk_data[samp] = [net_uptake.get(samp, {}).get(r, 0.0) for r in exchanges_sorted]
    
    net_secretion_df = pd.DataFrame(prod_data, index=exchanges_sorted)
    net_uptake_df = pd.DataFrame(uptk_data, index=exchanges_sorted)

    net_secretion_df.to_csv(res_path / 'inputDiet_net_secretion_fluxes.csv', index=True, index_label='Net secretion')
    net_uptake_df.to_csv(res_path / 'inputDiet_net_uptake_fluxes.csv', index=True, index_label='Net uptake')


    return net_secretion_df, net_uptake_df