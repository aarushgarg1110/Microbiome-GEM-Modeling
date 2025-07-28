"""
Diet Adapter Module for MiGEMox Pipeline

This module focuses on preprocessing and adapting dietary constraint files,
such as those from VMH, to be compatible with the AGORA-based microbiome models.
It handles normalization, addition of essential metabolites, and adjustment
of micronutrient bounds.
"""

import pandas as pd
import os

# adapt_vmh_diet_to_agora: Essential AGORA metabolites
ESSENTIAL_METS = [
        'EX_12dgr180(e)', 'EX_26dap_M(e)', 'EX_2dmmq8(e)', 'EX_2obut(e)', 'EX_3mop(e)', 'EX_4abz(e)',
        'EX_4hbz(e)', 'EX_ac(e)', 'EX_acgam(e)', 'EX_acmana(e)', 'EX_acnam(e)', 'EX_ade(e)', 'EX_adn(e)',
        'EX_adocbl(e)', 'EX_ala_D(e)', 'EX_ala_L(e)', 'EX_amet(e)', 'EX_amp(e)', 'EX_arab_D(e)', 'EX_arab_L(e)',
        'EX_arg_L(e)', 'EX_asn_L(e)', 'EX_btn(e)', 'EX_ca2(e)', 'EX_cbl1(e)', 'EX_cgly(e)', 'EX_chor(e)',
        'EX_chsterol(e)', 'EX_cit(e)', 'EX_cl(e)', 'EX_cobalt2(e)', 'EX_csn(e)', 'EX_cu2(e)', 'EX_cys_L(e)',
        'EX_cytd(e)', 'EX_dad_2(e)', 'EX_dcyt(e)', 'EX_ddca(e)', 'EX_dgsn(e)', 'EX_fald(e)', 'EX_fe2(e)',
        'EX_fe3(e)', 'EX_fol(e)', 'EX_for(e)', 'EX_gal(e)', 'EX_glc_D(e)', 'EX_gln_L(e)', 'EX_glu_L(e)',
        'EX_gly(e)', 'EX_glyc(e)', 'EX_glyc3p(e)', 'EX_gsn(e)', 'EX_gthox(e)', 'EX_gthrd(e)', 'EX_gua(e)',
        'EX_h(e)', 'EX_h2o(e)', 'EX_h2s(e)', 'EX_his_L(e)', 'EX_hxan(e)', 'EX_ile_L(e)', 'EX_k(e)', 'EX_lanost(e)',
        'EX_leu_L(e)', 'EX_lys_L(e)', 'EX_malt(e)', 'EX_met_L(e)', 'EX_mg2(e)', 'EX_mn2(e)', 'EX_mqn7(e)',
        'EX_mqn8(e)', 'EX_nac(e)', 'EX_ncam(e)', 'EX_nmn(e)', 'EX_no2(e)', 'EX_ocdca(e)', 'EX_ocdcea(e)',
        'EX_orn(e)', 'EX_phe_L(e)', 'EX_pheme(e)', 'EX_pi(e)', 'EX_pnto_R(e)', 'EX_pro_L(e)', 'EX_ptrc(e)',
        'EX_pydx(e)', 'EX_pydxn(e)', 'EX_q8(e)', 'EX_rib_D(e)', 'EX_ribflv(e)', 'EX_ser_L(e)', 'EX_sheme(e)',
        'EX_so4(e)', 'EX_spmd(e)', 'EX_thm(e)', 'EX_thr_L(e)', 'EX_thymd(e)', 'EX_trp_L(e)', 'EX_ttdca(e)',
        'EX_tyr_L(e)', 'EX_ura(e)', 'EX_val_L(e)', 'EX_xan(e)', 'EX_xyl_D(e)', 'EX_zn2(e)', 'EX_glu_D(e)',
        'EX_melib(e)', 'EX_chtbs(e)', 'EX_metsox_S_L(e)', 'EX_hdca(e)', 'EX_gam(e)', 'EX_indole(e)', 'EX_glcn(e)',
        'EX_coa(e)', 'EX_man(e)', 'EX_fum(e)', 'EX_succ(e)', 'EX_no3(e)', 'EX_ins(e)', 'EX_uri(e)', 'EX_drib(e)',
        'EX_pime(e)', 'EX_lac_L(e)', 'EX_glypro(e)', 'EX_urea(e)', 'EX_duri(e)', 'EX_h2(e)', 'EX_mal_L(e)',
        'EX_tre(e)', 'EX_orot(e)', 'EX_glymet(e)', 'EX_glyleu(e)', 'EX_pydx5p(e)', 'EX_so3(e)', 'EX_nh4(e)'
    ]

# adapt_vmh_diet_to_agora: Allow uptake of certain dietary compounds that are currently not mapped in the Diet Designer
UNMAPPED_METS = [
        'EX_asn_L(e)', 'EX_gln_L(e)', 'EX_crn(e)', 'EX_elaid(e)', 'EX_hdcea(e)', 'EX_dlnlcg(e)', 'EX_adrn(e)',
        'EX_hco3(e)', 'EX_sprm(e)', 'EX_carn(e)', 'EX_7thf(e)', 'EX_Lcystin(e)', 'EX_hista(e)', 'EX_orn(e)',
        'EX_ptrc(e)', 'EX_creat(e)', 'EX_cytd(e)', 'EX_so4(e)'
    ]

# adapt_vmh_diet_to_agora: Increase the uptake rate of micronutrients with too low defined uptake rates to sustain microbiota model growth (below 1e-6 mol/day/person).
# Lower bounds will be relaxed by factor 100 if allowed uptake is below 0.1 mmol*gDW-1*hr-1.
MICRONUTRIENTS = [
        'EX_adocbl(e)', 'EX_vitd2(e)', 'EX_vitd3(e)', 'EX_psyl(e)', 'EX_gum(e)', 'EX_bglc(e)', 'EX_phyQ(e)',
        'EX_fol(e)', 'EX_5mthf(e)', 'EX_q10(e)', 'EX_retinol_9_cis(e)', 'EX_pydxn(e)', 'EX_pydam(e)', 'EX_pydx(e)',
        'EX_pheme(e)', 'EX_ribflv(e)', 'EX_thm(e)', 'EX_avite1(e)', 'EX_pnto_R(e)', 'EX_na1(e)', 'EX_cl(e)',
        'EX_k(e)', 'EX_pi(e)', 'EX_zn2(e)', 'EX_cu2(e)', 'EX_btn(e)'
    ]

def adapt_vmh_diet_to_agora(vmh_diet_file: str, setup_type: str = 'Microbiota') -> pd.DataFrame:
    """
    Adapt VMH diet file to AGORA microbiota model constraints.
    
    Biological Process:
    1. Load and normalize VMH diet exchange reaction bounds
    2. Add essential metabolites required for AGORA model feasibility  
    3. Adjust micronutrient constraints for gut environment
    4. Transform reaction IDs based on model setup type
    
    Args:
        vmh_diet_file: Path to VMH diet file (.txt or .csv)
        setup_type: Model setup ('AGORA', 'Pairwise', or 'Microbiota')
        
    Returns:
        DataFrame with adapted diet constraints (rxn_id, lower_bound, upper_bound)
    """
    # Step 1: Load and preprocess diet data
    diet_df = _load_diet_data(vmh_diet_file)

    # Step 2: Add essential AGORA metabolites
    present_rxns = set(diet_df['rxn_id'])
    missing_mets = sorted(set(ESSENTIAL_METS) - present_rxns)
    additional_rows = [{'rxn_id': rxn, 'lower_bound': -0.1, 'upper_bound': 0} for rxn in missing_mets]
    diet_df = pd.concat([diet_df, pd.DataFrame(additional_rows)], ignore_index=True)

    # Step 3: Add unmapped compounds
    missing_unmapped = sorted(set(UNMAPPED_METS) - present_rxns)
    unmapped_rows = [{'rxn_id': rxn, 'lower_bound': -50.0, 'upper_bound': 0} for rxn in missing_unmapped]
    diet_df = pd.concat([diet_df, pd.DataFrame(unmapped_rows)], ignore_index=True)

    # Add choline if missing
    if 'EX_chol(e)' not in present_rxns:
        diet_df = pd.concat([diet_df, pd.DataFrame([{
            'rxn_id': 'EX_chol(e)', 'lower_bound': -41.251, 'upper_bound': 0
        }])], ignore_index=True)

    # Step 4: Relax micronutrient constraints
    diet_df = _relax_micronutrient_constraints(diet_df)

    # Step 5: Apply setup-specific transformations
    diet_df = _apply_setup_transformations(diet_df, setup_type)

    return diet_df[['rxn_id', 'lower_bound', 'upper_bound']]

def _load_diet_data(vmh_diet_file: str) -> pd.DataFrame:
    """Load VMH diet file and perform initial preprocessing."""
    diet_df = pd.read_csv(vmh_diet_file, header=0, sep='\t', names=['rxn_id', 'lower_bound'])
    
    # Store original bounds and flip sign for uptake (VMH convention)
    diet_df['original_lb'] = diet_df['lower_bound'].astype(float)
    diet_df['lower_bound'] = -diet_df['lower_bound'].astype(float)
    
    # Normalize exchange reaction IDs for AGORA compatibility
    diet_df['rxn_id'] = diet_df['rxn_id'].str.replace('[e]', '(e)', regex=False)
    
    # Apply known metabolite ID mappings
    id_mappings = {
        'EX_adpcbl(e)': 'EX_adocbl(e)',  # Adenosylcobalamin
        'EX_glc(e)': 'EX_glc_D(e)',      # D-Glucose
        'EX_sbt-d(e)': 'EX_sbt_D(e)'     # D-Sorbitol
    }
    diet_df['rxn_id'] = diet_df['rxn_id'].replace(id_mappings)
    
    return diet_df

def _relax_micronutrient_constraints(diet_df: pd.DataFrame) -> pd.DataFrame:
    """Relax constraints on micronutrients to prevent infeasibility in gut environment."""
    def relax_limit(row):
        rxn = row['rxn_id']
        lb = row['lower_bound']
        if rxn in MICRONUTRIENTS and abs(lb) <= 0.1:
            lb *= 100
        if rxn == 'EX_pnto_R(e)' and abs(lb) < 0.1:
            return -0.1
        if rxn in ['EX_fol(e)', 'EX_arab_L(e)', 'EX_xyl_D(e)', 'EX_amp(e)', 'EX_nh4(e)', 'EX_cobalt2(e)'] and abs(lb) < 1:
            return -1
        return lb

    diet_df['lower_bound'] = diet_df.apply(relax_limit, axis=1)
    return diet_df

def _apply_setup_transformations(diet_df: pd.DataFrame, setup_type: str) -> pd.DataFrame:
    """Apply model setup-specific ID transformations."""
    if setup_type == 'Pairwise':
        diet_df['rxn_id'] = diet_df['rxn_id'].str.replace('(e)', '[u]', regex=False)
    elif setup_type == 'Microbiota':
        # First, add an upper_bound column with default of 0: no secretion
        diet_df['upper_bound'] = 0.0
        
        # For each row that matches an original diet constraint, 
        # set upper bound to -0.8 * original_lower_bound
        # This allows some secretion but limits it to 80% of the uptake rate
        for index, row in diet_df.iterrows():
            # Only negative lower bounds (uptake reactions) get this treatment
            if row['lower_bound'] < 0 and pd.notnull(row['original_lb']):
                diet_df.at[index, 'upper_bound'] = -0.8 * row['original_lb']

        # Now rename all reactions to use Diet_EX_ prefix and [d] compartment
        diet_df['rxn_id'] = diet_df['rxn_id'].str.replace('EX_', 'Diet_EX_')
        diet_df['rxn_id'] = diet_df['rxn_id'].str.replace('(e)', '[d]')
    
    return diet_df[['rxn_id', 'lower_bound', 'upper_bound']]