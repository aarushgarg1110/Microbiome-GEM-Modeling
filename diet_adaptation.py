import cobra
import multiprocessing
from cobra.io import save_json_model, load_matlab_model, load_json_model
from multiprocessing import Pool
import os
import pandas as pd
import scipy
import docplex
import cplex
from tqdm import tqdm
import sys
import sklearn 
import matplotlib 
import joblib
from concurrent.futures import ProcessPoolExecutor, as_completed
from compy import compy

def get_individual_size_name(abun_file_path, mod_path):
    # Step 1: Read abundance CSV
    try:
        df = pd.read_csv(abun_file_path, index_col=0)
    except Exception as e:
        raise ValueError(f"Error reading file: {e}")

    samp_names = list(df.columns)
    organisms = list(df.index)

    # Step 2: Clean sample names to be valid Python identifiers
    clean_samp_names = []
    for name in samp_names:
        valid_name = name
        if not valid_name.isidentifier():
            import re
            valid_name = re.sub('[^0-9a-zA-Z_]', '_', valid_name)
            if not valid_name[0].isalpha():
                valid_name = 'sample_' + valid_name
        clean_samp_names.append(valid_name)

    # Step 3: Load models and extract [e] metabolites
    models = []
    ex_mets = set()

    for organism in organisms:
        model_file = os.path.join(mod_path, organism + '.mat')
        if not os.path.exists(model_file):
            raise FileNotFoundError(f"Model file not found: {model_file}")
        
        try:
            model = load_matlab_model(model_file)
        except Exception as e:
            raise ValueError(f"Error loading model {organism}: {e}")

        models.append(model)

        # Extract extracellular metabolites (assuming suffix '[e]')
        mets = [met.id for met in model.metabolites if met.id.endswith('[e]')]
        ex_mets.update(mets)

    return clean_samp_names, organisms, list(sorted(ex_mets))

def adapt_vmh_diet_to_agora(vmh_diet_file, setup_used='Microbiota'):
    """
    Adapt a VMH diet file to make it compatible with AGORA microbiota models.

    Args:
        vmh_diet_file (str): Path to the VMH diet file (.txt or .csv).
        setup_used (str): One of ['AGORA', 'Pairwise', 'Microbiota'].

    Returns:
        pd.DataFrame: Adapted diet with columns ['rxn_id', 'lower_bound', 'upper_bound'].
    """
    # Load diet
    diet_df = pd.read_csv(vmh_diet_file, header=0, sep='\t', names=['rxn_id', 'lower_bound'])
    diet_df['original_lb'] = diet_df['lower_bound'].astype(float)
    diet_df['lower_bound'] = -diet_df['lower_bound'].astype(float)  # Flip sign for uptake

    # Normalize exchange IDs
    diet_df['rxn_id'] = diet_df['rxn_id'].str.replace('[e]', '(e)', regex=False)
    diet_df['rxn_id'] = diet_df['rxn_id'].replace({
        'EX_adpcbl(e)': 'EX_adocbl(e)',
        'EX_glc(e)': 'EX_glc_D(e)',
        'EX_sbt-d(e)': 'EX_sbt_D(e)'
    })

    # Essential AGORA metabolites
    essential_mets = [
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

    # Add missing essential metabolites
    present_rxns = set(diet_df['rxn_id'])
    missing_mets = sorted(set(essential_mets) - present_rxns)
    additional_rows = [{'rxn_id': rxn, 'lower_bound': -0.1, 'upper_bound': 0} for rxn in missing_mets]
    diet_df = pd.concat([diet_df, pd.DataFrame(additional_rows)], ignore_index=True)

    # Add unmapped compounds
    unmapped = [
        'EX_asn_L(e)', 'EX_gln_L(e)', 'EX_crn(e)', 'EX_elaid(e)', 'EX_hdcea(e)', 'EX_dlnlcg(e)', 'EX_adrn(e)',
        'EX_hco3(e)', 'EX_sprm(e)', 'EX_carn(e)', 'EX_7thf(e)', 'EX_Lcystin(e)', 'EX_hista(e)', 'EX_orn(e)',
        'EX_ptrc(e)', 'EX_creat(e)', 'EX_cytd(e)', 'EX_so4(e)'
    ]
    missing_unmapped = sorted(set(unmapped) - present_rxns)
    unmapped_rows = [{'rxn_id': rxn, 'lower_bound': -50.0, 'upper_bound': 0} for rxn in missing_unmapped]
    diet_df = pd.concat([diet_df, pd.DataFrame(unmapped_rows)], ignore_index=True)

    # Add choline if missing
    if 'EX_chol(e)' not in present_rxns:
        diet_df = pd.concat([diet_df, pd.DataFrame([{
            'rxn_id': 'EX_chol(e)', 'lower_bound': -41.251, 'upper_bound': 0
        }])], ignore_index=True)

    # Adjust micronutrient uptake rates
    micronutrients = [
        'EX_adocbl(e)', 'EX_vitd2(e)', 'EX_vitd3(e)', 'EX_psyl(e)', 'EX_gum(e)', 'EX_bglc(e)', 'EX_phyQ(e)',
        'EX_fol(e)', 'EX_5mthf(e)', 'EX_q10(e)', 'EX_retinol_9_cis(e)', 'EX_pydxn(e)', 'EX_pydam(e)', 'EX_pydx(e)',
        'EX_pheme(e)', 'EX_ribflv(e)', 'EX_thm(e)', 'EX_avite1(e)', 'EX_pnto_R(e)', 'EX_na1(e)', 'EX_cl(e)',
        'EX_k(e)', 'EX_pi(e)', 'EX_zn2(e)', 'EX_cu2(e)', 'EX_btn(e)'
    ]

    def relax_limit(row):
        rxn = row['rxn_id']
        lb = row['lower_bound']
        if rxn in micronutrients and abs(lb) <= 0.1:
            lb *= 100
        if rxn == 'EX_pnto_R(e)' and abs(lb) < 0.1:
            return -0.1
        if rxn in ['EX_fol(e)', 'EX_arab_L(e)', 'EX_xyl_D(e)', 'EX_amp(e)', 'EX_nh4(e)', 'EX_cobalt2(e)'] and abs(lb) < 1:
            return -1
        return lb

    diet_df['lower_bound'] = diet_df.apply(relax_limit, axis=1)

    # Add setup-specific transformations
    if setup_used == 'Pairwise':
        diet_df['rxn_id'] = diet_df['rxn_id'].str.replace('(e)', '[u]', regex=False)
    elif setup_used == 'Microbiota':
        # First, add an upper_bound column with default of 0
        diet_df['upper_bound'] = 0
        
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

def run_fva(model, reaction_ids, fraction=0.9999):
    """
    Runs Flux Variability Analysis (FVA) on selected reactions.

    Args:
        model (cobra.Model): COBRA model object
        reaction_ids (list): List of reaction IDs to analyze
        fraction (float): Optimality threshold for FVA

    Returns:
        min_fluxes (dict), max_fluxes (dict)
    """
    result = cobra.flux_analysis.flux_variability_analysis(model, reaction_list=reaction_ids, fraction_of_optimum=fraction)
    return result['minimum'].to_dict(), result['maximum'].to_dict()

def create_rxn(rxn_identifier, name, subsystem, lb, ub) -> cobra.Reaction:
    rxn = cobra.Reaction(rxn_identifier)
    rxn.name = name
    rxn.subsystem = subsystem
    rxn.lower_bound = lb
    rxn.upper_bound = ub

    return rxn

# def process_single_sample(samp, model_dir, diet_constraints, exchanges, res_path, 
#                          lower_bm, upper_bm, solver, humanMets):
#     """
#     Process a single sample: load model, apply diet, run FVA, calculate net fluxes.
    
#     Args:
#         samp (str): Sample name
#         model_dir (str): Directory containing models
#         diet_constraints (pd.DataFrame): Diet constraints
#         exchanges (list): List of exchange reactions
#         res_path (str): Results path
#         lower_bm, upper_bm (float): Community biomass bounds
#         solver (str): Solver name
#         humanMets (dict): Human metabolites dict
        
#     Returns:
#         tuple: (samp, net_production_dict, net_uptake_dict)
#     """
#     try:
#         model_path = os.path.join(model_dir, f"{samp}_communitymodel_final.json")
#         model = load_json_model(model_path)
#         model.solver = solver

#         print(f"Processing {samp}: got model")

#         # Before applying diet constraints, update reaction IDs in the model
#         diet_rxns = [r.id for r in model.reactions if '[d]' in r.id and r.id.startswith('EX_')]
#         for rxn_id in diet_rxns:
#             new_id = rxn_id.replace('EX_', 'Diet_EX_')
#             if new_id not in model.reactions:
#                 model.reactions.get_by_id(rxn_id).id = new_id

#         # First: Set ALL Diet_EX_ reactions to lower bound 0 (like useDiet does)
#         for rxn in model.reactions:
#             if rxn.id.startswith('Diet_EX_'):
#                 rxn.lower_bound = 0

#         # Apply diet
#         for _, row in diet_constraints.iterrows():
#             rxn = row['rxn_id']
#             if rxn in model.reactions:
#                 model.reactions.get_by_id(rxn).lower_bound = float(row['lower_bound'])
#                 if pd.notnull(row['upper_bound']):
#                     model.reactions.get_by_id(rxn).upper_bound = float(row['upper_bound'])

#         print(f"Processing {samp}: diet applied")
        
#         # Constrain community biomass
#         if 'communityBiomass' in model.reactions:
#             model.reactions.communityBiomass.lower_bound = lower_bm
#             model.reactions.communityBiomass.upper_bound = upper_bm

#         for rxn in model.reactions:
#             if rxn.id.startswith('UFEt_') or rxn.id.startswith('DUt_') or rxn.id.startswith('EX_'):
#                 rxn.upper_bound = 1e6

#         # Include the human mets if not already included in the diet
#         for met_id, bound in humanMets.items():
#             rxn_id = f'Diet_EX_{met_id}[d]'
#             if rxn_id not in diet_constraints['rxn_id'].values:  
#                 rxn = create_rxn(rxn_id, 'Bacteroides_sp_2_1_33B_biomass525', '', bound, 10000.)
#                 model.add_reactions([rxn])
#             model.reactions.get_by_id(rxn.id).bounds = bound, 10000.

#         # close demand and limit sink reactions
#         for rxn in model.reactions:
#             if '_DM_' in rxn.id:
#                 rxn.lower_bound = 0
#             elif '_sink_' in rxn.id:
#                 rxn.lower_bound = -1 

#         # Objective: EX_microbeBiomass[fe]
#         model.objective = 'EX_microbeBiomass[fe]'
#         model.optimize()
#         print(f"Processing {samp}: model optimized")
        
#         # Save the diet-adapted model
#         save_dir = os.path.join(res_path, 'Diet')
#         os.makedirs(save_dir, exist_ok=True)
#         diet_model_path = os.path.join(save_dir, f"pymicrobiota_model_diet_{samp}.json")
#         cobra.io.save_json_model(model, diet_model_path)

#         print(f"Processing {samp}: starting fva")

#         # Run FVA
#         fecal_rxns = [r.id for r in model.reactions if '[fe]' in r.id]
#         diet_rxns = [r.id for r in model.reactions if '[d]' in r.id]

#         min_flux_fecal, max_flux_fecal = run_fva(model, fecal_rxns)
#         print(f"Processing {samp}: starting diet fva")
#         min_flux_diet, max_flux_diet = run_fva(model, diet_rxns)

#         net_production_samp = {}
#         net_uptake_samp = {}

#         # model.exchanges contains only EX_rxns in fecal compartment

#         # exchanges derived from exMets (all exchanged metabolites across all individual models) -> intersect it with rxns in this particular model
#         exchanges = set(fecal_rxns).intersection(set(exchanges))

#         for rxn in exchanges:
#             fecal = rxn
#             diet = rxn.replace('EX_', 'Diet_EX_').replace('[fe]', '[d]')
#             prod = min_flux_diet.get(diet, 0) + max_flux_fecal.get(fecal, 0)
#             uptk = max_flux_diet.get(diet, 0) + min_flux_fecal.get(fecal, 0)
#             net_production_samp[rxn] = prod
#             net_uptake_samp[rxn] = uptk

#         print(f"Processing {samp}: completed")
#         return samp, net_production_samp, net_uptake_samp
        
#     except Exception as e:
#         print(f"Error processing sample {samp}: {str(e)}")
#         raise e


# def simulate_microbiota_models(
#     samp_names, ex_mets, model_dir, diet_file, res_path,
#     lower_bm=0.4, upper_bm=1.0, solver='cplex'
# ):
#     """
#     Apply diet and simulate all microbiota models with FVA using multiprocessing.

#     Args:
#         samp_names (list): List of individual sample names
#         ex_mets (list): List of extracellular metabolites
#         model_dir (str): Folder containing sample .mat models
#         diet_file (str): Path to VMH diet file
#         res_path (str): Folder to save FVA results
#         lower_bm, upper_bm (float): Community biomass bounds

#     Returns:
#         net_production (dict), net_uptake (dict)
#     """
#     os.makedirs(res_path, exist_ok=True)
#     exchanges = [f"EX_{m.replace('[e]', '[fe]')}" for m in ex_mets if m != 'biomass[e]']

#     net_production = {}
#     net_uptake = {}

#     # define human-derived metabolites present in the gut: primary bile acids, amines, mucins, host glycans
#     humanMets = {
#         'gchola': -10, 'tdchola': -10, 'tchola': -10, 'dgchol': -10,
#         '34dhphe': -10, '5htrp': -10, 'Lkynr': -10, 'f1a': -1,
#         'gncore1': -1, 'gncore2': -1, 'dsT_antigen': -1, 'sTn_antigen': -1,
#         'core8': -1, 'core7': -1, 'core5': -1, 'core4': -1,
#         'ha': -1, 'cspg_a': -1, 'cspg_b': -1, 'cspg_c': -1,
#         'cspg_d': -1, 'cspg_e': -1, 'hspg': -1
#     }

#     # Adapt diet
#     diet_constraints = adapt_vmh_diet_to_agora(diet_file, setup_used='Microbiota')

#     print("Got constraints, starting parallel processing")

#     # Use ProcessPoolExecutor for parallel processing
#     with ProcessPoolExecutor() as executor:
#         # Submit all sample processing jobs
#         futures = [
#             executor.submit(
#                 process_single_sample, 
#                 samp, model_dir, diet_constraints, exchanges, res_path,
#                 lower_bm, upper_bm, solver, humanMets
#             ) 
#             for samp in samp_names
#         ]
        
#         # Collect results as they complete
#         for future in tqdm(as_completed(futures), total=len(futures), desc='Processing samples'):
#             samp, net_production_samp, net_uptake_samp = future.result()
#             net_production[samp] = net_production_samp
#             net_uptake[samp] = net_uptake_samp

#     print("All samples processed successfully")
#     return exchanges, net_production, net_uptake

def simulate_microbiota_models(
    samp_names, ex_mets, model_dir, diet_file, res_path,
    lower_bm=0.4, upper_bm=1.0, solver='cplex'
):
    """
    Apply diet and simulate all microbiota models with FVA.

    Args:
        samp_names (list): List of individual sample names
        ex_mets (list): List of extracellular metabolites
        model_dir (str): Folder containing sample .mat models
        diet_file (str): Path to VMH diet file
        res_path (str): Folder to save FVA results
        lower_bm, upper_bm (float): Community biomass bounds

    Returns:
        net_production (dict), net_uptake (dict)
    """
    os.makedirs(res_path, exist_ok=True)
    exchanges = [f"EX_{m.replace('[e]', '[fe]')}" for m in ex_mets if m != 'biomass[e]']

    net_production = {}
    net_uptake = {}

    # define human-derived metabolites present in the gut: primary bile acids, amines, mucins, host glycans
    humanMets = {
        'gchola': -10, 'tdchola': -10, 'tchola': -10, 'dgchol': -10,
        '34dhphe': -10, '5htrp': -10, 'Lkynr': -10, 'f1a': -1,
        'gncore1': -1, 'gncore2': -1, 'dsT_antigen': -1, 'sTn_antigen': -1,
        'core8': -1, 'core7': -1, 'core5': -1, 'core4': -1,
        'ha': -1, 'cspg_a': -1, 'cspg_b': -1, 'cspg_c': -1,
        'cspg_d': -1, 'cspg_e': -1, 'hspg': -1
    }

    # Adapt diet
    diet_constraints = adapt_vmh_diet_to_agora(diet_file, setup_used='Microbiota')

    print("got constraints")

    for samp in tqdm(samp_names):
        model_path = os.path.join(model_dir, f"{samp}_communitymodel_final.json")
        model = load_json_model(model_path)
        model.solver = solver

        print(f"Processing {samp}: got model")

        # Before applying diet constraints, update reaction IDs in the model
        diet_rxns = [r.id for r in model.reactions if '[d]' in r.id and r.id.startswith('EX_')]
        for rxn_id in diet_rxns:
            new_id = rxn_id.replace('EX_', 'Diet_EX_')
            if new_id not in model.reactions:
                model.reactions.get_by_id(rxn_id).id = new_id

        # First: Set ALL Diet_EX_ reactions to lower bound 0 (like useDiet does)
        for rxn in model.reactions:
            if rxn.id.startswith('Diet_EX_'):
                rxn.lower_bound = 0

        # Apply diet
        for _, row in diet_constraints.iterrows():
            rxn = row['rxn_id']
            if rxn in model.reactions:
                model.reactions.get_by_id(rxn).lower_bound = float(row['lower_bound'])
                if pd.notnull(row['upper_bound']):
                    model.reactions.get_by_id(rxn).upper_bound = float(row['upper_bound'])

        print(f"Processing {samp}: diet applied")
        
        # Constrain community biomass
        if 'communityBiomass' in model.reactions:
            model.reactions.communityBiomass.lower_bound = lower_bm
            model.reactions.communityBiomass.upper_bound = upper_bm

        for rxn in model.reactions:
            if rxn.id.startswith('UFEt_') or rxn.id.startswith('DUt_') or rxn.id.startswith('EX_'):
                rxn.upper_bound = 1e6

        # Include the human mets if not already included in the diet
        for met_id, bound in humanMets.items():
            rxn_id = f'Diet_EX_{met_id}[d]'
            if rxn_id not in diet_constraints['rxn_id'].values:  
                rxn = create_rxn(rxn_id, 'Bacteroides_sp_2_1_33B_biomass525', '', bound, 10000.)
                model.add_reactions([rxn])
            model.reactions.get_by_id(rxn.id).bounds = bound, 10000.

        # close demand and limit sink reactions
        for rxn in model.reactions:
            if '_DM_' in rxn.id:
                rxn.lower_bound = 0
            elif '_sink_' in rxn.id:
                rxn.lower_bound = -1 

        # Objective: EX_microbeBiomass[fe]
        model.objective = 'EX_microbeBiomass[fe]'
        model.optimize()
        print(f"Processing {samp}: model optimized")
        
        # Save the diet-adapted model
        save_dir = os.path.join(res_path, 'Diet')
        os.makedirs(save_dir, exist_ok=True)
        diet_model_path = os.path.join(save_dir, f"pymicrobiota_model_diet_{samp}.json")
        cobra.io.save_json_model(model, diet_model_path)

        print(f"Processing {samp}: starting fva")

        # Run FVA
        fecal_rxns = [r.id for r in model.reactions if '[fe]' in r.id]
        diet_rxns = [r.id for r in model.reactions if '[d]' in r.id]

        min_flux_fecal, max_flux_fecal = run_fva(model, fecal_rxns)
        print(f"Processing {samp}: starting diet fva")
        min_flux_diet, max_flux_diet = run_fva(model, diet_rxns)

        net_production[samp] = {}
        net_uptake[samp] = {}

        # model.exchanges contains only EX_rxns in fecal compartment

        # exchanges derived from exMets (all exchanged metabolites across all individual models) -> intersect it with rxns in this particular model
        exchanges = set(fecal_rxns).intersection(set(exchanges))

        # cut off very small values below solver sensitivity
        tol = 1e-07

        for rxn in exchanges:
            fecal = rxn
            diet = rxn.replace('EX_', 'Diet_EX_').replace('[fe]', '[d]')

            if abs(max_flux_fecal.get(fecal, 0)) < tol: max_flux_fecal.get(fecal, 0) == 0

            prod = min_flux_diet.get(diet, 0) + abs(max_flux_fecal.get(fecal, 0))
            uptk = max_flux_diet.get(diet, 0) + abs(min_flux_fecal.get(fecal, 0))
            print(rxn, min_flux_diet.get(diet, 0), max_flux_diet.get(diet, 0), min_flux_fecal.get(fecal, 0), max_flux_fecal.get(fecal, 0), prod, uptk)
            net_production[samp][rxn] = prod
            net_uptake[samp][rxn] = uptk

        print(f"Processing {samp}: completed")
        
    return exchanges, net_production, net_uptake

def collect_flux_profiles(samp_names, exchanges, net_production, net_uptake):
    """
    Convert net production/uptake dicts into Pandas DataFrames.

    Returns:
        net_secretion_df, net_uptake_df (pd.DataFrame)
    """
    prod_data = {samp: [net_production[samp].get(r, 0) for r in exchanges] for samp in samp_names}
    uptk_data = {samp: [net_uptake[samp].get(r, 0) for r in exchanges] for samp in samp_names}

    net_secretion_df = pd.DataFrame(prod_data, index=exchanges)
    net_uptake_df = pd.DataFrame(uptk_data, index=exchanges)

    return net_secretion_df, net_uptake_df
