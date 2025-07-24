"""
Microbiome Genome-Scale Modeling Toolbox (MiGEMox)

Optimized Python implementation of mgPipe community modeling pipeline.
Transforms individual AGORA species reconstructions into integrated community models
with diet and fecal compartments for host-microbiome interaction modeling.
"""

# Computational biology imports
import cobra
from cobra import Reaction, Metabolite
from cobra.io import load_matlab_model

# Numerical computing
import numpy as np
from scipy.sparse import csr_matrix, hstack
from scipy.io import savemat

# Data processing
import pandas as pd

# System utilities
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# Metabolite exchange bounds (mmol/gDW/h)
EXCHANGE_BOUNDS = (-1000.0, 10000.0) # Max uptake and secretion rates

# Transport reaction bounds (mmol/gDW/h) 
TRANSPORT_BOUNDS = (0.0, 10000.0) # Unidirectional transport

# Species inclusion threshold
ABUNDANCE_THRESHOLD = 1e-7

# Default Coupling Factor (Used in https://doi.org/10.4161/gmic.22370)
COUPLING_FACTOR = 400

def create_rxn(rxn_identifier: str, name: str, subsystem: str, bounds: tuple) -> cobra.Reaction:
    """
    Create a COBRA reaction with specified bounds and metadata.
    
    Args:
        rxn_identifier: Unique reaction ID
        name: Human-readable reaction name
        subsystem: Metabolic subsystem classification
        bounds: Lower and upper bounds of reaction (mmol/gDW/h)
        
    Returns:
        Configured COBRA reaction object
    """
    rxn = Reaction(rxn_identifier)
    rxn.name = name
    rxn.subsystem = subsystem
    lb, ub = bounds
    rxn.lower_bound = lb
    rxn.upper_bound = ub
    return rxn

def add_diet_fecal_compartments(model: cobra.Model) -> cobra.Model:
    """
    Add diet and fecal compartments to community model for host interaction.
    
    Biological System Modeled:
    - Diet compartment [d]: Nutrients from host dietary intake
    - Lumen compartment [u]: Shared microbial metabolite pool  
    - Fecal compartment [fe]: Metabolites excreted from host system
    
    Transport Chain: Diet[d] → DUt → Lumen[u] → UFEt → Fecal[fe] → EX_
    
    This creates the host-microbiome metabolite exchange interface essential
    for modeling dietary interventions and metabolite production.

    For every general metabolite in the lumen, 4 reactions will be added:
        (diet)
        EX_2omxyl[d]: 2omxyl[d] <=>
        DUt_2omxyl: 2omxyl[d] <=> 2omxyl[u]
        
        (fecal)
        UFEt_2omxyl: 2omxyl[u] <=> 2omxyl[fe]
        EX_2omxyl[fe]: 2omxyl[fe] <=>
    
    Args:
        model: Community model with species-tagged reactions
        
    Returns:
        Model with diet and fecal compartments and exchange and transport reactions
    """
    # Delete all EX_ reaction artifacts from the single cell models
    # E.g., EX_dad_2(e): dad_2[e] <=>, EX_thymd(e): thymd[e] <=>
    to_remove = [r for r in model.reactions if "_EX_" in r.id or "(e)" in r.id]
    model.remove_reactions(to_remove)

    # Create the diet and fecal compartments for reactions and metabolites
    # Get all of our general extracellular metabolites
    general_mets = []
    for reac in model.reactions:
        if "IEX" in reac.id:
            iex_reac = model.reactions.get_by_id(reac.id)
            # Pick only general (unlabeled) metabolites on the LHS
            for met in iex_reac.reactants:
                if "[u]" in met.id:
                    general_mets.append(met.id)
    general_mets = set(general_mets)

    # Create diet and fecal compartments, with new transport and exchange reactions
    existing_mets = {m.id for m in model.metabolites}
    existing_rxns = {r.id for r in model.reactions}
    
    for lumen_met in general_mets:
        base_name = lumen_met.split('[')[0]  # Remove [u] suffix

        # EX_2omxyl[d]: 2omxyl[d] <=>
        _add_exchange_reaction(base_name, existing_mets, model, EXCHANGE_BOUNDS, "d", "diet")
        # DUt_4hbz: 4hbz[d] --> 4hbz[u]
        _add_transport_reaction(f'DUt_{base_name}', existing_rxns, model, f'{base_name}[d]', lumen_met, TRANSPORT_BOUNDS, "diet to lumen")
        # EX_4abut[fe]: 4abut[fe] <=>
        _add_exchange_reaction(base_name, existing_mets, model, EXCHANGE_BOUNDS, "fe", "fecal")
        # UFEt_arabinoxyl: arabinoxyl[u] --> arabinoxyl[fe]
        _add_transport_reaction(f'UFEt_{base_name}', existing_rxns, model, lumen_met, f'{base_name}[fe]', TRANSPORT_BOUNDS, "lumen to fecal")

    return model

def _add_exchange_reaction(base_name, existing_met_ids, model, bounds, compartment, label):
    met_id = f'{base_name}[{compartment}]'
    if met_id not in existing_met_ids:
        reac_id = "EX_" + met_id
        reaction = create_rxn(reac_id, f"{met_id} {label} exchange", ' ', bounds)
        model.add_reactions([reaction])
        model.add_metabolites([Metabolite(met_id, compartment=compartment)])
        reaction = model.reactions.get_by_id(reac_id)
        reaction.add_metabolites({model.metabolites.get_by_id(met_id): -1})

def _add_transport_reaction(rxn_id, existing_rxn_ids, model, reactant_id, product_id, bounds, label):
    if rxn_id not in existing_rxn_ids:
        reaction = create_rxn(rxn_id, f"{rxn_id} {label} transport", ' ', bounds)
        model.add_reactions([reaction])
        reaction.reaction = f"{reactant_id} --> {product_id}"
        reaction.bounds = bounds

def com_biomass(model: cobra.Model, abun_path: str, sample_com: str):
    """
    Create weighted community biomass reaction based on species abundances.
    
    Biological Equation: Community_Biomass = Σ(abundance_i × species_biomass_i)
    
    This represents the total microbial biomass production weighted by each
    species' relative abundance in the sample, creating a community-level
    growth objective that reflects the natural composition.
    
    Args:
        model: Community model with individual species biomass reactions
        abundance_path: Path to species abundance CSV file
        sample_name: Column name in abundance file for this sample
        
    Returns:
        Model with community biomass reaction and transport to fecal compartment
    """

    # Deleting all previous community biomass equations
    biomass_reactions = [r for r in model.reactions if "Biomass" in r.id]
    model.remove_reactions(biomass_reactions)

    # Load abundance data and filter by threshold
    abun_df = pd.read_csv(abun_path)
    abun_df = abun_df[abun_df[sample_com] > ABUNDANCE_THRESHOLD]

    # Creating the community biomass reaction
    reaction = create_rxn("communityBiomass", "communityBiomass", ' ', (0., 10000.))
    model.add_reactions([reaction])
    community_biomass = model.reactions.communityBiomass

    # Build abundance-weighted biomass stoichiometry
    biomass_stoichiometry = {}
    
    for _, row in abun_df.iterrows():
        species_name = row["X"]
        abundance = float(row[sample_com])
        biomass_met_id = f"{species_name}_biomass[c]"
        if biomass_met_id in model.metabolites:
            biomass_stoichiometry[biomass_met_id] = -abundance
        else:
            print(f"⚠️ Biomass metabolite missing in model: {biomass_met_id}")
    
    community_biomass.add_metabolites(metabolites_to_add=biomass_stoichiometry, combine=True)

    # Adding the microbeBiomass metabolite
    model.add_metabolites([Metabolite("microbeBiomass[u]", formula=" ", \
                                      name="product of community biomass", compartment="u"),])
    community_biomass.add_metabolites({model.metabolites.get_by_id("microbeBiomass[u]"): 1})

    # Adding the exchange reaction compartment
    reac_name = "EX_microbeBiomass[fe]"
    reaction = create_rxn(reac_name, reac_name, ' ', (-10000., 10000.))
    model.add_reactions([reaction])
    model.add_metabolites([Metabolite("microbeBiomass[fe]", formula=" ", \
                                      name="product of community biomass", compartment="fe"),])
    new_fe_react = model.reactions.get_by_id("EX_microbeBiomass[fe]")
    new_fe_react.add_metabolites({model.metabolites.get_by_id("microbeBiomass[fe]"): -1})

    # Adding the UFEt reaction
    reaction = create_rxn("UFEt_microbeBiomass", "UFEt_microbeBiomass", ' ', TRANSPORT_BOUNDS)
    model.add_reactions([reaction])
    reaction.reaction = "microbeBiomass[u] --> microbeBiomass[fe]"
    reaction.bounds = TRANSPORT_BOUNDS

    return model

def tag_metabolite(met: cobra.Metabolite, species_name: str, compartment: str):
    '''
    Helper function for species_to_community for tagging metabolites
    '''
    met.compartment = compartment
    no_c_name = met.id.replace(f"[{compartment}]", "")
    met.id = f'{species_name}_{no_c_name}[{compartment}]'

def species_to_community(model: cobra.Model, species_model_name: str):
    """
    Takes a single cell AGORA GEM and changes its reaction and metabolite formatting so it 
    can be added to a community model in the Com_py pipeline.
  
        Tags everything intracellular and intracellular to extracellular with the species name:
            (intracellular)
            tag[c] -> tag[c]

            (transport)
            tag[c] -> tag[u]

            (IEX reactions)
            tagged[e] -> general[e]
  
    INPUTS:
        model: a .mat file of an AGORA single cell model
        species_model_name: the species name to be tagged (extracted in the Com_py pipeline)
  
    OUTPUTS:
        model: updated model with tagged reactions and metabolites
    """
    # Tagging all reactions and metabolites in the cell with species name
    # For each species model, iterate through its reactions and add the species tag

    # Extracting the species name from the species model name
    short_species_name = os.path.splitext(os.path.basename(species_model_name))[0]

    # Step 1: Remove all exchange reactions except for the biomass reaction
    ex_rxns = [rxn for rxn in model.reactions if "EX_" in rxn.id and "biomass" not in rxn.id]
    model.remove_reactions(ex_rxns)

    # Step 2: Tag metabolites in intra- and extracellular compartments of model
    for rxn in model.reactions:
        # Change the intracellualr reaction from [c] --> [c]
        if "[e]" in rxn.reaction or "[c]" in rxn.reaction:
            rxn.id = f'{short_species_name}_{rxn.id}'
            # tag each metabolite in rxn, if not already tagged
            for met in rxn.metabolites:
                if ("[c]" in met.id or "[p]" in met.id) and short_species_name not in met.id:
                    compartment = 'c' if '[c]' in met.id else 'p'
                    tag_metabolite(met, short_species_name, compartment)
                elif "[e]" in met.id and short_species_name not in met.id:
                    met.compartment = "u"
                    no_c_name = met.id.replace("[e]", "")
                    met.id = f'{short_species_name}_{no_c_name}[u]'

    # Step 3: Create inter-species metabolite exchange
    model = _create_inter_species_exchange(model, short_species_name)

    # Step 4: Ensure all components are properly tagged
    model = _finalize_species_tagging(model, short_species_name)

    return model

def _create_inter_species_exchange(model: cobra.Model, species_name: str) -> cobra.Model:
    """
    Create IEX reactions for species-specific ↔ general metabolite exchange.
    
    Biological Rationale: Allows species to contribute/consume shared metabolites
    in community lumen while maintaining species-specific uptake kinetics.
    """
    species_lumen_metabolites = [
        met for met in model.metabolites 
        if "[u]" in met.id and species_name in met.id
    ]
    
    for species_met in species_lumen_metabolites:
        general_met_id = species_met.id.replace(f"{species_name}_", "")
        
        # Create general metabolite if it doesn't exist
        if general_met_id not in [m.id for m in model.metabolites]:
            general_met = Metabolite(
                general_met_id, 
                compartment="u", 
                name=general_met_id.split("[")[0]
            )
            model.add_metabolites([general_met])
        
        # Create IEX reaction: general_met <=> species_met
        iex_rxn_id = f"{species_name}_IEX_{general_met_id}tr"
        iex_rxn = create_rxn(iex_rxn_id, f"{species_name}_IEX", " ", (-1000.0, 1000.0))
        model.add_reactions([iex_rxn])
        iex_rxn.reaction = f"{general_met_id} <=> {species_met.id}"
        iex_rxn.bounds = (-1000.0, 1000.0)
    
    return model

def _finalize_species_tagging(model: cobra.Model, species_name: str) -> cobra.Model:
    """Ensure all reactions and metabolites are properly tagged with species name."""
    # Tag any remaining untagged reactions
    for rxn in model.reactions:
        if not rxn.id.startswith(species_name):
            rxn.id = f"{species_name}_{rxn.id}"
    
    # Tag any remaining untagged [c] and [p] metabolites
    for met in model.metabolites:
        if ("[c]" in met.id or "[p]" in met.id) and species_name not in met.id:
            compartment = 'c' if '[c]' in met.id else 'p'
            tag_metabolite(met, species_name, compartment)
    
    return model

def prune_zero_abundance_species(model: cobra.Model, zero_abundance_species: list[str]) -> cobra.Model:
    """
    Remove all reactions and metabolites from species below abundance threshold.
    Biological Rationale: Species below detection limit (typically 0.01% relative
    abundance) don't contribute meaningful metabolic flux to community phenotype.

    Args:
        model: Community model containing all species
        zero_abundance_species: Species names below abundance threshold
        
    Returns:
        Pruned model containing only detected species
    """
    print('Pruning metabolites and Reactions from Zero-Abundance Species in Sample')
    zero_prefixes = {f"{species}_" for species in zero_abundance_species}
    metabolites_to_remove = [
        met for met in model.metabolites 
        if any(met.id.startswith(prefix) for prefix in zero_prefixes)
    ]
    # Destructive removal also removes associated reactions automatically
    model.remove_metabolites(metabolites_to_remove, destructive=True)
    print(f"Pruned {len(metabolites_to_remove)} Metabolites")
    return model

def build_global_model(abundance_df: pd.DataFrame, mod_dir: str) -> cobra.Model:
    """
    Loads all species found in the abundance table and builds a unified, unpruned community model.
    Species are combined into a single COBRA model, tagged and merged. Cleans the community as well
    by adding the diet and fecal compartments

    Parameters:
        abundance_df: species x samples abundance dataframe.
        mod_dir: path to folder containing AGORA .mat files.

    Returns:
        global_model: a unified COBRApy model containing all species.
    """

    print("Building global community model".center(40, '*'))
    all_species = abundance_df.index.tolist()
    first_path = os.path.join(mod_dir, all_species[0] + ".mat")
    print(f"Added first species model: {all_species[0]}".center(40, '*'))
    first_model = cobra.io.load_matlab_model(first_path)
    global_model = species_to_community(first_model, species_model_name=first_path)

    for species in all_species[1:]:
        species_path = os.path.join(mod_dir, species + ".mat")
        model = cobra.io.load_matlab_model(species_path)
        tagged_model = species_to_community(model, species_path)

        # Avoid duplicate reaction IDs
        existing_rxns = {r.id for r in global_model.reactions}
        new_rxns = [r for r in tagged_model.reactions if r.id not in existing_rxns]
        global_model.add_reactions(new_rxns)

    print("Finished adding GEM reconstructions to community".center(40, '*'))

    print("Adding diet and fecal compartments".center(40, '*'))
    clean_model = add_diet_fecal_compartments(model=global_model)
    print("Done adding diet and fecal compartments".center(40, '*'))

    global_C, global_d, global_dsense, global_ctrs = build_global_coupling_constraints(clean_model, all_species)

    return clean_model, global_C, global_d, global_dsense, global_ctrs

def build_sample_model(sample_name: str, global_model: cobra.Model, abundance_df: pd.DataFrame, 
                       abun_path: str, out_dir: str, diet_path: str = None,
                       global_C=None, global_d=None, global_dsense=None, global_ctrs=None):
    """
    Takes a deep copy of the global model and builds the sample-specific model:
        - prunes zero-abundance species
        - adds diet constraints
        - adds community biomass
        - saves as .json

    Parameters:
        sample_name: column name in abundance_df
        global_model: unpruned community model
        abundance_df: pandas dataframe of abundances
        abun_path: path to abundance CSV (needed by com_biomass)
        out_dir: directory to save the output model
        diet_path: path to the diet file (optional)

    Returns:
        Path to the saved model
    """
    model = global_model.copy()
    model.name = sample_name
    sample_abun = abundance_df[sample_name]

    # Prune zero-abundance species from the model
    zero_species = [sp for sp in sample_abun.index if sample_abun[sp] < 1e-7]
    present_species = [sp for sp in sample_abun.index if sample_abun[sp] >= 1e-7]

    model = prune_zero_abundance_species(model, zero_abundance_species=zero_species)

    # Add a community biomass reaction to the model
    print("Adding community biomass reaction".center(40, '*'))
    model = com_biomass(model=model, abun_path=abun_path, sample_com=sample_name)

    # Prune coupling constraints from the global model (C, dsense, d, ctrs)
    sample_C, sample_d, sample_dsense, sample_ctrs = prune_coupling_constraints_by_species(
        global_model, global_C, global_d, global_dsense, global_ctrs, present_species, model
    )

    # Ensuring the reversablity fits all compartments
    for reac in [r for r in model.reactions if "DUt" in r.id or "UFEt" in r.id]:
        reac.lower_bound = 0.

    # Setting the new community biomass as the objective
    # Setting EX_microbeBiomass[fe] as objective to match MATLAB mgPipe
    model.objective = "EX_microbeBiomass[fe]"

    os.makedirs(out_dir, exist_ok=True)
    model_dict = make_mg_pipe_model_dict(
        model, C=sample_C, d=sample_d, dsense=sample_dsense, ctrs=sample_ctrs
    )

    save_path = os.path.join(out_dir, f"microbiota_model_samp_{sample_name}.mat")
    savemat(save_path, {'model': model_dict}, do_compression=True, oned_as='column')

    return save_path

def build_global_coupling_constraints(model: cobra.Model, species_list: list[str], coupling_factor: float=400):
    """
    Build sparse coupling constraint matrix for species biomass relationships.
    
    Biological Constraint: v_reaction ≤ coupling_factor × v_biomass
    
    Rationale: Individual reaction rates cannot exceed species biomass production
    by more than the coupling factor. Prevents unrealistic flux distributions
    where species have high metabolic activity but low biomass.
    
    Args:
        model: Community model with all species reactions
        species_list: List of species identifiers  
        coupling_factor: Biomass coupling strength (default from config)
        
    Returns:
        Coupling matrix (C), bounds (d), constraint sense (dsense), names (ctrs)
    """
    rxn_id_to_index = {r.id: i for i, r in enumerate(model.reactions)}

    all_constraints = []
    all_d = []
    all_dsense = []
    all_ctrs = []

    for species in species_list:
        # Find species reactions and biomass reaction
        species_rxns = [r for r in model.reactions if r.id.startswith(species + '_')]
        biomass_rxns = [r for r in species_rxns if 'biomass' in r.id.lower()]

        if not biomass_rxns:
            continue

        biomass_rxn = biomass_rxns[0]
        biomass_idx = rxn_id_to_index[biomass_rxns[0].id]

        for rxn in species_rxns:
            if rxn.id == biomass_rxn.id:
                continue  # Don't couple biomass to itself

            rxn_idx = rxn_id_to_index[rxn.id]

            # Create constraint: v_rxn - 400*v_biomass <= 0
            constraint_row = np.zeros(len(model.reactions))
            constraint_row[rxn_idx] = 1.0  # coefficient for v_rxn
            constraint_row[biomass_idx] = -coupling_factor  # coefficient for v_biomass

            all_constraints.append(constraint_row)
            all_d.append(0.0)
            all_dsense.append('L')  # <= constraint
            all_ctrs.append(f"slack_{rxn.id}")

            # Also add reverse constraint: v_rxn + 400*v_biomass >= 0 (for reversible reactions)
            if rxn.lower_bound < 0:
                constraint_row_rev = np.zeros(len(model.reactions))
                constraint_row_rev[rxn_idx] = 1.0
                constraint_row_rev[biomass_idx] = coupling_factor

                all_constraints.append(constraint_row_rev)
                all_d.append(0.0)
                all_dsense.append('G')
                all_ctrs.append(f"slack_{rxn.id}_R")

    if all_constraints:
        C = csr_matrix(np.vstack(all_constraints))
        d = np.array(all_d).reshape(-1, 1)
        dsense = np.array(all_dsense, dtype='<U1')
        ctrs = np.array(all_ctrs, dtype=object)
    else:
        C = csr_matrix((0, len(model.reactions)))
        d = np.zeros((0, 1))
        dsense = np.array([], dtype='<U1')
        ctrs = np.array([], dtype=object)

    return C, d, dsense, ctrs

def prune_coupling_constraints_by_species(global_model, global_C, global_d, global_dsense, global_ctrs, 
                                        present_species, model):
    """
    Prune coupling constraints to only include those for species present in the sample.
    This mimics MATLAB's approach of removing coupling matrices for zero-abundance species.
    """

    present_species_set = set(present_species)
    slack_prefix = "slack_"
    keep_rows = []

    for i, ctr_name in enumerate(tqdm(global_ctrs)):
        # Extract species name from constraint name (e.g., "slack_Bacteroides_sp_2_1_33B_IEX_12ppd_S[u]tr")
        if ctr_name.startswith(slack_prefix):
            species_part = ctr_name[len(slack_prefix):]
            for species in present_species_set:
                if species_part.startswith(species):
                    keep_rows.append(i)
                    break
    
    if keep_rows:
          # Remap columns to match sample-specific model
        sample_rxn_ids = [r.id for r in model.reactions]
        global_rxn_ids = [r.id for r in global_model.reactions]

        global_rxn_idx_map = {rid: i for i, rid in enumerate(global_rxn_ids)}

        cols = []
        for rid in sample_rxn_ids:
            idx = global_rxn_idx_map.get(rid)
            if idx is not None:
                cols.append(global_C[keep_rows, idx])
            else:
                cols.append(csr_matrix((len(keep_rows), 1)))
        pruned_C = hstack(cols)
        
        pruned_d = global_d[keep_rows, :]
        pruned_dsense = global_dsense[keep_rows]
        pruned_ctrs = global_ctrs[keep_rows]
    else:
        pruned_C = csr_matrix((0, len(model.reactions)))
        pruned_d = np.zeros((0, 1))
        pruned_dsense = np.array([], dtype='<U1')
        pruned_ctrs = np.array([], dtype=object)
    
    return pruned_C, pruned_d, pruned_dsense, pruned_ctrs

def make_mg_pipe_model_dict(model, C=None, d=None, dsense=None, ctrs=None):
    """
    Enhanced version that properly handles all MATLAB-compatible fields
    """
    num_rxns = len(model.reactions)
    num_mets = len(model.metabolites)

    # Objective vector
    c = np.zeros((num_rxns, 1))
    obj_rxn_id = str(model.objective.expression).split('*')[1].split('-')[0].strip()
    for i, rxn in enumerate(model.reactions):
        if rxn.id == obj_rxn_id:
            c[i, 0] = 1.0
            break

    # S matrix
    S = cobra.util.create_stoichiometric_matrix(model)

    # Bounds
    lb = np.array([rxn.lower_bound for rxn in model.reactions], dtype=np.float64).reshape(-1, 1)
    ub = np.array([rxn.upper_bound for rxn in model.reactions], dtype=np.float64).reshape(-1, 1)

    # Met and rxn info
    rxns = np.array([[rxn.id] for rxn in model.reactions], dtype=object)
    rxnNames = np.array([[rxn.name] for rxn in model.reactions], dtype=object)
    metNames = np.array([[met.name] for met in model.metabolites], dtype=object)
    mets = np.array([[met.id] for met in model.metabolites], dtype=object)

    # Sense
    csense = np.array(['E'] * num_mets, dtype='U1')

    # Default coupling matrices
    if C is None: C = csr_matrix((0, num_rxns))
    if d is None: d = np.zeros((0, 1))
    if dsense is None: dsense = np.array([], dtype='<U1')
    if ctrs is None: ctrs = np.array([], dtype=object).reshape(-1, 1)

    C = csr_matrix((0, num_rxns)) if C is None else C
    d = np.zeros((0, 1)) if d is None else d
    dsense = np.array([], dtype='<U1') if dsense is None else dsense
    ctrs = np.array([], dtype=object).reshape(-1, 1) if ctrs is None else ctrs.reshape(-1, 1)

    # Model name
    model_name = np.array([model.name], dtype=object)
    osenseStr = np.array(['max'], dtype='U3')

    return {
        'rxns': rxns,
        'rxnNames': rxnNames,
        'mets': mets,
        'metNames': metNames,
        'S': S,
        'b': np.zeros((num_mets, 1)),
        'c': c,
        'lb': lb,
        'ub': ub,
        'metChEBIID': np.array([
            m.annotation['chebi'][0].replace('CHEBI:', '') if 'chebi' in m.annotation and isinstance(m.annotation['chebi'], list)
            else m.annotation['chebi'].replace('CHEBI:', '') if 'chebi' in m.annotation and isinstance(m.annotation['chebi'], str)
            else ''
            for m in model.metabolites
        ], dtype=object).reshape(-1, 1),
        'metCharges': np.array([m.charge if getattr(m, 'charge', None) is not None else np.nan for m in model.metabolites]).reshape(-1, 1),
        'metFormulas': np.array([m.formula if getattr(m, 'formula', None) is not None else np.nan for m in model.metabolites]).reshape(-1, 1),
        'rules': np.array([getattr(r, 'gene_reaction_rule', '') for r in model.reactions], dtype=object).reshape(-1, 1),
        'subSystems': np.array([getattr(r, 'subsystem', '') for r in model.reactions], dtype=object).reshape(-1, 1),
        'osenseStr': osenseStr,
        'csense': csense,
        'C': C,
        'd': d,
        'dsense': dsense,
        'ctrs': ctrs,
        'name': model_name
    }

def migemox(abun_filepath, mod_filepath, out_filepath, diet_filepath=None, workers=1):
    """
    Inspired by MICOM community building and mgpipe.m code.
    Main pipeline which inputs the GEMs data and accesses the different functions.

    INPUTS:
        abun_path: path to the species abundance .csv file
            Formatting for the species abundance:
                The columns should have the names of the .mat files of the species you want to load
                See file normCoverage_smaller.csv for template example   
        modpath: path to folder with all AGORA models 
            E.g. "~/data_input/AGORA201/'"
        respath: path where the community models will be output (defaults to same folder as )
            E.g. "~/results/"
        dietpath: path to the AGORA compatible diet (for community model) .csv file
            E.g. "~/data_input/AverageEU_diet_fluxes.csv"
  
    OUTPUTS:
        All sample community models to a specified local folder
    """

    print("Starting MiGeMox pipeline".center(40, '*'))
    print("Reading abundance file".center(40, '*'))

    sample_info = pd.read_csv(abun_filepath)
    sample_info.rename(columns={list(sample_info)[0]:"species"}, inplace=True)
    sample_info.set_index("species", inplace=True)

    global_model, global_C, global_d, global_dsense, global_ctrs = build_global_model(sample_info, mod_filepath)
    samples = sample_info.columns.tolist()

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(build_sample_model, s, global_model, sample_info, abun_filepath, out_filepath, 
                                   diet_filepath, global_C, global_d, global_dsense, global_ctrs)
                   for s in samples]
        for f in tqdm(futures, desc='Building sample models'):
            f.result()