import os
import csv
from tqdm import tqdm
import pandas as pd
import cobra.io
from cobra import Reaction, Metabolite
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.sparse import csr_matrix, hstack
from scipy.io import savemat

# Model configuration constants
EXCHANGE_LOWER_BOUND = -1000.0
EXCHANGE_UPPER_BOUND = 10000.0
TRANSPORT_LOWER_BOUND = 0.0
TRANSPORT_UPPER_BOUND = 10000.0
ABUNDANCE_THRESHOLD = 1e-7

def print_out(string, fixed_length=20, padding_str="-"):
    """
    Helper function that prints a string with a fixed length and padding on both sides.
    
    Does not support f strings.

    INPUTS:
        string: string to be printed
        fixed_length: length of the printed string
        padding_str: string to be used for padding

    OUTPUTS:
        None
    """
    padding = padding_str * 10
    formatted_string = "{:-^{}}".format(string, fixed_length - len(padding))
    padded_string = padding + formatted_string + padding
    print(padded_string)

def create_rxn(rxn_identifier, name, subsystem, lb, ub) -> cobra.Reaction:
    rxn = Reaction(rxn_identifier)
    rxn.name = name
    rxn.subsystem = subsystem
    rxn.lower_bound = lb
    rxn.upper_bound = ub

    return rxn

def clean_community(model):
    """
    Takes a combined community model and creates fecal and diet transport and exchange reactions.
    
    This function adds 4 types of reactions for every general metabolite in the lumen:
        (diet)
        EX_2omxyl[d]: 2omxyl[d] <=>
        DUt_2omxyl: 2omxyl[d] <=> 2omxyl[u]
        
        (fecal)
        UFEt_2omxyl: 2omxyl[u] <=> 2omxyl[fe]
        EX_2omxyl[fe]: 2omxyl[fe] <=>
        
    INPUTS:
        model: a .mat file of an AGORA single cell model
  
    OUTPUTS:
        model: updated model with fecal and diet compartments
    """
    # Delete all EX_ reaction artifacts from the single cell models
    # E.g., EX_dad_2(e): dad_2[e] <=>, EX_thymd(e): thymd[e] <=>
    to_remove = [r for r in model.reactions if "_EX_" in r.id or "(e)" in r.id]
    model.remove_reactions(to_remove)

    existing_met_ids = [met.id for met in model.metabolites]
    existing_rxn_ids = [rxn.id for rxn in model.reactions]

    # Create the diet and fecal compartments for reactions and metabolites
    # Get all of our general extracellular metabolites
    gen_mets = []
    for reac in model.reactions:
        if "IEX" in reac.id:
            iex_reac = model.reactions.get_by_id(reac.id)
            # Pick only general (unlabeled) metabolites on the LHS
            for met in iex_reac.reactants:
                if "[u]" in met.id:
                    gen_mets.append(met.id)
    gen_mets = set(gen_mets)

    # Creating the diet and fecal compartments
    for met_name in gen_mets:

        clean_met_name = met_name.split('[')[0] # strip [u] from end
        d_met_name = clean_met_name + "[d]" # met[d]
        dut_name = "DUt_" + clean_met_name # Rxn name = DUt_met
        fe_met_name = clean_met_name +  '[fe]' # met[fe]
        ufet_name = 'UFEt_' + clean_met_name # Rxn name = UFEt_met

        # Creating the [d] exchange reactions
        # E.g., EX_2omxyl[d]: 2omxyl[d] <=>
        # Ensure there are no duplicates if user runs function more than once
        if d_met_name not in existing_met_ids:
            reac_name = "EX_" + d_met_name
            reaction = create_rxn(reac_name, d_met_name + 'diet exchange', ' ', EXCHANGE_LOWER_BOUND, EXCHANGE_UPPER_BOUND)
            model.add_reactions([reaction])
            model.add_metabolites([Metabolite(d_met_name, compartment="d")])
            new_dietreact = model.reactions.get_by_id(reac_name)
            new_dietreact.add_metabolites({model.metabolites.get_by_id(d_met_name): -1})

        # Creating the [d] transport reactions to lumen
        # E.g., DUt_4hbz: 4hbz[d] --> 4hbz[u]
        # Ensure there are no duplicates if user runs function more than once
        if dut_name not in existing_rxn_ids:
            dut_formula = f"{d_met_name} --> {met_name}"
            reaction = create_rxn(dut_name, dut_name + 'diet to lumen', ' ', TRANSPORT_LOWER_BOUND, TRANSPORT_UPPER_BOUND)
            model.add_reactions([reaction])

            # Adding the correct d --> u formula to the reaction
            reaction.reaction = dut_formula
            reaction.bounds = (0., 10000.)

        # Creating the fe exchange reactions
        # E.g., EX_4abut[fe]: 4abut[fe] <=>
        # Ensure there are no duplicates if user runs function more than once
        if fe_met_name not in existing_met_ids:
            reac_name = "EX_" + fe_met_name
            reaction = create_rxn(reac_name, fe_met_name + "fecal exchange", ' ', EXCHANGE_LOWER_BOUND, EXCHANGE_UPPER_BOUND)
            model.add_reactions([reaction])
            model.add_metabolites([Metabolite(fe_met_name, compartment="fe")])
            new_fe_react = model.reactions.get_by_id(reac_name)
            new_fe_react.add_metabolites({model.metabolites.get_by_id(fe_met_name): -1})

        # Creating the [fe] transport reactions to lumen
        # E.g., UFEt_arabinoxyl: arabinoxyl[u] --> arabinoxyl[fe]
        # Ensure there are no duplicates if user runs function more than once
        if ufet_name not in existing_rxn_ids:
            ufet_formula = f"{met_name} --> {fe_met_name}"
            reaction = create_rxn(ufet_name, ufet_name + "diet to lumen", ' ', TRANSPORT_LOWER_BOUND, TRANSPORT_UPPER_BOUND)
            model.add_reactions([reaction])

            # Adding the corrected d --> u formula to the reaction
            reaction.reaction = ufet_formula
            reaction.bounds = (0., TRANSPORT_UPPER_BOUND)

    return model

def com_biomass(model, abun_path, sample_com):
    """
    Takes a combined community model and adds a community biomass formula to the model.
    
    INPUTS:
        model: a .mat file of an AGORA single cell model
        abun_path: path to the species abundance .csv file
        sample_com: the sample name string (internal to the com_py pipeline)
  
    OUTPUTS:
        model: updated model with community biomass equation
    """

    # Deleting all previous community biomass equations
    biomass_reactions = [r for r in model.reactions if "Biomass" in r.id]
    model.remove_reactions(biomass_reactions)

    # Extracting biomass metabolites from the different single cell models
    biomass_mets_list = []
    for mets in model.metabolites:
        if "biomass" in mets.id:
            biomass_mets_list.append(mets.id)

    # Sort species alphabetically
    biomass_mets_list = sorted(biomass_mets_list)

    # Reading in the abundance file, sort species alphabetically,
    # and remove species with abundances < 1e-7
    norm_abund = pd.read_csv(abun_path)
    norm_abund = norm_abund[norm_abund[sample_com] > ABUNDANCE_THRESHOLD]

    # Creating the community biomass reaction
    reaction = create_rxn("communityBiomass", "communityBiomass", ' ', 0., 10000.)
    # reaction.add_metabolites(com_biomass_dic)
    model.add_reactions([reaction])
    community_biomass = model.reactions.communityBiomass

    # Initialize biomass contribution dictionary
    com_biomass_dict = {}

    # Build biomass ID from species name
    for _, row in norm_abund.iterrows():
        species = row["X"]
        abundance = row[sample_com]
        biomass_id = f"{species}_biomass[c]"
        if biomass_id in model.metabolites:
            com_biomass_dict[biomass_id] = -float(abundance)
        else:
            print(f"⚠️ Biomass metabolite missing in model: {biomass_id}")
    
    community_biomass.add_metabolites(metabolites_to_add=com_biomass_dict, combine=True)

    # Adding the microbeBiomass metabolite
    model.add_metabolites([Metabolite("microbeBiomass[u]", formula=" ", \
                                      name="product of community biomass", compartment="u"),])
    community_biomass.add_metabolites({model.metabolites.get_by_id("microbeBiomass[u]"): 1})

    # Adding the exchange reaction compartment
    reac_name = "EX_microbeBiomass[fe]"
    reaction = create_rxn(reac_name, reac_name, ' ', -10000., EXCHANGE_UPPER_BOUND)
    model.add_reactions([reaction])

    # Adding a microbeBiomass [fe] metabolite
    model.add_metabolites([Metabolite("microbeBiomass[fe]", formula=" ", \
                                      name="product of community biomass", compartment="fe"),])

    # Adding the exchange reaction
    new_fe_react = model.reactions.get_by_id("EX_microbeBiomass[fe]")
    new_fe_react.add_metabolites({model.metabolites.get_by_id("microbeBiomass[fe]"): -1})

    # Adding the UFEt reaction
    ufet_formula = "microbeBiomass[u] --> microbeBiomass[fe]"
    reaction = create_rxn("UFEt_microbeBiomass", "UFEt_microbeBiomass", ' ', TRANSPORT_LOWER_BOUND, TRANSPORT_UPPER_BOUND)
    model.add_reactions([reaction])

    # Adding the correct d --> u formula to the reaction
    reaction.reaction = ufet_formula
    reaction.bounds = (0., TRANSPORT_UPPER_BOUND)

    return model

def tag_metabolite(model, metabolite_name, species_name, compartment: str):
    '''
    Helper function for species_to_community for tagging metabolites
    '''
    model.metabolites.get_by_id(metabolite_name).compartment = compartment
    no_c_name = metabolite_name.replace(f"[{compartment}]", "")
    updated_m_name = species_name + "_" + no_c_name + f"[{compartment}]"
    model.metabolites.get_by_id(metabolite_name).id = updated_m_name
    model.metabolites.get_by_id(updated_m_name).compartment = compartment

def species_to_community(model, species_model_name):
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

    for rxn in model.reactions:
        # Removing all exchange reactions except for the biomass reaction
        if "EX_" in rxn.id and "biomass" not in rxn.id:
            model.remove_reactions(model.reactions.get_by_id(rxn.id))

        # Change the intracellualr reaction from [c] --> [c]
        elif "[e]" not in rxn.reaction and "[c]" in rxn.reaction:
                reaction_bounds = rxn.bounds
                reaction_name = rxn.id
                updated_r_name =  short_species_name + "_" + reaction_name
                model.reactions.get_by_id(reaction_name).id = updated_r_name

                # For each metabolite in that reaction, if that met already has the
                # name tag skip it, if not add it
                for met in rxn.metabolites:
                    if short_species_name not in str(met.id):
                        if "[c]" in str(met.id):
                            tag_metabolite(model, met.id, short_species_name, 'c')
                        if "[p]" in str(met.id):
                            tag_metabolite(model, met.id, short_species_name, 'p')

                # Get the reversiblity to carry over into the new models
                model.reactions.get_by_id(updated_r_name).bounds = reaction_bounds

        # Working on the extracellular reactions from [c] --> [e]
        elif "[e]" in rxn.reaction:
            reaction_name = rxn.id
            updated_r_name =  short_species_name + "_" + reaction_name
            model.reactions.get_by_id(reaction_name).id = updated_r_name

            # For each metabolite in that reaction, if it is the cellular metabolite, tag it;
            # if it is the extracellular metabolite, tag it
            for met in rxn.metabolites:
                if "[c]" in met.id:
                    if short_species_name not in str(met.id):
                        tag_metabolite(model, met.id, short_species_name, 'c')
                elif "[e]" in met.id:
                    if short_species_name not in str(met.id):
                        metabolite_name = met.id
                        model.metabolites.get_by_id(metabolite_name).compartment = "u"
                        no_c_name = metabolite_name.replace("[e]", "")
                        updated_m_name = short_species_name + "_" + no_c_name + "[u]"
                        model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                        model.metabolites.get_by_id(updated_m_name).compartment = "u"
                elif "[p]" in met.id:
                    if short_species_name not in str(met.id):
                        tag_metabolite(model, met.id, short_species_name, 'p')

    # Adding the tagged intracellular reactions for the extracellular metabolites to the model
    # For all extracellular species specific metabolites, add an exchange reaction
    for met in model.metabolites:
        if "[u]" in met.id:
            if short_species_name in met.id:
                replacing_specname = short_species_name + "_"
                general_name = met.id.replace(replacing_specname, "")
                iex_formula = f"{general_name} <=> {met.id}"
                iex_reaction_name = short_species_name + "_IEX_" + general_name + "tr"

                if general_name not in [m.id for m in model.metabolites]:
                    model.add_metabolites([
                        Metabolite(general_name, compartment="u", name=general_name.split("[")[0])
                    ])

                reaction = create_rxn(iex_reaction_name, short_species_name + "_IEX", ' ', -1000., 1000.)
                model.add_reactions([reaction])
                reaction.reaction = iex_formula
                reaction.bounds = (-1000., 1000.)

    # Ensure each reaction has the species tag on it
    for rxn in model.reactions:
        if short_species_name not in rxn.id:
            reaction_name = rxn.id
            updated_r_name =  short_species_name + "_" + reaction_name
            model.reactions.get_by_id(reaction_name).id = updated_r_name

    # Ensure each [c] and [p] compartment metabolite has the species tag on it
    for met in model.metabolites:
        if "[c]" in met.id:
            if short_species_name not in met.id:
                tag_metabolite(model, met.id, short_species_name, 'c')
        elif "[p]" in met.id:
            if short_species_name not in met.id:
                tag_metabolite(model, met.id, short_species_name, 'p')

    # Woring on individual biomass reactions
    # Each species will have a [c] --> [c] reaction that will then be added to
    # the community biomass
    for rxn in model.reactions:
        if short_species_name not in rxn.id:
            reaction_name = rxn.id
            model.reactions.get_by_id(reaction_name).id = updated_r_name

    return model

def prune_zero_abundance(model, zero_abundance_species):
    print('Pruning metabolites and Reactions from Zero-Abundance Species in Sample')
    zero_prefixes = {f"{species}_" for species in zero_abundance_species}
    metabolites_to_remove = [
        met for met in model.metabolites 
        if any(met.id.startswith(prefix) for prefix in zero_prefixes)
    ]
    # Setting destructive to True also removes all associated reactions,
    # so no need to loop over reactions to remove (will prune 0)
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

    print_out("Building global community model", padding_str='*')
    all_species = abundance_df.index.tolist()
    first_path = os.path.join(mod_dir, all_species[0] + ".mat")
    print_out("Added model: " + all_species[0], padding_str="*")
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
        print_out("Added model: " + species.split("/")[-1], padding_str="*")

    print_out("Finished adding GEM reconstructions to community")

    print_out("Adding diet and fecal compartments")
    clean_model = clean_community(model=global_model)
    print_out("Done adding diet and fecal compartments")

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

    model = prune_zero_abundance(model, zero_abundance_species=zero_species)

    # Add a community biomass reaction to the model
    print_out("Adding community biomass reaction")
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

def build_global_coupling_constraints(model, species_list, coupling_factor=400):
    """
    Build coupling constraints equivalent to MATLAB's coupleRxnList2Rxn
    Creates constraints: v_i - coupling_factor * v_biomass <= 0
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
        biomass_idx = rxn_id_to_index[biomass_rxn.id]

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
    
    # Find which constraint rows to keep
    # keep_rows = []
    # for i, ctr_name in enumerate(tqdm(global_ctrs)):
    #     # Extract species name from constraint name (e.g., "slack_Bacteroides_sp_2_1_33B_IEX_12ppd_S[u]tr")
    #     species_name = ctr_name.split("slack_")[1]
    #     matched_name = next((name for name in present_species if species_name.startswith(name)), None)
    #     if matched_name:
    #         keep_rows.append(i)

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



def compy(abun_filepath, mod_filepath, out_filepath, diet_filepath=None, workers=1):
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
            E.g. "~/compy_results/"
        dietpath: path to the AGORA compatible diet (for community model) .csv file
            E.g. "~/data_input/AverageEU_diet_fluxes.csv"
  
    OUTPUTS:
        All sample community models to a specified local folder
    """

    print_out("Starting ComPy pipeline")
    print_out("Reading abundance file")

    sample_info = pd.read_csv(abun_filepath)
    sample_info.rename(columns={list(sample_info)[0]:"species"}, inplace=True)
    sample_info.set_index("species", inplace=True)

    # # Setting solver and model configurations
    # solver = [s for s in ["cplex", "gurobi", "osqp", "glpk"] \
    #           if s in cobra.util.solver.solvers][0]

    global_model, global_C, global_d, global_dsense, global_ctrs = build_global_model(sample_info, mod_filepath)
    samples = sample_info.columns.tolist()

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(build_sample_model, s, global_model, sample_info, abun_filepath, out_filepath, 
                                   diet_filepath, global_C, global_d, global_dsense, global_ctrs)
                   for s in samples]
        for f in tqdm(futures, desc='Building sample models'):
            f.result()