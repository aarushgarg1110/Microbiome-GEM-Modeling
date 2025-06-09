import os
import csv
from tqdm import tqdm
import pandas as pd
import cobra.io
from cobra import Reaction, Metabolite
from concurrent.futures import ProcessPoolExecutor


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
    while len([reac for reac in model.reactions if "_EX_" in reac.id or "(e)" in reac.id]) > 0:
        for reac in model.reactions:
            if "_EX_" in reac.id or "(e)" in reac.id:
                model.reactions.get_by_id(reac.id).remove_from_model()

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
        if d_met_name not in [met.id for met in model.metabolites]:
            reac_name = "EX_" + d_met_name
            reaction = create_rxn(reac_name, d_met_name + 'diet exchange', ' ', -1000., 10000.)
            model.add_reactions([reaction])
            model.add_metabolites([Metabolite(d_met_name, compartment="d")])
            new_dietreact = model.reactions.get_by_id(reac_name)
            new_dietreact.add_metabolites({model.metabolites.get_by_id(d_met_name): -1})

        # Creating the [d] transport reactions to lumen
        # E.g., DUt_4hbz: 4hbz[d] --> 4hbz[u]
        # Ensure there are no duplicates if user runs function more than once
        if dut_name not in [reac.id for reac in model.reactions]:
            dut_formula = f"{d_met_name} --> {met_name}"
            reaction = create_rxn(dut_name, dut_name + 'diet to lumen', ' ', 0., 10000.)
            model.add_reactions([reaction])

            # Adding the correct d --> u formula to the reaction
            reaction.reaction = dut_formula
            reaction.bounds = (0., 10000.)

        # Creating the fe exchange reactions
        # E.g., EX_4abut[fe]: 4abut[fe] <=>
        # Ensure there are no duplicates if user runs function more than once
        if fe_met_name not in [met.id for met in model.metabolites]:
            reac_name = "EX_" + fe_met_name
            reaction = create_rxn(reac_name, fe_met_name + "fecal exchange", ' ', -1000., 10000.)
            model.add_reactions([reaction])
            model.add_metabolites([Metabolite(fe_met_name, compartment="fe")])
            new_fe_react = model.reactions.get_by_id(reac_name)
            new_fe_react.add_metabolites({model.metabolites.get_by_id(fe_met_name): -1})

        # Creating the [fe] transport reactions to lumen
        # E.g., UFEt_arabinoxyl: arabinoxyl[u] --> arabinoxyl[fe]
        # Ensure there are no duplicates if user runs function more than once
        if ufet_name not in [reac.id for reac in model.reactions]:
            ufet_formula = f"{met_name} --> {fe_met_name}"
            reaction = create_rxn(ufet_name, ufet_name + "diet to lumen", ' ', 0., 10000.)
            model.add_reactions([reaction])

            # Adding the corrected d --> u formula to the reaction
            reaction.reaction = ufet_formula
            reaction.bounds = (0., 10000.)

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
    while len([reac for reac in model.reactions if "Biomass" in reac.id]) > 0:
        for reac in model.reactions:
            if "Biomass" in reac.id:
                model.reactions.get_by_id(reac.id).remove_from_model()

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
    norm_abund = norm_abund[norm_abund[sample_com] > 1e-7]

    # Creating the community biomass reaction
    reaction = create_rxn("communityBiomass", "community biomass", ' ', 0., 10000.)
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
    reaction = create_rxn(reac_name, reac_name + "fecal exchange", ' ', -10000., 10000.)
    model.add_reactions([reaction])

    # Adding a microbeBiomass [fe] metabolite
    model.add_metabolites([Metabolite("microbeBiomass[fe]", formula=" ", \
                                      name="product of community biomass", compartment="fe"),])

    # Adding the exchange reaction
    new_fe_react = model.reactions.get_by_id("EX_microbeBiomass[fe]")
    new_fe_react.add_metabolites({model.metabolites.get_by_id("microbeBiomass[fe]"): -1})

    # Adding the UFEt reaction
    ufet_formula = "microbeBiomass[u] --> microbeBiomass[fe]"
    reaction = create_rxn("UFEt_microbeBiomass", "UFEt_microbeBiomassdiet to lumen", ' ', 0., 10000.)
    model.add_reactions([reaction])

    # Adding the correct d --> u formula to the reaction
    reaction.reaction = ufet_formula
    reaction.bounds = (0., 10000.)

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
    short_species_name = species_model_name.split("/")[-1].split(".")[0]

    # Removing all exchange reactions except for the biomass reaction
    for rxn in model.reactions:
        if "EX_" in rxn.id and "biomass" not in rxn.id:
            model.remove_reactions(model.reactions.get_by_id(rxn.id))

    # Change the intracellualr reaction from [c] --> [c]
    for rxn in model.reactions:
        if "[e]" not in rxn.reaction:
            if "[c]" in rxn.reaction:
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
                # print(model.reactions.get_by_id(updated_r_name).bounds)

    # Workign on the extracellular reactions from [c] --> [e]
    for rxn in model.reactions:
        if "[e]" in rxn.reaction:
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
    metabolites_to_remove = []
    for met in tqdm(model.metabolites, desc="Checking metabolites"):
        for species in zero_abundance_species:
            if met.id.startswith(f"{species}_"):
                metabolites_to_remove.append(met)
                break
    model.remove_metabolites(metabolites_to_remove, destructive=True)
    print(f"Pruned {len(metabolites_to_remove)} Metabolites")
    
    reactions_to_remove = []
    for rxn in tqdm(model.reactions, desc='Checking reactions'):
        for species in zero_abundance_species:
            if rxn.id.startswith(f"{species}_"):
                reactions_to_remove.append(rxn.id)
                break
    model.remove_reactions(reactions_to_remove)
    print(f"Pruned {len(reactions_to_remove)} Reactions")

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

    return clean_model

def build_sample_model(sample_name: str, global_model: cobra.Model, abundance_df: pd.DataFrame, 
                       abun_path: str, out_dir: str, diet_path: str = None):
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
    sample_abun = abundance_df[sample_name]

    # Prune zero-abundance species from the model
    zero_species = [sp for sp in sample_abun.index if sample_abun[sp] < 1e-7]
    model = prune_zero_abundance(model, zero_abundance_species=zero_species)

    # Add a community biomass reaction to the model
    print_out("Adding community biomass reaction")
    model = com_biomass(model=model, abun_path=abun_path, sample_com=sample_name)

    # Ensuring the reversablity fits all compartments
    print_out("Setting reversible reaction bounds to -1000")
    for reac in [r for r in model.reactions if "DUt" in r.id or "UFEt" in r.id]:
        reac.lower_bound = 0.

    # Setting the given diet to the model
    print_out("Setting diet reaction bounds")
    if diet_path:
        diet_bounds = {}
        with open(diet_path) as f:
            next(f)
            for key, *vals in csv.reader(f, delimiter="\t"):
                diet_bounds[key.strip()] = float(vals[0])

        for rxn in model.reactions:
            if "[d]" in rxn.id:
                if rxn.id in diet_bounds:
                    rxn.lower_bound = -diet_bounds[rxn.id]
                    print_out("Metabolite added to diet: " + rxn.id)
                else:
                    rxn.lower_bound = -1000

    # Setting the new community biomass as the objective
    # Setting EX_microbeBiomass[fe] as objective to match MATLAB mgPipe
    model.objective = "EX_microbeBiomass[fe]"

    os.makedirs(out_dir, exist_ok=True)
    save_path = os.path.join(out_dir, f"{sample_name}_communitymodel_final.json")
    cobra.io.save_json_model(model, save_path)

    return save_path

def compy(abun_filepath, mod_filepath, out_filepath, diet_filepath=None):
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

    global_model = build_global_model(sample_info, mod_filepath)
    samples = sample_info.columns.tolist()

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(build_sample_model, s, global_model, sample_info, abun_filepath, out_filepath, diet_filepath)
                   for s in samples]
        for f in tqdm(futures, desc='Building sample models'):
            f.result()