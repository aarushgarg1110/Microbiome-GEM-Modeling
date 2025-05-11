import os
import csv

import pandas as pd
import cobra.io
from cobra import Reaction, Metabolite


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
            gen_mets.append((model.reactions.get_by_id(reac.id).reaction).split(" <=> ")[0])
    gen_mets = set(gen_mets)

    # Creating the diet compartment
    for met_name in gen_mets:
        # Get only the metabolite
        clean_met_name = met_name.split("[")[0]

        # Set to format metabolite[d]
        d_met_name = clean_met_name + "[d]"

        # Reaction name DUt_metabolite
        dut_name = "DUt_" + clean_met_name

        # Creating the [d] exchange reactions
        # E.g., EX_2omxyl[d]: 2omxyl[d] <=>
        # Ensure there are no duplicates if user runs function more than once
        if d_met_name not in [met.id for met in model.metabolites]:
            reac_name = "EX_" + d_met_name
            reaction = Reaction(reac_name)
            reaction.name = d_met_name + "diet exchange"
            reaction.subsystem = " "
            reaction.lower_bound = -1000.
            reaction.upper_bound = 1000.
            model.add_reactions([reaction])
            model.add_metabolites([Metabolite(d_met_name, formula=" ", name="", compartment="d")])
            new_dietreact = model.reactions.get_by_id(reac_name)
            new_dietreact.add_metabolites({model.metabolites.get_by_id(d_met_name): -1})

        # Creating the [d] transport reactions to lumen
        # E.g., DUt_4hbz: 4hbz[d] --> 4hbz[u]
        # Ensure there are no duplicates if user runs function more than once
        if dut_name not in [reac.id for reac in model.reactions]:
            dut_formula = f"{d_met_name} --> {met_name}"
            reac_name = dut_name
            reaction = Reaction(reac_name)
            reaction.name = dut_name + "diet to lumen"
            reaction.subsystem = " "
            reaction.lower_bound = 0.
            reaction.upper_bound = 1000.
            model.add_reactions([reaction])

            # Adding the correct d --> u formula to the reaction
            reaction.reaction = dut_formula

    # Creating the fecal ([fe]) compartment
    for met_name in gen_mets:
        # Get only the metabolite
        clean_met_name = met_name.split("[")[0]

        # metabolite[d]
        fe_met_name = clean_met_name + "[fe]"

        # reaction name UFEt_metabolite
        ufet_name = "UFEt_" + clean_met_name

        # Creating the fe exchange reactions
        # E.g., EX_4abut[fe]: 4abut[fe] <=>
        # Ensure there are no duplicates if user runs function more than once
        if fe_met_name not in [met.id for met in model.metabolites]:
            reac_name = "EX_" + fe_met_name
            reaction = Reaction(reac_name)
            reaction.name = fe_met_name + "fecal exchange"
            reaction.subsystem = " "
            reaction.lower_bound = -1000.
            reaction.upper_bound = 1000.
            model.add_reactions([reaction])
            model.add_metabolites([Metabolite(fe_met_name, formula=" ", name="", compartment="fe")])
            new_fe_react = model.reactions.get_by_id(reac_name)
            new_fe_react.add_metabolites({model.metabolites.get_by_id(fe_met_name): -1})

        # Creating the [fe] transport reactions to lumen
        # E.g., UFEt_arabinoxyl: arabinoxyl[u] --> arabinoxyl[fe]
        # Ensure there are no duplicates if user runs function more than once
        if ufet_name not in [reac.id for reac in model.reactions]:
            ufet_formula = f"{met_name} --> {fe_met_name}"
            reac_name = ufet_name
            reaction = Reaction(reac_name)
            reaction.name = ufet_name + "diet to lumen"
            reaction.subsystem = ' '
            reaction.lower_bound = 0.
            reaction.upper_bound = 1000.
            model.add_reactions([reaction])

            # Adding the corrected d --> u formula to the reaction
            reaction.reaction = ufet_formula

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
    # and remove species with abundances < 0.0000001
    norm_abund = pd.read_csv(abun_path)
    norm_abund = norm_abund.sort_values("X", ascending=True)
    norm_abund = norm_abund.reset_index(drop=True)
    norm_abund = norm_abund.loc[~((norm_abund[sample_com] < 0.0000000))]
    norm_abund_list = norm_abund[sample_com].tolist()

    # Creating the community biomass reaction
    reaction = Reaction("communityBiomass")
    reaction.name = "community biomass "
    reaction.subsystem = " "
    reaction.lower_bound = 0.
    reaction.upper_bound = 1000.
    # reaction.add_metabolites(com_biomass_dic)
    model.add_reactions([reaction])
    community_biomass = model.reactions.communityBiomass
    counter = 0

    # Adding the metabolites with their weighted abundances
    com_biomass_dict = {}
    for biomass in biomass_mets_list:
        com_biomass_dict[biomass] = -float(norm_abund_list[counter])
        counter += 1
    community_biomass.add_metabolites(metabolites_to_add=com_biomass_dict, combine=True)

    # Adding the microbeBiomass metabolite
    model.add_metabolites([Metabolite("microbeBiomass[u]", formula=" ", \
                                      name="product of community biomass", compartment="u"),])
    community_biomass.add_metabolites({model.metabolites.get_by_id("microbeBiomass[u]"): 1})

    # Adding the exchange reaction compartment
    reac_name = "EX_microbeBiomass[fe]"
    reaction = Reaction(reac_name)
    reaction.name = reac_name + "fecal exchange"
    reaction.subsystem = " "
    reaction.lower_bound = 0.
    reaction.upper_bound = 1000.
    model.add_reactions([reaction])

    # Adding a microbeBiomass [fe] metabolite
    model.add_metabolites([Metabolite("microbeBiomass[fe]", formula=" ", \
                                      name="product of community biomass", compartment="fe"),])

    # Adding the exchange reaction
    new_fe_react = model.reactions.get_by_id("EX_microbeBiomass[fe]")
    new_fe_react.add_metabolites({model.metabolites.get_by_id("microbeBiomass[fe]"): -1})

    # Adding the UFEt reaction
    ufet_formula = "microbeBiomass[u] --> microbeBiomass[fe]"
    reac_name = "UFEt_microbeBiomass"
    reaction = Reaction(reac_name)
    reaction.name = "UFEt_microbeBiomassdiet to lumen"
    reaction.subsystem = " "
    reaction.lower_bound = 0.
    reaction.upper_bound = 1000.
    model.add_reactions([reaction])

    # Adding the correct d --> u formula to the reaction
    reaction.reaction = ufet_formula

    return model


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
                            metabolite_name = met.id
                            model.metabolites.get_by_id(metabolite_name).compartment = "c"
                            no_c_name = metabolite_name.replace("[c]", "")
                            updated_m_name = short_species_name + "_" + no_c_name + "[c]"
                            model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                            model.metabolites.get_by_id(updated_m_name).compartment = "c"
                        if "[p]" in str(met.id):
                            metabolite_name = met.id
                            model.metabolites.get_by_id(metabolite_name).compartment = ""
                            no_c_name = metabolite_name.replace("[p]", "")
                            updated_m_name = short_species_name + "_" + no_c_name + "[p]"
                            model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                            model.metabolites.get_by_id(updated_m_name).compartment = "p"

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
                        metabolite_name = met.id
                        model.metabolites.get_by_id(metabolite_name).compartment = "c"
                        no_c_name = metabolite_name.replace("[c]", "")
                        updated_m_name = short_species_name + "_" + no_c_name + "[c]"
                        model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                        model.metabolites.get_by_id(updated_m_name).compartment = "c"
                if "[e]" in met.id:
                    if short_species_name not in str(met.id):
                        metabolite_name = met.id
                        model.metabolites.get_by_id(metabolite_name).compartment = "u"
                        no_c_name = metabolite_name.replace("[e]", "")
                        updated_m_name = short_species_name + "_" + no_c_name + "[u]"
                        model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                        model.metabolites.get_by_id(updated_m_name).compartment = "u"
                if "[p]" in met.id:
                    if short_species_name not in str(met.id):
                        metabolite_name = met.id
                        model.metabolites.get_by_id(metabolite_name).compartment = "p"
                        no_c_name = metabolite_name.replace("[p]", "")
                        updated_m_name = short_species_name + "_" + no_c_name + "[p]"
                        model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                        model.metabolites.get_by_id(updated_m_name).compartment = "p"

    # Adding the tagged intracellular reactions for the extracellular metabolites to the model
    # For all extracellular species specific metabolites, add an exchange reaction
    for met in model.metabolites:
        if "[u]" in met.id:
            if short_species_name in met.id:
                replacing_specname = short_species_name + "_"
                general_name = met.id.replace(replacing_specname, "")
                iex_formula = f"{general_name} <=> {met.id}"
                iex_reaction_name = short_species_name + "_IEX_" + general_name + "tr"
                reaction = Reaction(iex_reaction_name)
                reaction.name = short_species_name + "_IEX"
                reaction.lower_bound = 0.
                reaction.upper_bound = 1000.
                reaction.subsystem = ""
                model.add_reactions([reaction])
                reaction.reaction = iex_formula

    # Ensure each reaction has the species tag on it
    for rxn in model.reactions:
        if short_species_name not in rxn.id:
            reaction_name = rxn.id
            updated_r_name =  short_species_name + "_" + reaction_name
            model.reactions.get_by_id(reaction_name).id = updated_r_name

    # Ensure each [c] compartment metabolite has the species tag on it
    for met in model.metabolites:
        if "[c]" in met.id:
            if short_species_name not in met.id:
                metabolite_name = met.id
                model.metabolites.get_by_id(metabolite_name).compartment = "c"
                no_c_name = metabolite_name.replace("[c]", "")
                updated_m_name = short_species_name + "_" + no_c_name + "[c]"
                model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                model.metabolites.get_by_id(updated_m_name).compartment = "c"

    # Ensure each [p] compartment metabolite has the species tag on it
    for met in model.metabolites:
        if "[p]" in met.id:
            if short_species_name not in met.id:
                metabolite_name = met.id
                model.metabolites.get_by_id(metabolite_name).compartment = "p"
                no_p_name = metabolite_name.replace("[p]", "")
                updated_m_name = short_species_name + "_" + no_p_name + "[p]"
                model.metabolites.get_by_id(metabolite_name).id = updated_m_name
                model.metabolites.get_by_id(updated_m_name).compartment = "p"

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
    for met in model.metabolites:
        for species in zero_abundance_species:
            if met.id.startswith(f"{species}_"):
                metabolites_to_remove.append(met)
                print('Pruning Zero-Abundance metabolite:', met.id)
                break
    model.remove_metabolites(metabolites_to_remove)

    reactions_to_remove = []
    for rxn in model.reactions:
        for species in zero_abundance_species:
            if rxn.id.startswith(f"{species}_"):
                reactions_to_remove.append(rxn.id)
                print('Pruning Zero-Abundance reaction:', rxn.id)
                break
    model.remove_reactions(reactions_to_remove)
    return model

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

    # Loading the abundance input file and formatting it correctly
    sample_info = pd.read_csv(abun_filepath)
    sample_info.rename(columns={list(sample_info)[0]:"species"}, inplace=True)
    sample_info.set_index("species", inplace=True)
    # sample_list = sampleinfo.columns.tolist()

    # Setting solver and model configurations
    solver = [s for s in ["cplex", "gurobi", "osqp", "glpk"] \
              if s in cobra.util.solver.solvers][0]

    # Make a counter for the number of models
    model_count = 0
    model_count_max = len(sample_info.columns)

    # Joining the models
    # Iterate through each column in sample info to look at each sample
    for sample in sample_info.columns:
        print_out("Working on sample: " + sample, padding_str="*")

    # For one specific community sample, create a dictionary of the species and their abundances
        model_path_abun = {}

        # Add all species and abundances to model_path_abun_dic in format of "species name : abundace"
        for abundance in range(len(sample_info[sample])):
            species_name = sample_info.index[abundance]
            model_path_abun[mod_filepath + species_name + '.mat'] = sample_info[sample][abundance]


        # Create the first model in the sample, and add it to the final model
        first_species = list(model_path_abun.keys())[0]
        print_out("Added model: " + first_species, padding_str="*")
        first_model = cobra.io.load_matlab_model(first_species)
        final_model = species_to_community(model=first_model, species_model_name=first_species)

        # Loop through the remaining species model in this sample and merge them into
        # the (combined) final_model
        for species_model in list(model_path_abun.keys())[1:]:
            model = cobra.io.load_matlab_model(species_model)
            adding_model = species_to_community(model=model, species_model_name=species_model)
            rxn_ids = set(rxn.id for rxn in final_model.reactions)
            new = [rxn.id for rxn in adding_model.reactions if rxn.id not in rxn_ids]
            final_model.add_reactions(adding_model.reactions.get_by_any(new))
            print_out("Added model: " + species_model.split("/")[-1], padding_str="*")

        print_out("Finished adding GEM reconstructions to community")

        # Add the diet and fecal compartments to model
        print_out("Adding diet and fecal compartments")
        clean_model = clean_community(model=final_model)

        # maybe prune here before com_biomass?
        # Prune zero-abundance species from the model
        zero_abundance_species = [species.split("/")[-1].split(".")[0] for species, abundance in model_path_abun.items() if abundance < 0.0000001]
        clean_model = prune_zero_abundance(clean_model, zero_abundance_species)

        # Add a community biomass reaction to the model
        print_out("Adding community biomass reaction")
        clean_model = com_biomass(model=clean_model, abun_path=abun_filepath, sample_com=sample)

        # Ensuring the reversablity fits all compartments
        # for reac in clean_model.reactions:
            # clean_model.reactions.get_by_id(reac.id).lower_bound = -1000
        print_out("Setting reversible reaction bounds to -1000")
        for reac in [i for i in clean_model.reactions if "DUt" in i.id or "UFEt" in i.id]:
            clean_model.reactions.get_by_id(reac.id).lower_bound = 0.

        # Setting the given diet to the model
        print_out("Setting diet reaction bounds")
        if diet_filepath is not None:
            diet_dic = {}
            with open(diet_filepath) as f:
                next(f)
                for key, *values in csv.reader(f, delimiter="\t"):
                    diet_dic[key.strip()] = values[0]

            for reac in [i.id for i in clean_model.reactions if "[d]" in i.id]:
                if reac in diet_dic.keys():
                    clean_model.reactions.get_by_id(reac).lower_bound = -float(diet_dic[reac])
                    print_out("Metabolite added to diet: " + reac)
                else:
                    clean_model.reactions.get_by_id(reac).lower_bound = 0.

        # Setting the new community biomass as the objective
        clean_model.objective = "communityBiomass"
        if not os.path.exists(out_filepath):
            os.makedirs(out_filepath)

        # Set the path of the community model to be saved
        microbiome_model =  out_filepath + sample + "_communitymodel_final.json"

        # Saving the model
        cobra.io.save_json_model(clean_model, microbiome_model)

        model_count += 1
        print_out("Saved model " + str(model_count) + " of " + str(model_count_max) + ": " + \
                  microbiome_model.split("/")[-1], padding_str="*")

    print_out("Finished ComPy pipeline")
