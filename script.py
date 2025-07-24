import sys
sys.path.insert(0, './src')
from src.diet_adaptation import *
from src.migemox import migemox

migemox(abun_filepath="test_data_input/normCoverageReduced.csv",
        mod_filepath='test_data_input/AGORA103',
        out_filepath="Python_Models/Personalized",
        diet_filepath='test_data_input/AverageEU_diet_fluxes.txt',
        workers=2)

clean_samp_names, organisms, ex_mets = get_individual_size_name(
    abun_file_path='test_data_input/normCoverageReduced.csv',
    mod_path='test_data_input/AGORA103'
)

exchanges, net_production, net_uptake = simulate_microbiota_models(
    sample_names=clean_samp_names,
    ex_mets=ex_mets,
    model_dir='Python_Models/Personalized',
    diet_file='test_data_input/AverageEU_diet_fluxes.txt',
    res_path='Python_Models',
    biomass_bounds = (0.4, 1.0),
    solver='cplex',
    workers=2
)

net_secretion_df, net_uptake_df = collect_flux_profiles(
    samp_names=clean_samp_names,
    exchanges=sorted(exchanges),
    net_production=net_production,
    net_uptake=net_uptake
)

net_secretion_df.to_csv('python_net_secretion_fluxes.csv', index=True, index_label='Net secretion')
net_uptake_df.to_csv('python_net_uptake_fluxes.csv', index=True, index_label='Net uptake')

### downstream now
'''
Targeted analysis: Strain-level contributions to metabolites of interest
For metabolites of particular interest (e.g., for which the community-wide 
secretion potential was significantly different between disease cases and controls), 
the strains consuming and secreting the metabolites in each sample may be computed. 
This will yield insight into the contribution of each strain to each metabolite. 
Note that for metabolites for which the community wide secretion potential did 
not differ, the strains contributing to metabolites may still be significantly 
different.

The first step for the preparation of targeted analyses is the export of models 
that had already been constrained with the simulated dietary regime. They can 
be found in a subfolder called "Diet" in the results folder. Now, will set the 
input variable modPath to the folder with personalized models constrained with 
the simulated diet.

We will define a list of metabolites to analyze (default: all exchanged metabolites 
in the models). As an example, we will take acetate and formate. Enter the VMH 
IDs of target metabolites as a cell array.

The output 'minFluxes' shows the fluxes in the reverse direction through all 
internal exchange reactions that had nonzero flux for each analyzed metabolite. 
The output 'maxFluxes' shows the corresponding forward fluxes. 'fluxSpans' shows 
the distance between minimal and maximal fluxes for each internal exchange reaction 
with nonzero flux for each metabolite.
'''
mets = [x.split('EX_')[1].split('[')[0] for x in ex_mets]
from src.downstream_analysis.predict_microbe_contribution import *
minf, maxf, fluxs = predict_microbe_contributions('Matlab_Models\Diet', mets_list=mets, workers=2)