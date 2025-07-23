import sys
sys.path.insert(0, './src')
from src.diet_adaptation import *
from src.compy import compy

compy(abun_filepath="test_data_input/normCoverageReduced.csv",
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
mets = [x.split('EX_')[1].split('[')[0] for x in ex_mets]
from src.downstream_analysis.predict_microbe_contribution import *
minf, maxf, fluxs = predict_microbe_contributions('Matlab_Models\Diet', mets_list=mets, workers=2)