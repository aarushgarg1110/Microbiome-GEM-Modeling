import numpy as np
from scipy.sparse import csr_matrix

def make_dummy_model(num_mets, num_rxns):
    """
    Python equivalent of MATLAB's makeDummyModel.m
    Creates an empty model structure with all necessary fields
    """
    model_dict = {
        'mets': np.array([''] * num_mets, dtype=object),
        'rxns': np.array([''] * num_rxns, dtype=object),
        'c': np.zeros(num_rxns),
        'S': csr_matrix((num_mets, num_rxns)),
        'lb': np.zeros(num_rxns),
        'ub': np.zeros(num_rxns),
        'b': np.zeros(num_mets),
        'metNames': np.array([''] * num_mets, dtype=object),
        'metFormulas': np.array([''] * num_mets, dtype=object),
        'rxnNames': np.array([''] * num_rxns, dtype=object),
        'subSystems': np.array([''] * num_rxns, dtype=object),
        'genes': np.array([], dtype=object),
        'rxnGeneMat': csr_matrix((num_rxns, 0)),
        'rules': np.array([''] * num_rxns, dtype=object),
        'grRules': np.array([''] * num_rxns, dtype=object),
        'comments': np.array([''] * num_rxns, dtype=object),
        'citations': np.array([''] * num_rxns, dtype=object),
        'confidenceScores': np.array([''] * num_rxns, dtype=object),
        'ecNumbers': np.array([''] * num_rxns, dtype=object),
        'rxnKeggID': np.array([''] * num_rxns, dtype=object),
        'metChEBIID': np.array([''] * num_mets, dtype=object),
        'metCharges': np.zeros(num_mets),
        'metHMDB': np.array([''] * num_mets, dtype=object),
        'metInchiString': np.array([''] * num_mets, dtype=object),
        'metKeggID': np.array([''] * num_mets, dtype=object),
        'metSmile': np.array([''] * num_mets, dtype=object),
        'metPubChemID': np.array([''] * num_mets, dtype=object),
        'csense': np.array(['E'] * num_mets, dtype='<U1'),
        'osenseStr': 'max',
        # Initialize coupling constraint fields
        'C': csr_matrix((0, num_rxns)),
        'd': np.array([]).reshape(-1, 1),
        'dsense': np.array([], dtype='<U1'),
        'ctrs': np.array([], dtype=object),
        'name': ''
    }
    return model_dict