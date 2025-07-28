"""
I/O Utilities for MiGEMox Pipeline

This module provides functions for loading input data (e.g., abundance files),
and preparing model data for saving to disk. It centralizes file operations
to improve code organization and reusability.
"""

import pandas as pd
import os
import re
from cobra.io import load_matlab_model
from cobra.util import create_stoichiometric_matrix
from scipy.io import savemat
import numpy as np
from scipy.sparse import csr_matrix

def get_individual_size_name(abun_file_path: str, mod_path: str) -> tuple:
    """
    Reads abundance data from a CSV file and extracts sample names, organisms,
    and extracellular metabolites from associated model files or a dict of already
    loaded models.

    Args:
        abun_file_path: Path to the abundance CSV file.
        mod_path: Path to the directory containing organism model files (.mat).

    Returns:
        A tuple containing:
            - clean_samp_names: List of cleaned sample names (valid Python identifiers).
            - organisms: List of organism names from the abundance file.
            - ex_mets: List of sorted unique extracellular metabolites found in the models.

    Raises:
        ValueError: If there's an error reading the abundance file or loading a model.
        FileNotFoundError: If a model file is not found.
    """
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
            valid_name = re.sub('[^0-9a-zA-Z_]', '_', valid_name)
            if not valid_name[0].isalpha():
                valid_name = 'sample_' + valid_name
        clean_samp_names.append(valid_name)

    # Step 3: Load models and extract [e] metabolites
    ex_mets = set()
    for organism in organisms:
        model_file = os.path.join(mod_path, organism + '.mat')
        if not os.path.exists(model_file):
            raise FileNotFoundError(f"Model file not found: {model_file}")
        try:
            model = load_matlab_model(model_file)
        except Exception as e:
            raise ValueError(f"Error loading model {organism}: {e}")

        # Extract extracellular metabolites (assuming suffix '[e]')
        mets = [met.id for met in model.metabolites if met.id.endswith('[e]')]
        ex_mets.update(mets)

    return clean_samp_names, organisms, list(sorted(ex_mets))

def make_mg_pipe_model_dict(model, C=None, d=None, dsense=None, ctrs=None):
    """
    Creates a dictionary representation of a cobra model, including additional
    parameters (C, d, dsense, ctrs) for saving in a .mat file format compatible
    with mgPipe/Microbiome Modeling Toolbox.

    Args:
        model (cobra.Model): The cobra model object.
        C (scipy.sparse.csr_matrix, optional): Coupling matrix. Defaults to None.
        d (numpy.ndarray, optional): Right-hand side of coupling constraints. Defaults to None.
        dsense (numpy.ndarray, optional): Sense of coupling constraints ('E', 'L', 'G'). Defaults to None.
        ctrs (list, optional): Names of coupling constraints. Defaults to None.

    Returns:
        dict: A dictionary containing model components and additional parameters
              suitable for `scipy.io.savemat`.
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
    S = create_stoichiometric_matrix(model)

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