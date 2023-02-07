import numpy as np
from flask import request

from chemquest_website import app, db, fpe

@app.route('/generate_mcq')
def generate_mcq():
    '''Returns a set of 3 random similar molecules to the input SMILES. For use in MCQ questions.'''
    input_smiles = request.args.get('input_smiles')
    results = fpe.on_disk_similarity(input_smiles, threshold = 0.5, n_workers = 4)[1:] # remove the first element as it's the input molecule itself, with similarity 1
    choices = np.random.choice(results, 3, replace=False) # format: [[mol_id, similarity], ...]
    choices_smiles = list(db.molecules_original.aggregate([{ 
        "$facet": {
            "query0": [
                { "$match": {"mol_id": int(choices[0][0])} },
                { "$limit": 1 }
            ],
            "query1": [
                { "$match": {"mol_id": int(choices[1][0])} },
                { "$limit": 1 }
            ],
            "query2": [
                { "$match": {"mol_id": int(choices[2][0])} },
                { "$limit": 1 }
            ],
        }
    }]))[0] # format: {'query0': [{'_id': ..., 'SMILES': ..., 'mol_id': ...}], ...} the [0] is because it's returned as an iterable with 1 element
    choices_processed = [
        {"mol_id": int(choices[0][0]), "SMILES": choices_smiles["query0"][0]["SMILES"], "similarity": float(choices[0][1])},
        {"mol_id": int(choices[1][0]), "SMILES": choices_smiles["query1"][0]["SMILES"], "similarity": float(choices[1][1])},
        {"mol_id": int(choices[2][0]), "SMILES": choices_smiles["query2"][0]["SMILES"], "similarity": float(choices[2][1])},
    ]
    return choices_processed