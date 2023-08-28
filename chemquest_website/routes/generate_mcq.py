import random

from flask import request
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sqlalchemy import func

from chemquest_website import app, engine, meta

def get_options(source_smiles: str, num_options: int = 3) -> list[tuple[str, float]]:
    '''
    Returns a random list of molecules that are similar to the input SMILES string.
    Priority is given to molecules with the same molecular weight as the input molecule.
    '''
    # Get the molecular weight of the input molecule
    source_mol = Chem.MolFromSmiles(source_smiles)
    source_mol_weight = Descriptors.MolWt(source_mol)
    
    # Get a list of molecules with the same molecular weight as the input molecule
    molecules_table = meta.tables["molecules_original"]
    with engine.connect() as conn:
        rows = conn.execute(molecules_table.select().where(func.abs(molecules_table.c.molwt - source_mol_weight) < 0.5)).fetchall()

        # If there aren't enough molecules, get a wider range
        if len(rows) < num_options:
            rows = conn.execute(molecules_table.select().where(func.abs(molecules_table.c.molwt - source_mol_weight) < 10)).fetchall()

    smiles_to_compare = [row._mapping["smiles"] for row in rows]
    mols_to_compare = [Chem.MolFromSmiles(smiles) for smiles in smiles_to_compare]

    # Generate fingerprints
    fpgen = rdFingerprintGenerator.GetRDKitFPGenerator()
    source_fp = fpgen.GetFingerprint(source_mol)
    fps_to_compare = [fpgen.GetFingerprint(mol) for mol in mols_to_compare]

    # Compare similarity
    similarity = DataStructs.BulkTanimotoSimilarity(source_fp, fps_to_compare)

    # Sort by decreasing similarity
    smiles_and_similarities = list(zip(smiles_to_compare, similarity))
    smiles_and_similarities.sort(key=lambda x: x[1], reverse=True)

    if len(smiles_and_similarities) < num_options * 3:
        top_smiles = smiles_and_similarities[1:]
    else:
        top_smiles = smiles_and_similarities[1:num_options * 3]
    
    # Pick random options from the top few
    top_options = random.sample(top_smiles, num_options)
    return top_options




@app.route('/generate_mcq')
def generate_mcq():
    '''Returns a set of 3 random similar molecules to the input SMILES. For use in MCQ questions.'''
    input_smiles = request.args.get('input_smiles')
    
    # Get choices for the question
    choices = get_options(input_smiles)
    choices_processed = [{"SMILES": smiles, "similarity": similarity} for smiles, similarity in choices]

    # TODO: we need a better "similar molecule" algorithm. 
    # This algo tends to generate molecules that are much larger than the original, so it's easy to tell what the answer is.
    # results = fpe.on_disk_similarity(input_smiles, threshold = 0.5, n_workers = 4)[1:] # remove the first element as it's the input molecule itself, with similarity 1
    # if len(results) < 3:
    #     results = fpe.on_disk_similarity(input_smiles, threshold = 0.45, n_workers = 4)[1:] # quick fix for small molecules
    
    # choices = np.random.choice(results, 3, replace=False) # format: [[mol_id, similarity], ...]
    # choices_smiles = list(db.molecules_original.aggregate([{ 
    #     "$facet": {
    #         "query0": [
    #             { "$match": {"mol_id": int(choices[0][0])} },
    #             { "$limit": 1 }
    #         ],
    #         "query1": [
    #             { "$match": {"mol_id": int(choices[1][0])} },
    #             { "$limit": 1 }
    #         ],
    #         "query2": [
    #             { "$match": {"mol_id": int(choices[2][0])} },
    #             { "$limit": 1 }
    #         ],
    #     }
    # }]))[0] # format: {'query0': [{'_id': ..., 'SMILES': ..., 'mol_id': ...}], ...} the [0] is because it's returned as an iterable with 1 element
    # choices_processed = [
    #     {"mol_id": int(choices[0][0]), "SMILES": choices_smiles["query0"][0]["SMILES"], "similarity": float(choices[0][1])},
    #     {"mol_id": int(choices[1][0]), "SMILES": choices_smiles["query1"][0]["SMILES"], "similarity": float(choices[1][1])},
    #     {"mol_id": int(choices[2][0]), "SMILES": choices_smiles["query2"][0]["SMILES"], "similarity": float(choices[2][1])},
    # ]
    return choices_processed