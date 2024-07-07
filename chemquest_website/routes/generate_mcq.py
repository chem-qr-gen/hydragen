import random

from flask import request
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sqlalchemy import func

from chemquest_website import app, engine, meta
from chemquest_website.helpers.mcq_generator import generate_mcq_options

# def get_options(source_smiles: str, num_options: int = 3) -> list[tuple[str, float]]:
#     '''
#     Returns a random list of molecules that are similar to the input SMILES string.
#     Priority is given to molecules with the same molecular weight as the input molecule.
#     '''
#     # Get the molecular weight of the input molecule
#     source_mol = Chem.MolFromSmiles(source_smiles)
#     source_mol_weight = Descriptors.MolWt(source_mol)
    
#     # Get a list of molecules with the same molecular weight as the input molecule
#     molecules_table = meta.tables["molecules_original"]
#     with engine.connect() as conn:
#         rows = conn.execute(molecules_table.select().where(func.abs(molecules_table.c.molwt - source_mol_weight) < 0.5)).fetchall()

#         # If there aren't enough molecules, get a wider range
#         if len(rows) < num_options:
#             rows = conn.execute(molecules_table.select().where(func.abs(molecules_table.c.molwt - source_mol_weight) < 10)).fetchall()

#     smiles_to_compare = [row._mapping["smiles"] for row in rows]
#     mols_to_compare = [Chem.MolFromSmiles(smiles) for smiles in smiles_to_compare]

#     # Generate fingerprints
#     fpgen = rdFingerprintGenerator.GetRDKitFPGenerator()
#     source_fp = fpgen.GetFingerprint(source_mol)
#     fps_to_compare = [fpgen.GetFingerprint(mol) for mol in mols_to_compare]

#     # Compare similarity
#     similarity = DataStructs.BulkTanimotoSimilarity(source_fp, fps_to_compare)

#     # Sort by decreasing similarity
#     smiles_and_similarities = list(zip(smiles_to_compare, similarity))
#     smiles_and_similarities.sort(key=lambda x: x[1], reverse=True)

#     if len(smiles_and_similarities) < num_options * 3:
#         top_smiles = smiles_and_similarities[1:]
#     else:
#         top_smiles = smiles_and_similarities[1:num_options * 3]
    
#     # Pick random options from the top few
#     top_options = random.sample(top_smiles, num_options)
#     return top_options




@app.route('/generate_mcq')
def generate_mcq():
    '''Returns a set of 3 random similar molecules to the input SMILES. For use in MCQ questions.'''
    input_smiles = request.args.get('input_smiles')
    
    # Get choices for the question
    choices = generate_mcq_options(Chem.MolFromSmiles(input_smiles))
    choices_processed = [{"SMILES": smiles} for smiles in choices]
    
    return choices_processed