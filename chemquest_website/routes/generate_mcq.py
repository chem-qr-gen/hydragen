from flask import request
from rdkit import Chem

from chemquest_website import app
from chemquest_website.helpers.mcq_generator import generate_mcq_options


@app.route('/generate_mcq')
def generate_mcq():
    '''Returns a set of 3 random similar molecules to the input SMILES. For use in MCQ questions.'''
    input_smiles = request.args.get('input_smiles')
    
    return generate_mcq_options(Chem.MolFromSmiles(input_smiles))