from rdkit import Chem
from rdkit.Chem import rdmolops, rdChemReactions

class ChemicalTransformer:
    def __init__(self, transform_function: callable = None):
        if transform_function is None:
            self.transform = self.transform_placeholder
        else:
            self.transform = transform_function

    def transform_placeholder(self, mol: Chem.Mol):
        raise NotImplementedError

    @classmethod
    def from_substructs(cls, in_substruct: Chem.Mol | str, out_substruct: Chem.Mol | str):
        if isinstance(in_substruct, str):
            in_substruct = Chem.MolFromSmarts(in_substruct)

        if isinstance(out_substruct, str):
            out_substruct = Chem.MolFromSmarts(out_substruct)
        

        def transform_function(mol: Chem.Mol):
            results = rdmolops.ReplaceSubstructs(mol, in_substruct, out_substruct)
            # catch case where there's no matching substruct (function will return the original molecule)
            if len(results) == 1 and Chem.MolToSmiles(results[0]) == Chem.MolToSmiles(mol):
                return []
            return results
        
        return cls(transform_function)
    
    @classmethod
    def from_reaction(cls, reaction: rdChemReactions.ChemicalReaction | str):
        if isinstance(reaction, str):
            reaction = rdChemReactions.ReactionFromSmarts(reaction)

        def transform_function(mol: Chem.Mol):
            rxn_results = reaction.RunReactants((mol,))
            return [result[0] for result in rxn_results]
        
        return cls(transform_function)