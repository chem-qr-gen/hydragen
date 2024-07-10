from box import Box
from rdkit import Chem

# Functional group SMARTS obtained from https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html

funcgroups_mols = Box({
    "alkene": Chem.MolFromSmarts("[$([CX2](=C)=C)]"),
    "alkyne": Chem.MolFromSmarts("[$([CX2]#C)]"),
    "aldehyde": Chem.MolFromSmarts("[CX3H1](=O)[#6]"),
    "ketone": Chem.MolFromSmarts("[#6][CX3](=O)[#6]"),
    "carboxylic acid": Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]"),
    "ester": Chem.MolFromSmarts("[#6][CX3](=O)[OX2H0][#6]"),
    "amide": Chem.MolFromSmarts("[OX1]=CN"),
    "alcohol": Chem.MolFromSmarts("[OX2H][#6]"),
    "ether": Chem.MolFromSmarts("[OD2]([#6])[#6]"),
    "amine": Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]"),
    "imine": Chem.MolFromSmarts("[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]"),
    "nitrile": Chem.MolFromSmarts("[NX1]#[CX2]"),
    "nitro": Chem.MolFromSmarts("[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"),
    "halide": Chem.MolFromSmarts("[F,Cl,Br,I]"),
})

class FuncGroupCounter:
    
    def __init__(self, fgs: dict = funcgroups_mols):
        self.fgs = fgs

    def count(self, mol: Chem.Mol, fgs: str | list) -> dict:
        """
        Counts the number of different functional groups in a molecule.

        Args:
            mol (Chem.Mol): RDKit molecule object to count functional groups from.
            fgs (str | list): Functional group(s) to count.

        Returns:
            dict: Dictionary containing the count of each functional group.
        """
        if isinstance(fgs, str):
            return len(mol.GetSubstructMatches(self.fgs[fgs]))
        
        count = {}
        for fg in fgs:
            count[fg] = len(mol.GetSubstructMatches(self.fgs[fg]))
        return count
    
    def count_all(self, mol: Chem.Mol) -> dict:
        """
        Counts the number of all functional groups in a molecule.

        Args:
            mol (Chem.Mol): RDKit molecule object to count functional groups from.

        Returns:
            dict: Dictionary containing the count of each functional group.
        """
        return {fg: len(mol.GetSubstructMatches(self.fgs[fg])) for fg in self.fgs.keys()}


def calculate_difficulty(mol: Chem.Mol, fgs: dict = funcgroups_mols, base_difficulty: int = 750) -> int:
    """
    Calculates the difficulty (Elo rating) of a molecule based on the functional groups present.
    """

    fg_counter = FuncGroupCounter(fgs)
    fg_counts = fg_counter.count_all(mol)
    difficulty = base_difficulty
    for fg, count in fg_counts.items():
        if fg in ["imine", "nitro"]:
            difficulty += count * 200
        else:
            difficulty += count * 100
    return difficulty