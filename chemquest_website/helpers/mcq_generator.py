import random

from sqlalchemy import func
from rdkit import Chem, DataStructs
from rdkit.Chem import rdmolops, rdChemReactions, Descriptors, rdFingerprintGenerator

from chemquest_website import engine, meta
from chemquest_website.helpers.funcgroups import FuncGroupCounter, funcgroups_mols


# Functional group-based option generation
# TODO: add more of these
# converting amine/alcohol to halide seems like an obvious choice
# converting primary alcohols and amines to methyl secondary alcohols/amines (flip carbon and O/N)

def generate_option_halide(mol: Chem.Mol) -> list[dict]:
    halide_matches = mol.GetSubstructMatches(funcgroups_mols.halide)

    if not halide_matches:
        return []
    
    # replace halide with an alcohol
    alcohol_mols = rdmolops.ReplaceSubstructs(mol, funcgroups_mols.halide, Chem.MolFromSmiles("O"))
    options = [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["halide-alcohol"]}
         for m in alcohol_mols
    ]

    # replace halide with an amine
    amine_mols = rdmolops.ReplaceSubstructs(mol, funcgroups_mols.halide, Chem.MolFromSmiles("N"))
    options += [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["halide-amine"]}
         for m in amine_mols
    ]

    return options


def generate_option_carbonyl(mol: Chem.Mol) -> list[dict]:
    carbonyl_matches = mol.GetSubstructMatches(funcgroups_mols.ketone) + mol.GetSubstructMatches(funcgroups_mols.aldehyde)

    if not carbonyl_matches:
        return []
    
    # replace carbonyl with a methyl group
    alkyl_mols = rdmolops.ReplaceSubstructs(mol, Chem.MolFromSmarts("C=O"), Chem.MolFromSmiles("C(C)"))
    options = [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["carbonyl-methyl"]}
         for m in alkyl_mols
    ]

    return options


def generate_option_ether(mol: Chem.Mol) -> list[dict]:
    ether_matches = mol.GetSubstructMatches(funcgroups_mols.ether)

    if not ether_matches:
        return []
    
    # replace ether with an alcohol
    rxn = rdChemReactions.ReactionFromSmarts("[C:1][O:2][C:3]>>[C:1]([O:2])[C:3]")
    products = rxn.RunReactants((mol,))
    options = [
        {"smiles": Chem.MolToSmiles(p[0]), "explanation": ["ether-alcohol"]}
         for p in products
    ]

    return options


def generate_option_ester(mol: Chem.Mol) -> list[dict]:
    ester_matches = mol.GetSubstructMatches(funcgroups_mols.ester)

    if not ester_matches:
        return []
    
    # flip ester
    rxn1 = rdChemReactions.ReactionFromSmarts("[C:1](=[O:2])[O:3][C:4]>>[C:1][O:2][C:4](=[O:3])")
    products = rxn1.RunReactants((mol,))
    options = [
        {"smiles": Chem.MolToSmiles(p[0]), "explanation": ["ester-flip"]}
         for p in products
    ]

    # replace ester with a ketone + alcohol
    rxn2 = rdChemReactions.ReactionFromSmarts("[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[C:4][O:3]")
    products = rxn2.RunReactants((mol,))
    options += [
        {"smiles": Chem.MolToSmiles(p[0]), "explanation": ["ester-ketone-alcohol"]}
         for p in products
    ]

    # replace ester with an amide
    rxn3 = rdChemReactions.ReactionFromSmarts("[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[N][C:4]")
    products = rxn3.RunReactants((mol,))
    options += [
        {"smiles": Chem.MolToSmiles(p[0]), "explanation": ["ester-amide"]}
         for p in products
    ]

    return options


def generate_option_carboxylic_acid(mol: Chem.Mol) -> list[dict]:
    acid_matches = mol.GetSubstructMatches(funcgroups_mols.carboxylic_acid)

    if not acid_matches:
        return []
    
    # replace acid with an ester
    methyl_ester_mols = rdmolops.ReplaceSubstructs(mol, funcgroups_mols.carboxylic_acid, Chem.MolFromSmiles("C(=O)OC"))
    options = [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["acid-ester"]}
         for m in methyl_ester_mols
    ]
    ethyl_ester_mols = rdmolops.ReplaceSubstructs(mol, funcgroups_mols.carboxylic_acid, Chem.MolFromSmiles("C(=O)OCC"))
    options += [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["acid-ester"]}
         for m in ethyl_ester_mols
    ]

    # replace acid with an amide
    amide_mols = rdmolops.ReplaceSubstructs(mol, funcgroups_mols.carboxylic_acid, Chem.MolFromSmiles("C(=O)N"))
    options += [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["acid-amide"]}
         for m in amide_mols
    ]

    return options

def generate_option_nitrile(mol: Chem.Mol) -> list[dict]:
    nitrile_matches = mol.GetSubstructMatches(funcgroups_mols.nitrile)

    if not nitrile_matches:
        return []
    
    # replace nitrile with a terminal alkyne
    alkyne_mols = rdmolops.ReplaceSubstructs(mol, funcgroups_mols.nitrile, Chem.MolFromSmiles("C#C"))
    options = [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["nitrile-alkyne"]}
         for m in alkyne_mols
    ]

    return options



# Overall MCQ option generation

def generate_mcq_options(mol: Chem.Mol, num_options: int = 3) -> list[dict]:
    options_raw = []

    options_raw += generate_option_halide(mol)
    options_raw += generate_option_carbonyl(mol)
    options_raw += generate_option_ether(mol)
    options_raw += generate_option_ester(mol)
    options_raw += generate_option_carboxylic_acid(mol)

    options = merge_duplicate_options(options_raw)

    if len(options) >= num_options:
        sample = random.sample(options, num_options)
        return sample
    
    else:
        extra_options = generate_backup_options(mol, num_options)
        options_merged = merge_duplicate_options(extra_options + options)
        sample = random.sample(options_merged, num_options)
        return sample
    

def generate_backup_options(mol: Chem.Mol, num_options: int) -> list[dict]:
    """
    Generates MCQ options for mass spectrometry questions based on molecular weight and a similarity algorithm.
    This is the "backup" for when the functional group-based generation does not provide enough options.
    Does not allow for explanations of wrong answers.
    """
    source_mol_weight = Descriptors.MolWt(mol)

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
    source_fp = fpgen.GetFingerprint(mol)
    fps_to_compare = [fpgen.GetFingerprint(m) for m in mols_to_compare]

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
    # remove the similarity values
    top_options = [
        {"smiles": smiles, "explanation": ["similarity-algo"]}
         for smiles, _ in top_options]
    return top_options


def merge_duplicate_options(options_raw: list[dict]) -> list[dict]:
    """
    Merges duplicate options in the list of raw options. Combines the eligible explanations, and discards low-priority explanations (e.g. similarity-algo).
    """

    options = []
    for option in options_raw:
        # add option if it's not already in the list
        if option["smiles"] not in [o["smiles"] for o in options]:
            options.append(option)
        else:
            # skip if the option is a similarity algorithm option (an actual explanation is more valuable)
            if option["explanation"] == ["similarity-algo"]:
                continue
            
            # find the existing option and merge the explanations
            existing_option = next(o for o in options if o["smiles"] == option["smiles"])
            existing_option["explanation"] += option["explanation"]

    return options