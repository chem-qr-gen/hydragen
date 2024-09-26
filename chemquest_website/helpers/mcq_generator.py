import random

from sqlalchemy import func
from rdkit import Chem, DataStructs
from rdkit.Chem import rdmolops, rdChemReactions, Descriptors, rdFingerprintGenerator

from chemquest_website import engine, meta
from chemquest_website.helpers.funcgroups import FuncGroupCounter, funcgroups_mols


# Functional group-based option generation
# TODO: add more of these

def generate_option_fluoride(mol: Chem.Mol) -> list[dict]:
    fluoride_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("F"))
    if not fluoride_matches:
        return []
    
    # replace halide with an alcohol
    alcohol_mols = rdmolops.ReplaceSubstructs(mol, Chem.MolFromSmarts("F"), Chem.MolFromSmiles("O"))
    options = [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["halide-alcohol"]}
         for m in alcohol_mols
    ]

    # replace halide with an amine
    amine_mols = rdmolops.ReplaceSubstructs(mol, Chem.MolFromSmarts("F"), Chem.MolFromSmiles("N"))
    options += [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["halide-amine"]}
         for m in amine_mols
    ]

    return options


def generate_option_chloride(mol: Chem.Mol) -> list[dict]:
    chloride_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("Cl"))
    if not chloride_matches:
        return []
    
    # replace chloride with a thiol
    thiol_mols = rdmolops.ReplaceSubstructs(mol, Chem.MolFromSmarts("Cl"), Chem.MolFromSmiles("S"))
    options = [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["halide-thiol"]}
         for m in thiol_mols
    ]

    return options


def generate_option_halide(mol: Chem.Mol) -> list[dict]:
    return generate_option_fluoride(mol) + generate_option_chloride(mol)


def generate_option_alcohol(mol: Chem.Mol) -> list[dict]:
    alcohol_matches = mol.GetSubstructMatches(funcgroups_mols.alcohol)

    if not alcohol_matches:
        return []
    
    # replace alcohol with a fluoride
    rxn1 = rdChemReactions.ReactionFromSmarts("[C:1][O:2]>>[C:1][F:2]")
    products1 = rxn1.RunReactants((mol,))
    options = [
        {"smiles": Chem.MolToSmiles(m[0]), "explanation": ["alcohol-halide"]}
         for m in products1
    ]

    # replace alcohol with an amine
    rxn2 = rdChemReactions.ReactionFromSmarts("[C:1][O:2]>>[C:1][N:2]")
    products2 = rxn2.RunReactants((mol,))
    options += [
        {"smiles": Chem.MolToSmiles(m[0]), "explanation": ["alcohol-amine"]}
         for m in products2
    ]

    # add methyl to form an ether
    rxn3 = rdChemReactions.ReactionFromSmarts("[C:1][O:2]>>[C:1][O:2]C")
    products3 = rxn3.RunReactants((mol,))
    options += [
        {"smiles": Chem.MolToSmiles(m[0]), "explanation": ["alcohol-ether"]}
         for m in products3
    ]

    return options


def generate_option_primary_o(mol: Chem.Mol) -> list[dict]:
    primary_o_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CH2]O"))

    if not primary_o_matches:
        return []
    
    # flip the alcohol or ether to form a different compound
    rxn = rdChemReactions.ReactionFromSmarts("[C:1][O:2]>>[O:2][C:1]")
    products = rxn.RunReactants((mol,))
    options = [
        {"smiles": Chem.MolToSmiles(p[0]), "explanation": ["primary-o-flip"]}
         for p in products
    ]

    return options


def generate_option_primary_amine(mol: Chem.Mol) -> list[dict]:
    amine_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[N;H2;X3]"))

    if not amine_matches:
        return []
    
    # replace amine with an alcohol
    alcohol_mols = rdmolops.ReplaceSubstructs(mol, Chem.MolFromSmarts("[N;H2;X3]"), Chem.MolFromSmiles("O"))
    options = [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["amine-alcohol"]}
         for m in alcohol_mols
    ]

    # replace amine with a halide
    halide_mols = rdmolops.ReplaceSubstructs(mol, Chem.MolFromSmarts("[N;H2;X3]"), Chem.MolFromSmiles("F"))
    options += [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["amine-halide"]}
         for m in halide_mols
    ]

    return options


def generate_option_secondary_amine(mol: Chem.Mol) -> list[dict]:
    amine_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[N;H1;X3]"))

    if not amine_matches:
        return []
    
    # flip the amine to form a different compound
    rxn = rdChemReactions.ReactionFromSmarts("[C:1][N:2][C:3]>>[N:2][C:1][C:3]")
    products = rxn.RunReactants((mol,))
    options = [
        {"smiles": Chem.MolToSmiles(p[0]), "explanation": ["secondary-amine-flip"]}
         for p in products
    ]

    # replace amine with an ether
    ether_mols = rdmolops.ReplaceSubstructs(mol, Chem.MolFromSmarts("[N;H1;X3]"), Chem.MolFromSmiles("O"))
    options += [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["secondary-amine-ether"]}
         for m in ether_mols
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

def generate_option_thiol(mol: Chem.Mol) -> list[dict]:
    thiol_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[#16X2H]"))

    if not thiol_matches:
        return []
    
    # replace thiol with a chloride
    chloride_mols = rdmolops.ReplaceSubstructs(mol, Chem.MolFromSmarts("[#16X2H]"), Chem.MolFromSmiles("Cl"))
    options = [
        {"smiles": Chem.MolToSmiles(m), "explanation": ["thiol-halide"]}
         for m in chloride_mols
    ]

    return options


# Overall MCQ option generation

def generate_mcq_options(mol: Chem.Mol, num_options: int = 3) -> list[dict]:
    options_raw = []

    options_raw += generate_option_halide(mol)
    options_raw += generate_option_alcohol(mol)
    options_raw += generate_option_primary_o(mol)
    options_raw += generate_option_primary_amine(mol)
    options_raw += generate_option_secondary_amine(mol)
    options_raw += generate_option_carbonyl(mol)
    options_raw += generate_option_ether(mol)
    options_raw += generate_option_ester(mol)
    options_raw += generate_option_carboxylic_acid(mol)
    options_raw += generate_option_nitrile(mol)
    options_raw += generate_option_thiol(mol)

    options = merge_duplicate_options(mol, options_raw)

    if len(options) >= num_options:
        sample = random.sample(options, num_options)
        return sample
    
    else:
        extra_options = generate_backup_options(mol, num_options)
        options_merged = merge_duplicate_options(mol, extra_options + options)
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


def merge_duplicate_options(original_mol: Chem.Mol, options_raw: list[dict]) -> list[dict]:
    """
    Merges duplicate options in the list of raw options. Combines the eligible explanations, and discards low-priority explanations (e.g. similarity-algo).
    """

    options = []
    for option in options_raw:
        # remove the option if it's the same as the original molecule
        if option["smiles"] == Chem.MolToSmiles(original_mol):
            continue

        # add option if it's not already in the list
        if option["smiles"] not in [o["smiles"] for o in options]:
            options.append(option)
        else:
            # skip if the option is a similarity algorithm option (an actual explanation is more valuable)
            if option["explanation"] == ["similarity-algo"]:
                continue
            
            # find the existing option and merge the explanations
            existing_option = next(o for o in options if o["smiles"] == option["smiles"])

            # replace the similarity algorithm explanation if the new option has a better explanation
            if existing_option["explanation"] == ["similarity-algo"]:
                existing_option["explanation"] = option["explanation"]
            elif option["explanation"] not in existing_option["explanation"]:
                existing_option["explanation"] += option["explanation"]

    return options