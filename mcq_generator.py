import os.path
import pandas as pd
from FPSim2 import FPSim2Engine
from FPSim2.io import create_db_file

def init(db_filename = "db/molecules_original.csv", fp_filename = "db/fingerprints.h5") -> pd.DataFrame:
    '''Initialises the chemical fingerprints database if it doesn't already exist, and returns the pandas dataframe for later use.'''
    
    df = pd.read_csv(db_filename, dtype={"SMILES": "str", "mol_id": "int32"})

    if not os.path.isfile(fp_filename):
        print("mcq_generator: No fingerprint db found, generating one now.")
        create_db_file(df.values.tolist(), fp_filename, "Avalon")
        print("mcq_generator: Fingerprint db generated.")

    return df

def start_engine(fp_filename = "db/fingerprints.h5") -> FPSim2Engine:
    return FPSim2Engine(fp_filename)

def get_similar_compounds(df: pd.DataFrame, fpe: FPSim2Engine, input_smiles: str, threshold: float = 0.5, num_results: int = 10) -> list:
    results = fpe.similarity(input_smiles, threshold, n_workers = 4) # format: [mol_id, similarity]
    results_smiles = pd.DataFrame.from_records(
                         [(df.iloc[i[0] - 1]["SMILES"], i[0], i[1]) for i in results], 
                         columns = ["SMILES", "mol_id", "Similarity"]
                     )
    return results_smiles[1:num_results]
