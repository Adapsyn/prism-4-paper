"""
Extract 1024-bit ECFP6 fingerprints for known natural product structures
from the NP Atlas.
"""

import os
import pandas as pd
from rdkit.Chem import AllChem
from tqdm import tqdm

# set working directory
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# import functions
from functions import clean_mol, get_bit_vector

# read data
atlas = pd.read_table(git_dir + "/data/smiles/NPatlas/np_atlas_2018_09.tsv.gz")
atlas = atlas.loc[:, ['NPAID', 'SMILES']]

# get fingerprints
res = pd.DataFrame()
for idx, row in tqdm(atlas.iterrows()):
    mol_id = row['NPAID']
    smiles = row['SMILES']
    mol = clean_mol(smiles)
    if mol is None:
        continue 
    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, 3, nBits=1024)
    vec = get_bit_vector(fp)
    row = pd.concat([
        pd.DataFrame({ 'id': mol_id, 'smiles': smiles }, index=[0]),
        pd.DataFrame(vec.reshape(-1, len(vec)))], 
                 axis=1)
    res = res.append(row)

# write
res.to_csv(git_dir + '/data/smiles/NPatlas/np_atlas_fingerprints.csv.gz',
           index=False, compression='gzip')
