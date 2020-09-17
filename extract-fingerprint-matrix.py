"""
Calculate a feature matrix for machine learning on the 'tungsten' cluster set
by calculating 1024-bit hashed ECFP6 fingerprints for each PRISM predicted
structure, then averaging each bit over the entire set. 
"""

import json
import os
import pandas as pd
from rdkit.Chem import AllChem
from tqdm import tqdm

# set directory 
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# read platinum set clusters
paths = pd.read_table(git_dir + "/data/platinum/raw_paths.txt")
mols = pd.read_table(git_dir + "/data/platinum/all_mols.txt")
plat = pd.merge(mols, paths, how='left', on='id')

# import functions
from functions import clean_mol, get_bit_vector

# for each cluster: 
fps = pd.DataFrame()
for index, row in tqdm(plat.iterrows()):
    cluster_id = row['id']
    cluster = row['cluster']
    
    # generate bit fingerprints for PRISM predictions
    prism_dir = git_dir + "/data/predictions/prism"
    prism_file = prism_dir + "/" + cluster + ".json"
    if os.path.isfile(prism_file):
        # read all SMILES from JSON
        with open(prism_file) as f:
            root = json.load(f)
            prism = root['prism_results']
            prism_clusters = prism['clusters']
            pred_smiles = []
            for prism_cluster in prism_clusters:
                pathways = prism_cluster['biosynthetic_pathways']
                smiles = [pathway['smiles'] for pathway in pathways]
                pred_smiles.extend(smiles)
            # get fingerprints as bit vectors
            for smiles in pred_smiles:
                    mol = None
                    try:
                        mol = clean_mol(smiles)
                    except ValueError:
                        continue
                    fp = AllChem.GetMorganFingerprintAsBitVect(
                            mol, 3, nBits=1024)
                    vec = get_bit_vector(fp)
                    row = pd.concat([
                            pd.DataFrame({'id': cluster_id,
                                          'cluster': cluster,
                                          'smiles': smiles }, index=[0]),
                            pd.DataFrame(vec.reshape(-1, len(vec)))], axis=1)
                    fps = fps.append(row)

# write
output_file = git_dir + "/data/platinum/PRISM_fingerprints.csv.gz"
output_dir = os.path.dirname(output_file)
if not os.path.isdir(output_dir): 
    os.mkdir(output_dir)

fps.to_csv(output_file, index=False, compression='gzip')

# average fingerprints for each cluster 
avg = fps.drop(['cluster', 'smiles'], axis=1).groupby(['id']).mean()
# restore id column
avg = avg.reset_index()

# write
output_file = git_dir + "/data/platinum/PRISM_fingerprints_mean.csv.gz"
avg.to_csv(output_file, index=False, compression='gzip')
