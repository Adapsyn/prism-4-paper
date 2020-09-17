"""
Extract 1024-bit ECFP6 fingerprints for PRISM predicted structures
for clusters in a large set of complete bacterial genomes.
"""

import json
import os
import pandas as pd
from rdkit.Chem import AllChem
from tqdm import tqdm

# set directory 
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# import functions
from functions import clean_mol, get_bit_vector

# set up I/O tuples
io = [# complete genomes
      [git_dir + "/data/genomes/complete-genomes-NCBI-11282018-derep.txt",
       git_dir + "/data/genomes/prism_grid/output",
       git_dir + "/data/analysis/genomes/fingerprints_prism.csv"],
       # metagenome-assembled genomes
      [git_dir + "/data/genomes/PRJNA348753-genomes.txt",
       git_dir + "/data/MAGs/prism_lite/structure_pred",
       git_dir + "/data/analysis/MAGs/fingerprints_prism.csv"]]
## run #2: just get complete genomes
io = [io[0]]

for io_tuple in io: 
    print("analyzing I/O pair: " + io_tuple[1])
    
    # read all genomes from file
    genomes = []
    derep_file = io_tuple[0]
    with open(derep_file) as f:
        lines = f.readlines()
        genomes = [line.strip('\n').replace(".gz", "") for line in lines]
    
    # process each genome
    res = pd.DataFrame()
    for genome in tqdm(genomes):
        print(genome + "...")
        output_file = io_tuple[1] + "/" + genome + ".json"
        if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
            f = open(output_file)
            root = json.load(f)
            prism = root['prism_results']
            prism_clusters = prism['clusters']
            for cluster_idx, cluster in enumerate(prism_clusters):
                pathways = cluster['biosynthetic_pathways']
                pred_smiles = [pathway['smiles'] for pathway in pathways]
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
                            pd.DataFrame({'genome': genome, 
                                          'cluster': cluster_idx, 
                                          'smiles': smiles }, index=[0]),
                            pd.DataFrame(vec.reshape(-1, len(vec)))], axis=1)
                    res = res.append(row)  
    
    # write
    output_file = io_tuple[2] + ".gz"
    # confirm directory exists
    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir): 
        os.makedirs(output_dir)
    
    res.to_csv(output_file, index=False, compression='gzip')
