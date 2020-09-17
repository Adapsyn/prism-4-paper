"""
Calculate the Tc between PRISM and antiSMASH predicted natural products and
the entirety of the NP Atlas.
"""

import os
import pandas as pd
from tqdm import tqdm

# set working directory
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# import functions
from functions import clean_mols, get_ecfp6_fingerprints, get_tanimoto

# set up I/O pairs
io = [[git_dir + "/data/analysis/MAGs/fingerprints_intersect_random.csv.gz",
       git_dir + "/data/analysis/MAGs/intersect_Tc_NPAtlas.csv.gz"],
      [git_dir + "/data/analysis/genomes/fingerprints_intersect_random.csv.gz",
       git_dir + "/data/analysis/genomes/intersect_Tc_NPAtlas.csv.gz"]]

# read molecules from NP Atlas
atlas = pd.read_table(git_dir + "/data/smiles/NPatlas/np_atlas_2018_09.tsv.gz")
np_smiles = atlas['SMILES'].values
np_mols = clean_mols(np_smiles)
np_fps = get_ecfp6_fingerprints(np_mols)
 
for io_tuple in io: 
    input_file = io_tuple[0]
    output_file = io_tuple[1] 
    print("processing input file: " + os.path.basename(input_file))
        
    # read input 
    dat = pd.read_csv(input_file)
    
    # create results container
    res = pd.DataFrame()
    
    # process antiSMASH SMILES
    print(".. processing antiSMASH smiles ...")
    antismash_smiles = dat['smiles_as'].values
    antismash_mols = clean_mols(antismash_smiles)
    antismash_fps = get_ecfp6_fingerprints(antismash_mols)
    for i, query_fp in enumerate(tqdm(antismash_fps)): 
        tcs = get_tanimoto([query_fp], np_fps)
        rows = pd.DataFrame({'query_idx': i + 1,
                            'target_idx': list(range(1, len(tcs) + 1)),
                            'method': 'antiSMASH 5',
                            'Tc': tcs})
        res = res.append(rows)
        
    # process PRISM SMILES
    print(".. processing PRISM smiles ...")
    prism_smiles = dat['smiles_pr'].values
    prism_mols = clean_mols(prism_smiles)
    prism_fps = get_ecfp6_fingerprints(prism_mols)
    for i, query_fp in enumerate(tqdm(prism_fps)): 
        tcs = get_tanimoto([query_fp], np_fps)
        rows = pd.DataFrame({'query_idx': i + 1,
                            'target_idx': list(range(1, len(tcs) + 1)),
                            'method': 'PRISM 4',
                            'Tc': tcs})
        res = res.append(rows)
    
    # write     
    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir): 
        os.mkdir(output_dir)
    
    res.to_csv(output_file, index=False, compression='gzip')
