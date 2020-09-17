"""
Count functional groups in true and predicted products.
"""

import itertools
import json
import os
import pandas as pd
import random
import sys
from random import sample
from rdkit.Chem import RDConfig
from tqdm import tqdm

# import contrib modules
sys.path.append(os.path.join(RDConfig.RDContribDir, 'IFG'))
from ifg import identify_functional_groups

# set directory 
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# import functions
from functions import clean_mols
flatten = lambda x: list(itertools.chain.from_iterable(x))

# read true 'platinum' set structures
paths = pd.read_table(git_dir + "/data/platinum/raw_paths.txt")
mols = pd.read_table(git_dir + "/data/platinum/all_mols.txt")
plat = pd.merge(mols, paths, how='left', on='id')

# get unique inputs
clusters = plat.cluster.unique()

# for each cluster: 
res = pd.DataFrame()
random.seed(0)
for cluster in tqdm(clusters): 
    print(cluster)
    # get real SMILES
    true_smiles = plat['smiles'][plat['cluster'] == cluster]
    # create molecules
    true_mols = clean_mols(true_smiles)
    # randomly sample one
    true_mol = sample(true_mols, 1)[0]
    # get functional groups
    fgs = identify_functional_groups(true_mol)
    row = pd.DataFrame({'cluster': cluster, 
                        'type': 'True',
                        'functional_group': [str(fg) for fg in fgs]})
    res = res.append(row)
    
    # test PRISM predictions
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
            
            # get functional groups
            pred_mols = [mol for mol in clean_mols(pred_smiles) if \
                         mol is not None]
            if len(pred_mols) > 0:
                # randomly sample one
                pred_mol = sample(pred_mols, 1)[0]
                fgs = identify_functional_groups(pred_mol)
                row = pd.DataFrame({'cluster': cluster, 
                                    'type': 'PRISM 4',
                                    'functional_group': [str(fg) for fg in \
                                                         fgs]})
                res = res.append(row)
    
    # test antiSMASH 5 predictions
    filename = cluster.rsplit('.', 1)[0]
    ## does the JSON file exist?
    as_dir = git_dir + "/data/predictions/antismash/platinum/" + filename
    as_file = as_dir + "/" + filename + ".json"
    ## get all SMILES
    pred_smiles = []
    if os.path.isfile(as_file):
        # extract all clusters from the JSON file
        f = open(as_file)
        lines = [line.strip('\n') for line in f.readlines()]
        string = ''.join(lines)
        root = json.loads(string)
        records = root['records']
        for record in records:
            modules = record['modules']
            for module_name in modules.keys():
                module = modules[module_name]
                if 'region_predictions' in module.keys():
                    preds = module['region_predictions']
                    for cluster_idx in preds.keys():
                        clust = preds[cluster_idx]
                        for prediction in clust:
                            smiles = prediction['smiles']
                            pred_smiles.append(smiles)
    
    # get functional groups
    pred_mols = [mol for mol in clean_mols(pred_smiles) if mol is not None]
    if len(pred_mols) > 0:
        # randomly sample one
        pred_mol = sample(pred_mols, 1)[0]
        fgs = identify_functional_groups(pred_mol)
        row = pd.DataFrame({'cluster': cluster, 
                            'type': 'antiSMASH 5',
                            'functional_group': [str(fg) for fg in fgs]})
        res = res.append(row)
    
    # test NP.searcher predictions
    npsearcher_dir = git_dir + "/data/predictions/npsearcher"
    npsearcher_files = os.listdir(npsearcher_dir)
    pred_smiles = []
    for npsearcher_file in npsearcher_files:
        if npsearcher_file.startswith(cluster): 
            with open(npsearcher_dir + "/" + npsearcher_file) as f:
                lines = f.readlines()
                if len(lines) > 1:
                    smiles_line = lines[1]
                    pred_smiles.extend([smiles.replace('[X]', '[*]') for \
                                        smiles in smiles_line.split('.')])
    
    # get functional groups
    if len(pred_smiles) > 0:
        if len(pred_smiles) > 100:
            pred_smiles = pred_smiles[0:100]
    pred_mols = [mol for mol in clean_mols(pred_smiles) if mol is not None]
    if len(pred_mols) > 0:
        # randomly sample one
        pred_mol = sample(pred_mols, 1)[0]
        fgs = identify_functional_groups(pred_mol)
        row = pd.DataFrame({'cluster': cluster, 
                            'type': 'NP.searcher',
                            'functional_group': [str(fg) for fg in fgs]})
        res = res.append(row)
    
# write
output_file = git_dir + "/data/analysis/titanium/functional_groups.csv.gz"
# confirm directory exists
output_dir = os.path.dirname(output_file)
if not os.path.isdir(output_dir): 
    os.mkdir(output_dir)

res.to_csv(output_file, index=False, compression='gzip')
print('finished analysis')
