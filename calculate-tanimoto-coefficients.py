"""
Calculate Tanimoto coefficients between chemical structure predictions from 
PRISM, NP.searcher, and antiSMASH and the true structure 
biosynthesized by each cluster.
"""

import itertools
import json
import os
import pandas as pd
from tqdm import tqdm

# set directory 
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# read true 'platinum' set structures
paths = pd.read_table(git_dir + "/data/platinum/raw_paths.txt")
mols = pd.read_table(git_dir + "/data/platinum/all_mols.txt")
plat = pd.merge(mols, paths, how='left', on='id')

# import functions
from functions import clean_mols, get_ecfp6_fingerprints, get_tanimoto
flatten = lambda x: list(itertools.chain.from_iterable(x))

# get unique inputs
clusters = plat.cluster.unique()

# for each cluster: 
tc = pd.DataFrame()
for cluster in tqdm(clusters): 
    # get real SMILES
    true_smiles = plat['smiles'][plat['cluster'] == cluster]
    # create molecules
    true_mols = clean_mols(true_smiles)
    # get fingerprints
    true_fps = get_ecfp6_fingerprints(true_mols)
    
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
            # get Tanimoto coefficients
            pred_mols = clean_mols(pred_smiles)
            pred_fps = get_ecfp6_fingerprints(pred_mols)
            tcs = get_tanimoto(true_fps, pred_fps)
            true_col = [[y for x in pred_smiles] for y in true_smiles]
            pred_col = [[x for x in pred_smiles] for y in true_smiles]
            # create data frame
            res = pd.DataFrame({'cluster': cluster, 
                                'true_smiles': flatten(true_col),
                                'pred_smiles': flatten(pred_col),
                                'Tc': tcs, 
                                'method': 'PRISM 4' })
            # append to master sheet
            tc = tc.append(res)
    else:
        # assign Tc of zero
        res = pd.DataFrame({'cluster': cluster, 
                            'true_smiles': true_smiles,
                            'pred_smiles': [None for smiles in true_smiles],
                            'Tc': [0 for smiles in true_smiles], 
                            'method': 'PRISM 4' })
        tc = tc.append(res)
    
    # test antiSMASH 5 predictions
    filename = cluster.rsplit('.', 1)[0]
    ## does the JSON file exist?
    as_dir = git_dir + "/data/predictions/antismash/" + filename
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
    
    # get Tanimoto coefficients
    if len(pred_smiles) > 0:
        pred_mols = clean_mols(pred_smiles)
        pred_fps = get_ecfp6_fingerprints(pred_mols)
        tcs = get_tanimoto(true_fps, pred_fps)
        true_col = [[y for x in pred_smiles] for y in true_smiles]
        pred_col = [[x for x in pred_smiles] for y in true_smiles]
        # create data frame
        res = pd.DataFrame({'cluster': cluster,
                            'true_smiles': flatten(true_col),
                            'pred_smiles': flatten(pred_col),
                            'Tc': tcs,
                            'method': 'antiSMASH 5' })
        # append to master sheet
        tc = tc.append(res)
    else:
        # assign Tc of zero
        res = pd.DataFrame({'cluster': cluster,
                            'true_smiles': true_smiles,
                            'pred_smiles': [None for smiles in true_smiles],
                            'Tc': [0 for smiles in true_smiles],
                            'method': 'antiSMASH 5' })
        tc = tc.append(res)
    
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
    
    if len(pred_smiles) > 0:
        if len(pred_smiles) > 100:
            pred_smiles = pred_smiles[0:100]
        pred_mols = clean_mols(pred_smiles)
        pred_fps = get_ecfp6_fingerprints(pred_mols)
        # get Tanimoto coefficients
        tcs = get_tanimoto(true_fps, pred_fps)
        true_col = [[y for x in pred_smiles] for y in true_smiles]
        pred_col = [[x for x in pred_smiles] for y in true_smiles]
        # create data frame
        res = pd.DataFrame({'cluster': cluster, 
                            'true_smiles': flatten(true_col),
                            'pred_smiles': flatten(pred_col),
                            'Tc': tcs, 
                            'method': 'NP.searcher' })
        # append to master sheet
        tc = tc.append(res)
    else:
        # assign Tc of zero
        res = pd.DataFrame({'cluster': cluster, 
                            'true_smiles': true_smiles,
                            'pred_smiles': [None for smiles in true_smiles],
                            'Tc': [0 for smiles in true_smiles], 
                            'method': 'NP.searcher' })
        tc = tc.append(res)

# write
output_file = git_dir + "/data/analysis/titanium/tanimoto_coefficients.csv.gz"
# confirm directory exists
output_dir = os.path.dirname(output_file)
if not os.path.isdir(output_dir): 
    os.mkdir(output_dir)

tc.to_csv(output_file, index=False, compression='gzip')

print('finished analysis')
