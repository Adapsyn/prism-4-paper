"""
Extract combinatorial data from PRISM JSON files. 
"""

import json
import os
import pandas as pd

# set directory 
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# read true 'platinum' set structures
paths = pd.read_table(git_dir + "/data/platinum/raw_paths.txt")
mols = pd.read_table(git_dir + "/data/platinum/all_mols.txt")
plat = pd.merge(mols, paths, how='left', on='id')

# get unique inputs
clusters = plat.cluster.unique()

# for each cluster: 
combns = pd.DataFrame()
for cluster in clusters: 
    # load PRISM JSON file
    prism_dir = git_dir + "/data/predictions/prism"
    prism_file = prism_dir + "/" + cluster + ".json"
    # create lists for all clusters
    orf_permns = []
    n_graphs = []
    n_propeptides = []
    n_sugars = []
    n_combinatorial_plans = []
    if os.path.isfile(prism_file):
        # read all SMILES from JSON
        with open(prism_file) as f:
            root = json.load(f)
            prism = root['prism_results']
            prism_clusters = prism['clusters']
            for cluster_idx, prism_cluster in enumerate(prism_clusters):
                combn_dat = prism_cluster['combinatorial_data']
                orf_permns.append(combn_dat['num_orf_permutations'])
                n_graphs.append(combn_dat['num_graphs'])
                n_propeptides.append(combn_dat['num_propeptides'])
                n_sugars.append(combn_dat['num_sugars'])
                n_combinatorial_plans.append(
                        combn_dat['num_combinatorial_plans'])
    
    # create data frame
    if len(n_combinatorial_plans) > 0:
        row = pd.DataFrame({'cluster': cluster,
                            'n_combinatorial_plans': n_combinatorial_plans,
                            'n_graphs': n_graphs,
                            'n_orf_permutations': orf_permns,
                            'n_sugars': n_sugars,
                            'n_propeptides': n_propeptides })
        combns = combns.append(row)
    else: 
        row = pd.DataFrame({'cluster': cluster,
                            'n_combinatorial_plans': None,
                            'n_graphs': None,
                            'n_orf_permutations': None,
                            'n_sugars': None,
                            'n_propeptides': None }, index=[0])
        combns = combns.append(row)

# write
output_file = git_dir + "/data/analysis/titanium/combinatorial_data.csv"
# confirm directory exists
output_dir = os.path.dirname(output_file)
if not os.path.isdir(output_dir): 
    os.mkdir(output_dir)

combns.to_csv(output_file, index=False)
