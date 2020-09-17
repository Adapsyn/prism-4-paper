"""
Write the total number of clusters identified in each titanium set file for
PRISM, antiSMASH, and NP.searcher. 
"""

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

# get unique inputs
clusters = plat.cluster.unique()

# for each cluster: 
res = pd.DataFrame()
for cluster in tqdm(clusters): 
    # PRISM
    prism_dir = git_dir + "/data/predictions/prism"
    prism_file = prism_dir + "/" + cluster + ".json"
    n_prism = 0
    if os.path.isfile(prism_file):
        # read all SMILES from JSON
        with open(prism_file) as f:
            root = json.load(f)
            prism = root['prism_results']
            prism_clusters = prism['clusters']
            n_prism = len(prism_clusters)
    
    # antiSMASH
    filename = cluster.rsplit('.', 1)[0]
    as_dir = git_dir + "/data/predictions/antismash/" + filename
    as_file = as_dir + "/" + filename + ".json"
    n_antismash = 0
    if os.path.isfile(as_file):
        f = open(as_file)
        lines = [line.strip('\n') for line in f.readlines()]
        string = ''.join(lines)
        root = json.loads(string)
        records = root['records']
        for record in records:
            for feature in record['features']:
                feature_type = feature['type']
                if feature_type == 'region':
                    n_antismash += 1
    
    # NP.searcher
    np_dir = git_dir + "/data/predictions/npsearcher"
    np_files = os.listdir(np_dir)
    n_npsearcher = 0
    for np_file in np_files: 
        if np_file.startswith(cluster):
            n_npsearcher += 1
    
    # create row      
    row = pd.DataFrame({'file': cluster, 
                        'method': ['PRISM 4', 'antiSMASH 5', 'NP.searcher'],
                        'clusters': [n_prism, n_antismash, n_npsearcher] })
    res = res.append(row)

# write
output_file = git_dir + "/data/analysis/titanium/cluster_detection.csv"
# confirm directory exists
output_dir = os.path.dirname(output_file)
if not os.path.isdir(output_dir): 
    os.mkdir(output_dir)

res.to_csv(output_file, index=False)
