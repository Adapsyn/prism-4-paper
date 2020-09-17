"""
Extract cluster types and families from PRISM JSON files. 
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
types = pd.DataFrame()
for cluster in clusters: 
    # load PRISM JSON file
    prism_dir = git_dir + "/data/predictions/prism"
    prism_file = prism_dir + "/" + cluster + ".json"
    cluster_families = []
    cluster_types = []
    if os.path.isfile(prism_file):
        # read all SMILES from JSON
        with open(prism_file) as f:
            root = json.load(f)
            prism = root['prism_results']
            prism_clusters = prism['clusters']
            for prism_cluster in prism_clusters:
                cluster_families.extend(prism_cluster['family'])
                cluster_types.extend(prism_cluster['type'])
        
    # create data frame
    type_str = 'NA'
    if len(cluster_types) > 0:
         type_str = '|'.join(cluster_types)
    family_str = 'NA'
    if len(cluster_families) > 0:
         family_str = '|'.join(cluster_families)
    row = pd.DataFrame({'cluster': cluster, 'type': type_str,
                        'family': family_str }, index=[0])
    # append to master sheet
    types = types.append(row)

# write
output_file = git_dir + "/data/analysis/titanium/cluster_types.csv"
# confirm directory exists
output_dir = os.path.dirname(output_file)
if not os.path.isdir(output_dir): 
    os.mkdir(output_dir)

types.to_csv(output_file, index=False)
