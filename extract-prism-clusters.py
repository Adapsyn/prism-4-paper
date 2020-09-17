"""
Extract the ranges (contig, start, and end) of clusters detected by PRISM 
in a given set of genome sequences, for comparison with antiSMASH.
"""

import json
import os
import pandas as pd
from tqdm import tqdm

# set working directory
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# set up I/O tuples
io = [# complete genomes
      [git_dir + "/data/genomes/complete-genomes-NCBI-11282018-derep.txt",
       git_dir + "/data/genomes/prism_grid/output",
       git_dir + "/data/analysis/genomes/clusters_prism.csv"],
       # metagenome-assembled genomes
      [git_dir + "/data/genomes/PRJNA348753-genomes.txt",
       git_dir + "/data/MAGs/prism",
       git_dir + "/data/analysis/MAGs/clusters_prism.csv"]]

for io_tuple in io: 
    print("analyzing I/O pair: " + io_tuple[1])
    
    # read all genomes from file
    genomes = []
    derep_file = io_tuple[0]
    with open(derep_file) as f:
        lines = f.readlines()
        genomes = [line.strip('\n') for line in lines]
    
    # process each genome
    res = pd.DataFrame()
    for genome in tqdm(genomes):
        output_file = io_tuple[1] + "/" + genome.replace(".gz", "") + ".json"
        if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
            f = open(output_file)
            root = json.load(f)
            prism = root['prism_results']
            prism_clusters = prism['clusters']
            for cluster_idx, cluster in enumerate(prism_clusters):
                start = cluster['start']
                end = cluster['end']
                families = cluster['family']
                types = cluster['type']
                cluster_family = '|'.join(families)
                cluster_type = '|'.join(types)
                contig = cluster['contig_name']
                # append to results
                row = pd.DataFrame({'genome': genome, 'cluster': cluster_idx, 
                                    'type': cluster_type, 
                                    'family': cluster_family, 'contig': contig,
                                    'start': start, 'end': end }, index=[0])
                res = res.append(row)
        else:
            # genome output does not exist; append 'None'
            row = pd.DataFrame({'genome': genome, 'cluster': None, 
                                'type': None, 'family': None,
                                'contig': None, 'start': None, 'end': None }, 
                                index=[0])
            res = res.append(row)
    
    # write
    output_file = io_tuple[2] + ".gz"
    # confirm directory exists
    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir): 
        os.makedirs(output_dir)
    
    res.to_csv(output_file, index=False, compression='gzip')
