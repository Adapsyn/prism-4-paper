"""
Extract the feature matrices expected by trained models from PRISM JSON files
corresponding to large genomic collections.
"""

import json
import os
import pandas as pd
from tqdm import tqdm

# set directory 
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# read and concatenate features (domains/substrates)
dom = pd.read_csv(git_dir + "/data/platinum/PRISM_domains.csv")
sub = pd.read_csv(git_dir + "/data/platinum/PRISM_substrates.csv")
X = pd.merge(dom, sub, on='id').set_index('id')

# set up I/O tuples
io = [# complete genomes
      [git_dir + "/data/genomes/complete-genomes-NCBI-11282018-derep.txt",
       git_dir + "/data/genomes/prism_grid/output",
       git_dir + "/data/analysis/genomes/features_prism.csv"],
       # metagenome-assembled genomes
      [git_dir + "/data/genomes/PRJNA348753-genomes.txt",
       git_dir + "/data/MAGs/prism_lite/structure_pred",
       git_dir + "/data/analysis/MAGs/features_prism.csv"]]

for io_tuple in io: 
    print("analyzing I/O pair: " + io_tuple[1])
    
    # read all genomes from file
    genomes = []
    derep_file = io_tuple[0]
    with open(derep_file) as f:
        lines = f.readlines()
        genomes = [line.strip('\n').replace(".gz", "") for line in lines]
    
    # create containers to avoid pd.concat
    all_domains = {}
    all_substrates = {}
    
    # process each genome
    for genome in tqdm(genomes):
        print(genome + "...")
        output_file = io_tuple[1] + "/" + genome + ".json"
        if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
            f = open(output_file)
            root = json.load(f)
            prism = root['prism_results']
            prism_clusters = prism['clusters']
            for cluster_idx, cluster in enumerate(prism_clusters):
                # get domains and substrates
                domains = []
                substrates = []
                orfs = [orf for orf in cluster['orfs'] if \
                        len(orf['domains']) > 0]
                for orf in orfs: 
                    for domain in orf['domains']: 
                        domain_name = domain['family'] + ":" + \
                            domain['full_name']
                        domains.append(domain_name)
                        if len(domain['substrates']) > 0:
                            substrate_name = domain['full_name'] + ":" + \
                                domain['substrates'][0]['name']
                            substrates.append(substrate_name)
                # create row
                cluster_name = genome + "|" + str(cluster_idx)
                all_domains[cluster_name] = domains
                all_substrates[cluster_name] = substrates

    # convert to data frame
    cluster_genomes = [x.split("|")[0] for x in all_domains.keys()]
    cluster_idxs = [x.split("|")[1] for x in all_domains.keys()]
    feat = pd.DataFrame({'genome': cluster_genomes,
                         'cluster': cluster_idxs})
    
    # add columns one at a time
    for col_name in list(X):
        col_value = [all_domains[key].count(col_name) for key in \
                     all_domains.keys()]
        feat[col_name] = col_value
    
    # write
    output_file = io_tuple[2] + ".gz"
    # confirm directory exists
    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir): 
        os.makedirs(output_dir)
    
    feat.to_csv(output_file, index=False, compression='gzip')
