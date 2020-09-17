"""
Extract the ranges (contig, start, and end) of clusters detected by antiSMASH 5
in a given set of genome sequences, for comparison with PRISM.
"""

import json
import os
import pandas as pd
import re

# set working directory
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# set up I/O tuples
io = [# complete genomes
      [git_dir + "/data/genomes/complete-genomes-NCBI-11282018-derep.txt",
       git_dir + "/data/revisions/antismash-5/antismash_results",
       git_dir + "/data/analysis/genomes/clusters_antismash5.csv"],
       # metagenome-assembled genomes
      [git_dir + "/data/genomes/PRJNA348753-genomes.txt",
       git_dir + "/data/revisions/antismash-5/antismash_results",
       git_dir + "/data/analysis/MAGs/clusters_antismash5.csv"]]

for io_tuple in io:
    # read all genomes from file
    genomes = []
    derep_file = io_tuple[0]
    with open(derep_file) as f:
        lines = f.readlines()
        genomes = [line.strip('\n') for line in lines]

    # process each genome
    res = pd.DataFrame()
    for genome in genomes:
        print("  " + genome + " ...")
        # check output file
        genome_base = re.sub('\.fna', '', genome)
        genome_base = re.sub('\.gz', '', genome_base)
        output_file = io_tuple[1] + "/" + genome_base + '/' + genome_base + \
            '.json'
        if os.path.isfile(output_file):
            # extract all clusters from the JSON file
            f = open(output_file)
            lines = [line.strip('\n') for line in f.readlines()]
            string = ''.join(lines)
            root = json.loads(string)
            records = root['records']
            for record in records:
                contig = record['description']
                for feature in record['features']:
                    feature_type = feature['type']
                    if feature_type == 'region':
                        """
                        NAR 2019: "The regions in antiSMASH 5 correspond to the
                        entities called 'clusters' in antiSMASH 1 - 4 and now
                        constitute what is displayed on a page of the results
                        webpage."
                        """
                        cluster = feature
                        cluster_idx = cluster['qualifiers']['region_number']
                        location = cluster['location'].split(':')
                        start = re.sub('\[', '', location[0])
                        end = re.sub('\]', '', location[1])
                        cluster_type = '|'.join(cluster['qualifiers']['product'])
                        # append to results
                        row = pd.DataFrame({'genome': genome,
                                            'cluster': cluster_idx,
                                            'type': cluster_type,
                                            'contig': contig,
                                            'start': start,
                                            'end': end }, index=[0])
                        res = res.append(row)
        else:
            # genome output does not exist; append 'None'
            row = pd.DataFrame({'genome': genome, 'cluster': None, 'type': None,
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
