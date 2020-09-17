"""
Extract 1024-bit ECFP6 fingerprints for antiSMASH predicted structures
for clusters in a large set of complete bacterial genomes.
"""

import json
import os
import pandas as pd
import re
from rdkit.Chem import AllChem

# set working directory
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# import functions
from functions import clean_mol, get_bit_vector

# set up I/O tuples
io = [# complete genomes
      [git_dir + "/data/genomes/complete-genomes-NCBI-11282018-derep.txt",
       git_dir + "/data/genomes/antismash",
       git_dir + "/data/analysis/genomes/fingerprints_antismash.csv"],
       # metagenome-assembled genomes
      [git_dir + "/data/genomes/PRJNA348753-genomes.txt",
       git_dir + "/data/MAGs/antismash",
       git_dir + "/data/analysis/MAGs/fingerprints_antismash.csv"]]

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
                modules = record['modules']
                for module_name in modules.keys():
                    module = modules[module_name]
                    if 'region_predictions' in module.keys():
                        preds = module['region_predictions']
                        for cluster_idx in preds.keys():
                            cluster = preds[cluster_idx]
                            for prediction in cluster:
                                smiles = prediction['smiles']
                                try:
                                    mol = clean_mol(smiles)
                                    if mol is None:
                                        continue
                                    fp = AllChem.GetMorganFingerprintAsBitVect(
                                            mol, 3, nBits=1024)
                                    vec = get_bit_vector(fp)
                                    row = pd.concat([
                                            pd.DataFrame({'genome': genome,
                                                          'contig': contig,
                                                          'cluster': cluster_idx,
                                                          'smiles': smiles },
                                                         index=[0]),
                                            pd.DataFrame(vec.reshape(-1, len(vec)))],
                                                     axis=1)
                                    res = res.append(row)
                                except ValueError:
                                    continue

    # write
    output_file = io_tuple[2] + ".gz"
    # confirm directory exists
    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    res.to_csv(output_file, index=False, compression='gzip')
