"""
Use trained models to predict biological activity for encoded natural products
in large genomic collections.
"""

import os
import pandas as pd
from joblib import load

# set directory 
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# read fingerprints
fp1 = pd.read_csv(git_dir + "/data/analysis/genomes/fingerprints_prism.csv.gz")
fp2 = pd.read_csv(git_dir + "/data/analysis/MAGs/fingerprints_prism.csv.gz")

# average fingerprints for each cluster
fp1['cluster'] = fp1['cluster'].apply(str)
X1 = fp1.assign(cluster_id=fp1['genome'] + '|' + fp1['cluster']).\
    drop(['genome', 'cluster', 'smiles'], axis=1).\
    groupby(['cluster_id']).\
    mean()
fp2['cluster'] = fp2['cluster'].apply(str)
X2 = fp2.assign(cluster_id=fp2['genome'] + '|' + fp2['cluster']).\
    drop(['genome', 'cluster', 'smiles'], axis=1).\
    groupby(['cluster_id']).\
    mean()

# run trained models one at a time 
targets = ['Bacteria', 'Fungi', 'Prokaryote', 'Virus', 'Cancer', 
           'Immunomodulator']
for target in targets:
    print("predicting target: " + target)
    model_file = git_dir + "/data/analysis/antibiotic/models/" + \
        target + ".joblib"
    model = load(model_file)
    
    # predict
    y1 = model.predict_proba(X1)[:, 1]
    y2 = model.predict_proba(X2)[:, 1]
    
    # write data frames
    df1 = pd.DataFrame({'cluster_id': X1.index, 'score': y1})
    df2 = pd.DataFrame({'cluster_id': X2.index, 'score': y2})
    output1 = git_dir + "/data/analysis/genomes/activity/" + target + ".csv"
    output2 = git_dir + "/data/analysis/MAGs/activity/" + target + ".csv"
    
    # confirm directories exist
    output_dir1 = os.path.dirname(output1)
    if not os.path.isdir(output_dir1): 
        os.makedirs(output_dir1)
    
    output_dir2 = os.path.dirname(output2)
    if not os.path.isdir(output_dir2): 
        os.makedirs(output_dir2)
        
    df1.to_csv(output1, index=False)
    df2.to_csv(output2, index=False)
