"""
Train final models for BGC activity using gradient boosting machines and all
features on the complete tungsten set.
"""

import os
import pandas as pd
from joblib import dump
from sklearn.svm import SVC
from sklearn.impute import SimpleImputer

# set directory 
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# read platinum and the associated activity matrix
paths = pd.read_table(git_dir + "/data/platinum/raw_paths.txt")
act = pd.read_csv(git_dir + "/data/platinum/activity_matrix.csv")
## drop unnecessary activities
act = act.loc[:, ['id', 'Bacteria', 'Fungi', 'Prokaryote', 'Virus', 'Cancer', 
                  'Immunomodulator']]
plat = pd.merge(paths, act, how='left', on='id')
## set cluster
clusters = [os.path.basename(file) for file in plat['fasta']]
plat = plat.assign(cluster=clusters)

# read fingerprints
fps = pd.read_csv(git_dir + "/data/platinum/PRISM_fingerprints_mean.csv.gz")
X = fps.set_index('id')
X = X.reindex(plat['id'])

# impute missing values
imputer = SimpleImputer(strategy='constant', fill_value=0)
imputer.fit(X)
X = imputer.transform(X)

# build models
targets = ['Bacteria', 'Fungi', 'Prokaryote', 'Virus', 'Cancer', 
           'Immunomodulator']
for target in targets:
    print("building models for target: " + target)
    y = plat[target].values
        
    # set up model
    svc = SVC(C=0.01, gamma=0.1, kernel='poly', degree=3, coef0=10.0,
              probability=True, random_state=0)
        
    # fit the models
    svc.fit(X, y)
    
    # save trained model
    output_file = git_dir + "/data/analysis/antibiotic/models/" + \
        target + ".joblib"
    # confirm directory exists
    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir): 
        os.makedirs(output_dir)
    
    dump(svc, output_file)

