"""
Compare precision-recall curves for predicting BGC activity with different
sets of features for the tungsten set compounds, using 10-fold CV.
"""

import os
import pandas as pd
import sys
from sklearn.impute import SimpleImputer
from sklearn.model_selection import KFold
from sklearn.svm import SVC
from tqdm import tqdm

# get index from CLI
target_idx = int(sys.argv[1]) - 1

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

# read features
## Pfam domains
pfam = pd.read_csv(git_dir + "/data/platinum/Pfam_domains.csv")
## fingerprints
fps = pd.read_csv(git_dir + "/data/platinum/PRISM_fingerprints_mean.csv.gz")

# define features
features = {'pfam': pfam, 'fingerprints': fps}

# build models
targets = ['Bacteria', 'Fungi', 'Virus', 'Cancer', 'Immunomodulator']
target = targets[target_idx]
print("building models for target: " + target)
pred = pd.DataFrame()

# also iterate over feature combinations
for feature in tqdm(features.keys()):
    X = features[feature]
    X = X.set_index('id')
    X = X.reindex(plat['id'])
    y = plat[target].values

    # impute missing values
    imputer = SimpleImputer(strategy='constant', fill_value=0)
    imputer.fit(X)
    X = imputer.transform(X)

    # set up model
    model = SVC(C=0.01, gamma=0.1, kernel='poly', degree=3, coef0=10.0,
                probability=True, random_state=0)

    # fit the model
    cv = KFold(n_splits=10)
    for train_idx, test_idx in cv.split(X):
        X_train = X[train_idx]
        y_train = y[train_idx]
        X_test = X[test_idx]
        y_test = y[test_idx]
        # fit model and predict on y
        model.fit(X_train, y_train)
        probs = model.predict_proba(X_test)[:, 1]
        # append
        row = pd.DataFrame({'feature': feature,
                            'target': target,
                            'split': test_idx,
                            'label': y_test,
                            'score': probs})
        pred = pred.append(row)

# write results
output = git_dir + "/data/analysis/antibiotic/classifier_10foldCV_" + target \
    + ".csv"
# confirm directory exists
output_dir = os.path.dirname(output)
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

pred.to_csv(output, index=False)
