import itertools
import os
import pandas as pd
from tqdm import tqdm

# set directory
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# import functions
from functions import clean_mols, get_ecfp6_fingerprints, get_tanimoto
flatten = lambda x: list(itertools.chain.from_iterable(x))

# read true and predicted structures
tc = pd.read_csv(git_dir + '/data/analysis/titanium/tanimoto_coefficients.csv.gz')
## prefilter missing Tcs
tc = tc.dropna(subset=['pred_smiles'])
# get unique clusters and methods
clusters = tc['cluster'].unique()

# set up output data frame
mism = pd.DataFrame()

# iterate over clusters
for cluster in tqdm(clusters):
    # get the query structures
    pred_smiles = tc[(tc['method'] == 'PRISM 4') & \
                     (tc['cluster'] == cluster)].pred_smiles.unique()
    true_smiles = tc[(tc['method'] == 'PRISM 4') & \
                     (tc['cluster'] != cluster)].true_smiles.unique()

    # get Tanimoto coefficients
    true_mols = clean_mols(true_smiles)
    true_fps = get_ecfp6_fingerprints(true_mols)
    pred_mols = clean_mols(pred_smiles)
    pred_fps = get_ecfp6_fingerprints(pred_mols)
    tcs = get_tanimoto(true_fps, pred_fps)
    true_col = [[y for x in pred_smiles] for y in true_smiles]
    pred_col = [[x for x in pred_smiles] for y in true_smiles]
    # create data frame
    res = pd.DataFrame({'cluster': cluster,
                        'true_smiles': flatten(true_col),
                        'pred_smiles': flatten(pred_col),
                        'Tc': tcs })

    # append to master sheet
    mism = mism.append(res)

# write
output_file = git_dir + "/data/analysis/titanium/tanimoto_coefficients_mismatched.csv.gz"
# confirm directory exists
output_dir = os.path.dirname(output_file)
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

mism.to_csv(output_file, index=False, compression='gzip')

print('finished analysis')
