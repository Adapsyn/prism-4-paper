"""
Calculate a set of chemical descriptors for two lists of predicted 
natural products.
"""

import os
import pandas as pd
from collections import Counter
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Lipinski
from rdkit.Chem.AllChem import CalcNumAtomStereoCenters
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.GraphDescriptors import BertzCT

# set working directory
git_dir = os.path.expanduser("~/git/prism-4-paper")
os.chdir(git_dir)

# import functions
from functions import clean_mol, in_Ro5
from scores import readNPModel, calculateNPScore, readSAModel, calculateSAScore

# read NP-likeness and synthetic accessibility models
np_mod = readNPModel(git_dir + "/data/rdkit/publicnp.model.gz")
sa_mod = readSAModel(git_dir + "/data/rdkit/fpscores.pkl.gz")

# define function to calculate all metrics
def calculate_metrics(mol):
    # calculate chemical descriptors
    ## % of sp3 carbons
    pct_sp3 = Lipinski.FractionCSP3(mol)
    ## H bond donors/acceptors
    h_acceptor = Lipinski.NumHAcceptors(mol)
    h_donor = Lipinski.NumHDonors(mol)
    ## number of rotable bonds
    n_bonds = mol.GetNumBonds()
    if n_bonds > 0:
        rot_bonds = Lipinski.NumRotatableBonds(mol) / n_bonds
    else: 
        rot_bonds = 0
    ## number of rings, aromatic and aliphatic
    n_rings = Lipinski.RingCount(mol)
    n_rings_ali = Lipinski.NumAliphaticRings(mol)
    n_rings_aro = Lipinski.NumAromaticRings(mol)
    ## number of stereocentres
    Chem.AssignStereochemistry(mol)
    n_stereo = CalcNumAtomStereoCenters(mol)
    ## polarity
    tpsa = Chem.CalcTPSA(mol)
    ## hydrophobicity
    logP = Descriptors.MolLogP(mol)
    ## molecular weight
    mw = Descriptors.MolWt(mol)
    ## in Lipinski space?
    Ro5 = in_Ro5(mol)
    ## % heteroatoms
    n_atoms = len(mol.GetAtoms())
    pct_hetero = Lipinski.NumHeteroatoms(mol) / n_atoms
    ## number of each atom
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    atom_counts = Counter(symbols)
    ## Murcko scaffolds
    murcko = Chem.MolToSmiles(
            MurckoScaffold.GetScaffoldForMol(mol))
    ## NP-likeness 
    try: 
        np_score = calculateNPScore(mol, np_mod)
    except ValueError:
        np_score = None
    ## synthetic accessibility 
    try:
        sa_score = calculateSAScore(mol, sa_mod)
    except ValueError:
        sa_score = None
    ## topological complexity
    bertz_idx = BertzCT(mol)
    # create dict
    metrics = {'% sp3 carbons': pct_sp3, 
               'H bond acceptors': h_acceptor,
               'H bond donors': h_donor,
               '% rotatable bonds': rot_bonds,
               'Rings': n_rings,
               'Rings, aliphatic': n_rings_ali,
               'Rings, aromatic': n_rings_aro,
               'Stereocentres': n_stereo,
               'Topological polar surface area': tpsa,
               'LogP': logP,
               'Molecular weight': mw,
               'Lipinski rule of 5': Ro5,
               '% heteroatoms': pct_hetero,
               'Murcko scaffold': murcko,
               'NP-likeness score': np_score,
               'Synthetic accessibility score': sa_score,
               'Bertz topological complexity': bertz_idx}
    # append atom counts
    for key in atom_counts.keys(): 
        metrics['Atoms with symbol ' + key] = atom_counts[key]
    return(metrics)

# set up I/O pairs
io = [[git_dir + "/data/analysis/MAGs/fingerprints_intersect_random.csv.gz",
       git_dir + "/data/analysis/MAGs/intersect_descriptors.csv.gz"],
      [git_dir + "/data/analysis/genomes/fingerprints_intersect_random.csv.gz",
       git_dir + "/data/analysis/genomes/intersect_descriptors.csv.gz"]]
# analyze complete genomes only
io = [io[1]]
 
for io_tuple in io: 
    input_file = io_tuple[0]
    output_file = io_tuple[1] 
    print("processing input file: " + os.path.basename(input_file))
        
    # read input 
    dat = pd.read_csv(input_file)
    
    # create results container
    res = pd.DataFrame()
    
    # process antiSMASH SMILES
    for i, smiles in enumerate(dat['smiles_as']):
        mol = clean_mol(smiles)
        # calculate metrics
        metrics = calculate_metrics(mol)
        # create data frame
        row = pd.DataFrame({'cluster_idx': i+1,
                            'method': 'antiSMASH 5',
                            'smiles': smiles,
                            'metric': list(metrics.keys()), 
                            'value': list(metrics.values())})
        # append to results 
        res = res.append(row)
    
    # process PRISM 4 SMILES
    for i, smiles in enumerate(dat['smiles_pr']):
        mol = clean_mol(smiles)
        # calculate metrics
        metrics = calculate_metrics(mol)
        # create data frame
        row = pd.DataFrame({'cluster_idx': i+1,
                            'method': 'PRISM 4',
                            'smiles': smiles,
                            'metric': list(metrics.keys()), 
                            'value': list(metrics.values())})
        # append to results 
        res = res.append(row)
    
    # write     
    output_dir = os.path.dirname(output_file)
    if not os.path.isdir(output_dir): 
        os.mkdir(output_dir)
    
    res.to_csv(output_file, index=False, compression='gzip')
