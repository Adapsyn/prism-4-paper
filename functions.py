import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
from rdkit.DataStructs import FingerprintSimilarity, ConvertToNumpyArray

def clean_mol(smiles):
    """
    Construct a molecule from a SMILES string, removing stereochemistry and
    explicit hydrogens, and setting aromaticity. 
    """
    mol = Chem.MolFromSmiles(str(smiles), sanitize=64)
    if mol is None:
        raise ValueError("Invalid SMILES")
    Chem.RemoveStereochemistry(mol)
    Chem.SanitizeMol(mol)
    mol = Chem.RemoveHs(mol)
    return mol 

def clean_mols(all_smiles):
    """
    Construct a list of molecules from a list of SMILES strings, replacing
    invalid molecules with None in the list. 
    """
    mols = []
    for smiles in all_smiles: 
        try:
            mol = clean_mol(smiles)
            mols.append(mol)
        except ValueError:
            mols.append(None)
    return mols

def in_Ro5(mol):
    """
    Test whether a molecule is in Lipinski "Rule of 5" space, meaning 
    - 5 or fewer H bond donors
    - 10 or fewer H bond acceptors
    - MW < 500 Da
    - logP < 5
    """
    
    h_donor = Lipinski.NumHDonors(mol)
    h_accept = Lipinski.NumHAcceptors(mol)
    mw = Descriptors.MolWt(mol)
    logP = Descriptors.MolLogP(mol)
    
    Ro5 = h_donor <= 5 and h_accept <= 10 and mw <= 500 and logP < 5
    return(Ro5)

def get_ecfp6_fingerprints(mols): 
    """
    Get ECFP6 fingerprints for a list of molecules which may include `None`s,
    gracefully handling `None` values by returning a `None` value in that 
    position. 
    """
    fps = []
    for mol in mols:
        if mol is None:
            fps.append(None)
        else:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024)
            fps.append(fp)
    return(fps)
    
def get_bit_vector(fp):
    arr = np.zeros((1,))
    ConvertToNumpyArray(fp, arr)
    return(arr)

def get_tanimoto(list1, list2):
    tcs = []
    for fp1 in list1:
        for fp2 in list2: 
            if fp1 is None or fp2 is None:
                tcs.append(None)
            else:
                tc = FingerprintSimilarity(fp1, fp2)
                tcs.append(tc)
    return(tcs)

