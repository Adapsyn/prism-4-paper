"""
Natural product-likeness and synthetic accessibility scores, as provided in the
rdkit 'contrib' module.
""" 

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.six.moves import cPickle
from rdkit.six import iteritems
import math
import gzip
import pickle
from collections import namedtuple


"""
Natural product-likeness
""" 

def readNPModel(filepath='publicnp.model.gz'):
  """Read the NP scoring model."""
  fscore = pickle.load(gzip.open(filepath))
  return fscore

def scoreMolWConfidence(mol, fscore):
  """Next to the NP Likeness Score, this function outputs a confidence value
  between 0..1 that descibes how many fragments of the tested molecule
  were found in the model data set (1: all fragments were found).
  Returns namedtuple NPLikeness(nplikeness, confidence)"""

  if mol is None:
    raise ValueError('invalid molecule')
  fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
  bits = fp.GetNonzeroElements()

  # calculating the score
  score = 0.0
  bits_found = 0
  for bit in bits:
    if bit in fscore:
      bits_found += 1
      score += fscore[bit]

  score /= float(mol.GetNumAtoms())
  confidence = float(bits_found / len(bits))

  # preventing score explosion for exotic molecules
  if score > 4:
    score = 4. + math.log10(score - 4. + 1.)
  elif score < -4:
    score = -4. - math.log10(-4. - score + 1.)
  NPLikeness = namedtuple("NPLikeness", "nplikeness,confidence")
  return NPLikeness(score, confidence)

def calculateNPScore(mol, fscore):
  """Calculates the Natural Product Likeness of a molecule.
  Returns the score as float in the range -5..5."""
  return scoreMolWConfidence(mol, fscore).nplikeness

"""
Synthetic accessibility
"""

def readSAModel(filepath='fpscores.pkl.gz'):
  fscores = cPickle.load(gzip.open(filepath))
  outDict = {}
  for i in fscores:
    for j in range(1, len(i)):
      outDict[i[j]] = float(i[0])
  return outDict

def numBridgeheadsAndSpiro(mol, ri=None):
  nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
  nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
  return nBridgehead, nSpiro

def calculateSAScore(m, fscores):
  # fragment score
  fp = rdMolDescriptors.GetMorganFingerprint(m, 2)
  fps = fp.GetNonzeroElements()
  score1 = 0.
  nf = 0
  for bitId, v in iteritems(fps):
    nf += v
    sfp = bitId
    score1 += fscores.get(sfp, -4) * v
  score1 /= nf

  # features score
  nAtoms = m.GetNumAtoms()
  nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
  ri = m.GetRingInfo()
  nBridgeheads, nSpiro = numBridgeheadsAndSpiro(m, ri)
  nMacrocycles = 0
  for x in ri.AtomRings():
    if len(x) > 8:
      nMacrocycles += 1

  sizePenalty = nAtoms**1.005 - nAtoms
  stereoPenalty = math.log10(nChiralCenters + 1)
  spiroPenalty = math.log10(nSpiro + 1)
  bridgePenalty = math.log10(nBridgeheads + 1)
  macrocyclePenalty = 0.
  # ---------------------------------------
  # This differs from the paper, which defines:
  #  macrocyclePenalty = math.log10(nMacrocycles+1)
  # This form generates better results when 2 or more macrocycles are present
  if nMacrocycles > 0:
    macrocyclePenalty = math.log10(2)

  score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - \
      bridgePenalty - macrocyclePenalty

  # correction for the fingerprint density
  # not in the original publication, added in version 1.1
  # to make highly symmetrical molecules easier to synthetise
  score3 = 0.
  if nAtoms > len(fps):
    score3 = math.log(float(nAtoms) / len(fps)) * .5

  sascore = score1 + score2 + score3

  # need to transform "raw" value into scale between 1 and 10
  min = -4.0
  max = 2.5
  sascore = 11. - (sascore - min + 1) / (max - min) * 9.
  # smooth the 10-end
  if sascore > 8.:
    sascore = 8. + math.log(sascore + 1. - 9.)
  if sascore > 10.:
    sascore = 10.0
  elif sascore < 1.:
    sascore = 1.0

  return sascore
