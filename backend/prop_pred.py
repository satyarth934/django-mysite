#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Script to to call Graphdot property prediction'''
'''
Usage: prop_pred.py <\"smiles list\"> <\"property dict\">")

Arguments:
    smiles list - double-quoted string with python list of smiles mols
        Example: the first 5 are diesel fuel like molecules with known cetane 
        value, and the final one is a natural PKS product.
            CCCCCCCCCCCCCCCC Cetane (cetane=100)
            CC1=CC=CC2=CC=CC=C12 METHYLNAPHTHALENE (cetane=0)
            CC(CC(C)(C)C)CC(C)(C)CC(C)(C)C Isocetane (cetane=15)
            CC(C)CC(C)(C)C isooctane (octane=100)
            O=C(O)CCCCCCCCCCC\C=C/CCCCCCCC Erucic_acid (cetane=76)
            CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O)C)C)C Narbonolide
        Input string:
"['CCCCCCCCCCCCCCCC','CC1=CC=CC2=CC=CC=C12','CC(CC(C)(C)C)CC(C)(C)CC(C)(C)C','CC(C)CC(C)(C)C','O=C(O)CCCCCCCCCCC\\C=C/CCCCCCCC', 'CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O)C)C)C']"
        NOTE: double-quote is necessary to escape UNIX special characters
        and guarantee it will be treated as a single argument.

    property dict - double-quoted string with python dict of properties
        Allowed values (must have at least 1) for GraphDot trained models are:
          CN - cetane number 
          FP - flash point
          H1 - receptor pKd
          M2 - receptor pKd
          MP - melting point
          RON - research octane number 
          YSI - yield sooting index
        Example:
"{'RON':None, 'CN':None}" 

Return value:
    updated property dict with lists for each corresponding smiles input:
        - For each property key, a list of (float) predicted properties
        - Each property key + "_dev", a list of (float) property std deviations
    Example result:
{'CN': [[97.19695487945678, 3.8973963406438017, 13.056091798458741, 17.639384434056694, 75.31386408981265, 14.385125373410148]], 'CN_dev': [[3.154927075147002, 3.307434092849797, 3.63323350711247, 3.8969748875584402, 3.1467412092260005, 16.129204600920364]], 'RON': [[-21.449288379830932, 131.35693140662647, 79.00334388027292, 100.71953609770654, -0.6946304184010046, 95.24283286199183]], 'RON_dev': [[8.504497884581749, 8.17245571030485, 4.812979820250313, 3.2754367790940386, 4.238229817268251, 9.561777267423928]]}
'''

# Imports
from datetime import datetime
import random
import sys

import numpy as np
import networkx as nx

import re
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.metrics import r2_score

from graphdot import Graph
from graphdot.metric.maximin import MaxiMin
import graphdot.microkernel as uX
from graphdot.kernel.marginalized import MarginalizedGraphKernel
from graphdot.kernel.fix import Normalization
from graphdot.kernel.marginalized.starting_probability import Uniform
from graphdot.model.gaussian_process import GaussianProcessRegressor
from graphdot.model.gaussian_field.gfr import GaussianFieldRegressor
from graphdot.model.gaussian_field.weight import RBFOverDistance

from rdkit import Chem
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import AllChem

import pandas as pd


# Constants
#   Property codes - TODO, should probably be in a shared file
properties = ['CN', 'FP', 'H1', 'M2', 'MP', 'RON', 'YSI']
#   Delimeter tokens to produce results that can be parsed
begin_token = "---Begin Property Predictions---"
end_token = "---End Property Predictions---"
#   True to enable informational debug output
# debug = True
debug = False

# Output start time
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Start Time = "+current_time)

# Convert command line args to smiles list and property dict
nargs = len(sys.argv)
if nargs!=3:
    raise Exception("Usage: prop_pred.py <\"smiles list\"> <\"property dict\">")

if debug: 
    for i, arg in enumerate(sys.argv):
        print("argv["+str(i)+"] : "+sys.argv[i])

# argv[1] is the smiles list
if debug: 
    print("\nSmiles list: ")
smiles_query = eval(sys.argv[1])
if debug:
    print(smiles_query)
if type(smiles_query)!=list:
    raise Exception("Smiles list is not formatted correctly for 'eval'!")
if len(smiles_query)==0:
    raise Exception("Smiles list is empty!")

# argv[2] is the property dict
if debug:
    print("\nProp dict: ")
prop = eval(sys.argv[2])
if debug:
    print(prop)
if type(prop)!=dict:
    raise Exception("Propertry dict is not formatted correctly for 'eval'!")
if len(prop)==0:
    raise Exception("Property dict is empty!")

# copy return values into a new dictionary
retval = dict()
for p in prop:
    if not (p in properties):
        raise Exception("Property dict contains unknown code: " + p)
    retval[p] = list()
    retval[p+'_dev'] = list()


'''Graphdot modeling of biofuel and drug molecular properties'''
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Done with imports, Time =", current_time)

def make_kernel():
    return Normalization(MarginalizedGraphKernel(
        node_kernel=uX.Additive(
            aromatic=uX.Constant(0.5, c_bounds=(0.01, 1.0)) * uX.KroneckerDelta(0.5, (0.1, 0.9)),
            atomic_number=uX.Constant(0.5, c_bounds=(0.01, 1.0)) * uX.KroneckerDelta(0.5, (0.1, 0.9)),
            charge=uX.Constant(0.5, c_bounds=(0.01, 1.0)) * uX.SquareExponential(1.0),
            chiral=uX.Constant(0.5, c_bounds=(0.01, 1.0)) * uX.KroneckerDelta(0.5, (0.1, 0.9)),
            hcount=uX.Constant(0.5, c_bounds=(0.01, 1.0)) * uX.SquareExponential(1.0),
            hybridization=uX.Constant(0.5, c_bounds=(0.01, 1.0)) * uX.KroneckerDelta(0.5, (0.1, 0.9)),
            ring_list=uX.Constant(0.5, c_bounds=(0.001, 1.0)) * uX.Convolution(uX.KroneckerDelta(0.5, (0.1, 0.9))).normalized
        ).normalized,
        edge_kernel=uX.Additive(
            aromatic=uX.Constant(0.5, (0.1, 1.0)) * uX.KroneckerDelta(0.5, (0.1, 0.9)),
            conjugated=uX.Constant(0.5, (0.1, 1.0)) * uX.KroneckerDelta(0.5, (0.1, 0.9)),
            order=uX.Constant(0.5, (0.1, 1.0)) * uX.KroneckerDelta(0.5, (0.1, 0.9)),
            ring_stereo=uX.Constant(0.5, (0.1, 1.0)) * uX.KroneckerDelta(0.5, (0.1, 0.9)),
            stereo=uX.Constant(0.5, (0.1, 1.0)) * uX.KroneckerDelta(0.5, (0.1, 0.9))
        ).normalized,
        p=Uniform(1.0, p_bounds='fixed'),
        q=0.05
    ))

prop_gpr = dict()
path = '/global/cfs/cdirs/m3513/molinv/'
if debug:
    print("Loading models from: "+path+"models")

for p in prop:
    print('Creating ' + p + ' predictor ...')
    gpr = GaussianProcessRegressor(
        kernel=make_kernel(), alpha=1e-2, optimizer=True, normalize_y=True)
    gpr.load(path+'models', p+'_Graphdot_GPR.pickle')
    prop_gpr[p] = gpr

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Done with GPR model initialization, Time =", current_time)

if debug:
    print(smiles_query)

molout = []

for smilesin in smiles_query:
    if debug:
        print(smilesin)
    inputmol = MolFromSmiles(smilesin)
    Chem.SanitizeMol(inputmol)
    molout.append(inputmol)

# graphs_query = [Graph.from_rdkit(productmol)]
# graphs_query = [Graph.from_rdkit(mol) for mol in mols_query] 
print('Starting property prediction')
graphs_query = [Graph.from_rdkit(mol) for mol in molout]

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Done with mol prep, Time =", current_time)

# print(graphs_query)
for p in prop:
    if debug:
        print("Property: "+p+' predictions:')
    z, s = prop_gpr[p].predict(graphs_query, return_std=True)
    if debug:
        print(z.tolist())
        print(s.tolist())
    retval[p] = z.tolist()
    retval[p+'_dev'] = s.tolist()


now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Prediction done, Time =", current_time)

# This is the results section, surrounded it with tokens for parsing
print("\n" + begin_token)
print(retval)
print(end_token+"\n")

# Output end time
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("End Time = "+current_time)

