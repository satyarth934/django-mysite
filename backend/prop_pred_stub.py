#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Test script to stub out property prediction'''
'''
Usage: prop_pred_stub.py <\"smiles list\"> <\"property dict\">")

Arguments:
    smiles list - double-quoted string with python list of smiles mols
        Example: 
"['C=CC[C@@H](CCCC)CCC(=O)O','C=CC[C@H](CCC(=O)O)CC(C)C','C=CC[C@@H](CC)C[C@H](CC)C(=O)O']"
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
    Example:
{'CN_dev': [0.5194887091410073, 0.6983061179806591, 0.44325813269693826], 'RON_dev': [0.4430912933346375, 0.7839214386250254, 0.6446953053267908], 'CN': [15.52399098508925, 11.442535904463284, 10.676358281348532], 'RON': [10.127942411390203, 16.890252657766883, 19.81320505963582]}
'''

# Imports
from datetime import datetime
import random
import sys

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
    raise Exception("Usage: prop_pred_stub.py <\"smiles list\"> <\"property dict\">")

if debug: 
    for i, arg in enumerate(sys.argv):
        print("argv["+str(i)+"] : "+sys.argv[i])

# argv[1] is the smiles list
if debug: 
    print("\nSmiles list: ")
smiles = eval(sys.argv[1])
if debug:
    print(smiles)
if type(smiles)!=list:
    raise Exception("Smiles list is not formatted correctly for 'eval'!")
if len(smiles)==0:
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

for p in prop:
    if debug:
        print("Property: " + p)
    for c in smiles:
        if debug:
            print("  Smiles: " + c)
        val = 10*(random.random()+1)
        retval[p].append(val)
        if debug:
            print("  value = "+"{:.2}".format(val))
        dev = random.random()
        if debug:
            print("  std deviation = "+"{:.2}".format(dev))
        retval[p+'_dev'].append(dev)

# This is the results section, surrounded it with tokens for parsing
print("\n" + begin_token)
print(retval)
print(end_token+"\n")

# Output end time
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("End Time = "+current_time)

