#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Call superfacility API call for GraphDot property prediction'''

# Imports
import prop_sfapi as prop
import time
from datetime import datetime

# Output start time
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Starting test: " + current_time)

# Create the PropertyPredictor for client running on cori
# debug = 1 # produces more output
debug = 2 # produces a lot of output
# debug = 0 # produces minimal output

# 1. Choose one of these for where the client is running (for correct keys)
# client_sys = "spin"
# client_sys = "cori"
client_sys = "perl"
# 2. Choose pairs of these to indicate what/where to run
target_job = "perl" # call real code on Perlmutter GPUs
system = "perlmutter"
# target_job = "cori_test" # call test code on Cori Haswell nodes
# system = "cori"
# target_job = "perl_test" # call test code on Perlmutter AMD CPU (no GPU)
# system = "perlmutter"

# pp = prop.PropertyPredictor(client_sys, target_sys, debug)

path = "/global/cfs/cdirs/m3513/molinv/prod" # path to where batch job runs
client_id = "BIOARC_PERL_ID" # this environment variable = SFAPI client id
client_pem = "BIOARC_PERL_PEM" # this environment variable = private PEM key
client_env = True
pp = prop.PropertyPredictor(path, client_id, client_pem, client_env, \
        target_job, debug)

if debug > 0:
    print(pp)

# Open an SF API session - stays open until this is GC'ed
# TODO - check return code?
pp.open_session()

# Check the status of the target machine
# TODO - check return code?
status = pp.check_status()
print("System: "+system+", status: "+status)
if status=='unavailable':
    print("Target system: "+system+" status='unavailable', exiting")
    quit()

# Prepare the test input
#   list of smiles strings, double-quoted for command line
# mols_string = "\"['C=CC[C@@H](CCCC)CCC(=O)O','C=CC[C@H](CCC(=O)O)CC(C)C','C=CC[C@@H](CC)C[C@H](CC)C(=O)O']\""
mols_string = "\"['CCCCCCCCCCCCCCCC', 'CC1=CC=CC2=CC=CC=C12', 'CC(CC(C)(C)C)CC(C)(C)CC(C)(C)C', 'CC(C)CC(C)(C)C', 'O=C(O)CCCCCCCCCCC\\C=C/CCCCCCCC', 'CC[C@@H]1[C@@H](/C=C/C(=O)[C@@H](C[C@@H]([C@@H]([C@H](C(=O)[C@H](C(=O)O1)C)C)O)C)C)C']\""
smiles = eval(mols_string)
#   dict of property code strings, double-quoted for command line
props_string = "\"{'RON':None, 'CN':None}\"" 
props = eval(props_string)

# Submit the input and wait (several seconds) for a job id
job_id = pp.submit_query(smiles, props)

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("  Time " + current_time)
print("Job ID: " + job_id)

# Poll on the job status - simulate a user clicking "refresh" every 2 sec
print("\nPolling every 15 sec until job completes ...")
while True:
    status = pp.job_status(job_id)
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("  Time " + current_time + ": " + status)
    if debug > 0:
        squeue_status = pp.job_status_squeue(job_id)
        print("    Job status: squeue "+squeue_status+", sacct "+status)
    if status=='COMPLETED':
        break
    elif status=='FAILED':
        break
    else:
        time.sleep(15)

print("\nJob ID: " + job_id + " status : " + status)

# Process the job output file with the query results
if status=='COMPLETED' or status=='FAILED':
    result = pp.get_query_results(job_id)
    if debug > 0:
        print("Job ID: " + job_id + " output:")
    print(result)

# Output the end time
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Ending test: " + current_time)

