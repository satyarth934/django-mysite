#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Testing superfacility API call for "bioarc" collaboration account
calling the GraphDot property predictor, or stubs for testing

For more information on the SF API:
    https://docs.nersc.gov/services/sfapi/
    https://api.nersc.gov/api/v1.2/#/
'''

from authlib.integrations.requests_client import OAuth2Session
from authlib.oauth2.rfc7523 import PrivateKeyJWT
import os 
import json 
from pathlib import Path 
import requests
import time 
from datetime import datetime

class PropertyPredictor:

    '''
    Constructor to setup the sfapi client and system/job to run:
      path - NOTE: this is where the job will run and output goes
      client_id - string, SFAPI client id (or env var name with)
      client_pem - string, SFAPI client private PEM key (or env var name with)
      client_env - if True, get the sfapi client info from environment vars
      target_job - 'cori_test', 'perl_test', (stub tests) or 'perl' (real)
      debug - 0 (minimal), 1 (some), 2 (verbose) debug output
    '''
    def __init__(self, path, client_id, client_pem, client_env=True,
                 target_job="perl", debug=0):
        self.token_url = "https://oidc.nersc.gov/c2id/token"
        self.sfapi_url = "https://api.nersc.gov/api/v1.2"
        self.poll_sec = 2
        self.debug = debug
        # (ignore) client_sys - 'cori', 'perl', or 'spin' (where client is running from)
        # Setup for where the client is calling from
        # if client_sys == "cori": # Cori login nodes
        #     self.client_id = "65cdtjh2bdw6m"
        #     self.private_key = Path(self.path+"/private_cori.pem").read_text()
        # elif client_sys == "perlmutter": # Perlmutter login nodes
        #     self.client_id = "dzgc2ttfsaduw"
        #     self.private_key = Path(self.path+"/private_cori.pem").read_text()
        # elif client_sys == "spin": # Spin instance nodes
        #     self.client_id = "fljm3ujlzjnby"
        #     self.private_key = Path(self.path+"/private_spin.pem").read_text()
        # else: # Not supported
        #     raise Exception("Client system = " + client_sys + " not supported.")
        # self.client_sys = client_sys

        # Setup for what code is running where
        if target_job == "perl": # Run real jobs on Perlmutter
            self.target_sys = "perlmutter"
            self.test = False
        elif target_job == "perl_test": # Run test script on Perlmutter
            self.target_sys = "perlmutter"
            self.test = True
        elif target_job == "cori_test": # Run test script in Cori Haswell
            self.target_sys = "cori"
            self.test = True
        else: # Not supported
            raise Exception("Target system = " + target_sys + " not supported.")

        if len(path)==0:
            raise Exception("Empty path to job project? Required arg!")
        if debug > 0:
            print("\n Using job path: " + path) 
        self.path = path

        # Check that the environment variables exists and has non-zero len
        if client_env:
            client_id_env = client_id # copy env var names
            client_id = os.getenv(client_id_env)
        if (client_id is None) or (len(client_id)==0):
            raise Exception("Empty SFAPI client id empty: " + \
                    client_id_env)
        self.client_id = client_id
        print("Client id:\n"+self.client_id+"\n")
        if client_env:
            client_pem_env = client_pem
            client_pem = os.getenv(client_pem_env)
        if (client_pem is None) or (len(client_pem)==0):
            raise Exception("Empty SFAPI client PEM empty: " + \
                    client_pem_env)
        self.client_pem = client_pem
        print("Client PEM:\n"+self.client_pem+"\n")

    def __str__(self):
        vals = dict()
        vals['path'] = self.path
        vals['token_url'] = self.token_url
        vals['sfapi_url'] = self.sfapi_url
        vals['target_sys'] = self.target_sys
        vals['test'] = self.test
        vals['client_id'] = self.client_id
        vals['debug'] = self.debug
        return str(vals) 

    # def __del__(self):
    #     print("PropertyPredictor destroyed")
    #     # TODO: Might need to close session, etc.

    '''
    Open a session with the sfapi client_id and private key file
    '''
    def open_session(self):
        # TODO: check if a session is already open 
        self.session = OAuth2Session(self.client_id, self.client_pem,
            PrivateKeyJWT(self.token_url), grant_type="client_credentials",
            token_endpoint=self.token_url)
        self.session.fetch_token()
        if self.debug > 0:
            print('Got SFAPI access token!')
        if self.debug > 1:
            print(self.session.token['access_token'])

    '''Check the target system status
    Return value: string - 'unavailable' or 'available'
    '''
    def check_status(self):
        if self.debug > 0:
            print('\nSystem status API call:')
        url = self.sfapi_url+"/status/"+self.target_sys
        retval = self.session.get(url)
        if self.debug > 0:
            print(retval.json())
        status = retval.json()['status']
        return status

    def submit_query(self, smiles, props):
        if self.debug:
            print('\nSubmit job API call:')
        mols_string = "\"" + str(smiles) + "\""
        props_string = "\"" + str(props) + "\""
        if self.debug > 0:
            print("Molecule list :" + mols_string)
            print("Property dict :" + props_string)

        # run - different commands to run: test stubs or real thing
        # header - different batch parameters and module loads
        # command - different command parameters

        if self.test==True:
            run="/prop_pred_stub.py"
        else:
            run="/prop_pred.py"

        if self.target_sys=="perlmutter":
            header = "#!/bin/bash\n#SBATCH -A m3513_g\n#SBATCH -C gpu\n#SBATCH -q regular\n#SBATCH -t 0:10:00\n#SBATCH -n 1\n#SBATCH -o " + self.path + "/%j.log\n#SBATCH --ntasks-per-node=1\n#SBATCH -c 128\n#SBATCH --gpus-per-task=1\n\n\nmodule load python\nsource /global/homes/u/u6336/venv/graphdot/bin/activate\n\nexport SLURM_CPU_BIND=\"cores\"\n"
            command = "\nsrun -n 1 -c 64 --cpu_bind=cores "+self.path+run
        else: # self.target_sys=="cori", stub only
            header = "#!/bin/bash\n#SBATCH -A m3513\n#SBATCH -N 1\n#SBATCH -C haswell\n#SBATCH -q regular\n#SBATCH -t 00:05:00\n#SBATCH -o " + self.path + "/%j.log\n\n"
            command = "\nsrun -n 1 -c 64 --cpu_bind=cores "+self.path+run

        # Build up the submit script as a string
        submit_script = header+command+" "+mols_string+" "+props_string
        if self.debug > 1:
            print("Submit script")
            print(submit_script)
        url = self.sfapi_url+"/compute/jobs/"+self.target_sys
        retval = self.session.post(url,
                data = {"job": submit_script, "isPath": False})
        if self.debug > 1:
            print(retval.json())

        # Get the SFAPI task id, wait for it to return a job id
        task_id = retval.json()['task_id']
        if self.debug > 0:
            print('\nTask ' + task_id + ' status API call:')
        url = self.sfapi_url+"/tasks/"+task_id
        while True:
            r = self.session.get(url)
            status = r.json()['status']
            now = datetime.now()
            if self.debug > 0:
                current_time = now.strftime("%H:%M:%S")
                print('    ' + current_time + ': ' + status)
                print(r.json())
            if status != 'new':
                break
            else:
                time.sleep(self.poll_sec)

        # TODO - check that there's no error and status in result

        # Get the job id and return it
        resultstr = r.json()['result']
        if self.debug > 1:
            print(resultstr)
        resultstr = resultstr.replace('null','None')
        result = eval(resultstr)
        if result['error']!=None:
            job_id = "ERROR"
            if self.debug > 0:
                print('\nTask ' + task_id + ' failed with error!')

        job_id = result['jobid']
        return job_id

    '''Check the job_id status using sfapi for sacct'''
    def job_status(self, job_id):
        retval = dict()
        if self.debug > 0:
            print('\nSacct job ' + job_id + ' status API call:')
        url = self.sfapi_url+"/compute/jobs/"+self.target_sys +"/"+job_id+"?sacct=true"
        if self.debug > 0:
            print(url)
        result = self.session.get(url)
        if self.debug > 1:
            print(result.json())

        output = result.json()['output']
        error = result.json()['error']
        if self.debug > 0: 
            print("Error : " + str(error))
        if error!=None:
            status = "ERROR"
        else:
            status = output[0]['state']
        if self.debug > 0:
            print("Status of " + job_id + " : " + status)

        return status

    '''Check the job_id status using sfapi for squeue'''
    # NOTE: this will not find jobs that have complet
    def job_status_squeue(self, job_id):
        if self.debug > 0:
            print('\nSqueue job ' + job_id + ' status API call:')
        url = self.sfapi_url+"/compute/jobs/"+self.target_sys+"/"+job_id
        if self.debug > 0:
            print(url)
        result = self.session.get(url)
        if self.debug > 1:
            print(result.json())
        output = result.json()['output']
        if len(output) > 0:
            status = output[0]['state']
        else:
            status = "NONE"
        return status

    '''Parse the output file string for the property dictionary results'''
    def parse_results_file(self, contents):
        # Find the tags in the output lines, return the property dictionary
        start_token= "---Begin Property Predictions---"
        end_token = "---End Property Predictions---"
        start = contents.find(start_token)+len(start_token)
        end = contents.find(end_token)
        return eval(contents[start:end])

    '''Get the slurm job output file, return the property dictionary results'''
    def get_query_results(self, job_id):
        # filepath = 'global/homes/u/u6336/repos/bioarc/slurm-3280437.out'
        filepath = self.path+"/"+job_id+".log"
        if self.debug > 0:
            print('\nDownload file API call:')
            print(filepath)
        getcmd = self.sfapi_url+"/utilities/download/"+self.target_sys+"/"+filepath
        if self.debug > 0:
            print(getcmd)
        result = self.session.get(getcmd).json()
        if self.debug > 0:
            print(result)
        if result['status']=="ERROR":
            return "ERROR"
        contents = result['file']
        return self.parse_results_file(contents)

