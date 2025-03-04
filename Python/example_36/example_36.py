## ################################################################################## ##
##                                                                                    ##
## DDDDD   RRRR    EEEEE    AA    MM   MM      SSSSSS  UU   UU   II   TTTTTTTT  EEEEE ##
## DDDDDD  RRRR    EEEEE   AAAA   MM   MM      SSSSS   UU   UU   II   TTTTTTTT  EEEEE ##
## DD  DD  RR RR   EE     AA  AA  MMM MMM      SS      UU   UU   II      TT     EE    ##
## DD  DD  RR RR   EEE    AA  AA  MMMMMMM ---- SS      UU   UU   II      TT     EEE   ##
## DD  DD  RRRRR   EEE    AAAAAA  MMM MMM ---- SSSSSS  UU   UU   II      TT     EEE   ##
## DD  DD  RR RR   EE     AAAAAA  MM   MM          SS  UU   UU   II      TT     EE    ##
## DDDDDD  RR  RR  EEEEE  AA  AA  MM   MM       SSSSS  UUUUUUU   II      TT     EEEEE ##
## DDDDD   RR  RR  EEEEE  AA  AA  MM   MM      SSSSSS  UUUUUUU   II      TT     EEEEE ##
##                                                                                    ##
## ################################################################################## ##
##                                                                                    ##
## Example 36: Inverse modeling of the parameters of a two class flocculation model   ##
##             using temperature data and limits of acceptability                     ##
##                                                                                    ##
## ---------------------------------------------------------------------------------- ##

import numpy as np
import os, sys
import pandas as pd

current_dir = os.getcwd()                                       # Get the current working directory
from flocsedtran_1DV import model_properties

parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

# Problem settings defined by user
DREAMPar = {'d': 2,                        # Dimension of the problem
            'lik': 23}                     # Likelihood 23 (limits acceptability)

# Provide information on parameter space and initial sampling
Par_info = {'initial': 'latin',                 # Latin hypercube sampling
            'boundhandling': 'reflect',         # Explicit boundary handling
            'names': ['M', '\\tau_{\\rm c}'],   # Parameter names
            'min': [0, 1.8],                    # Minimum parameter values
            'max': [1, 2.5]}                    # Maximum parameter values

# Define the name of the function for posterior exploration
Func_name = 'flocsedtran_1DV.flocsedtran_1DV'

# Load experimental data (assuming the file is available in the working directory)
data = pd.read_csv('MeasData_BLK2018.txt', header=None, delimiter='\t').to_numpy()

# Experimental data and times
n, K = data.shape                       # Number of rows and columns
ts_obs = 3600 * data[:, 0]              # Observation times (converted from hours to seconds)
data = data[:, 2:]                      # Isolate the relevant columns for the measurement data

# Define training data
Meas_info = {'S': data.flatten()}       # Flatten the data into a 1D array (temperature data)

# Setup numerical simulation model (assuming model_properties is implemented in Python)
plugin = model_properties()             # Call the model_properties function to get the plugin structure
plugin['ts_obs'] = ts_obs               # Add observation times to the plugin

# Define method to use
method = 'dream_zs'

# DREAM settings based on method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10                # Number of Markov chains
    DREAMPar['T'] = 2500              # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3                 # Number of Markov chains
    DREAMPar['T'] = 3000              # Number of generations

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (200) * DREAMPar['d']  # Discrete steps for parameters

# Optional settings
options = {'modout': 'yes',           # Return model simulations?
           'epsilon': [2],            # Limits of acceptability
           'parallel': 'no',          # Run each chain on a different core
           'save': 'yes',             # Save workspace during the run
           'print': 'yes'}            # Print figures to screen

# Run eDREAM package
if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, [], [], plugin)

    # Save results
    np.savez('results.npz', chain = chain, output = output, FX = FX, Z = Z, logL = logL, 
            method = method, Func_name = Func_name, DREAMPar = DREAMPar, 
            Par_info = Par_info, Meas_info = Meas_info, options = options, plugin = plugin)
