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
## Example 30: Predator prey interactions                                             ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##                                                                                    ##
## ################################################################################## ##

## DOWNLOAD AT:
## http://izt.ciens.ucv.ve/ecologia/Archivos/ECO_POB#202007/ECOPO2_2007/Cariboni#20et#20al#202007.pdf

import numpy as np
import os, sys
import scipy.io
import matplotlib.pyplot as plt

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from DREAM_Suite_functions import genparset                  # Import functions

# Problem settings defined by user
DREAMPar = {'d': 4,     # Number of parameters
            'lik': 11}  # Gaussian likelihood

# Provide information for parameter space and initial sampling
Par_info = {'initial': 'latin',  		    # Latin hypercube sampling
	 	    'boundhandling': 'reflect',  	# Explicit boundary handling
            'names': ['r','\\alpha','m','\\theta'],
            'min': [0.8, 0.2, 0.6, 0.05],  	# Minimum parameter values
            'max': [1.8, 1.0, 1.0, 0.15]}  	# Maximum parameter values

# Define name of the function for posterior exploration
Func_name = 'pred_prey'                     # [= allowed] 'pred_prey.pred_prey'

# Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'mtdream_zs'

# Set number of chains and generations based on the method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  		# Markov chains
    DREAMPar['T'] = 2000  		# generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3  			# Markov chains
    DREAMPar['T'] = 5000  		# generations
elif method == 'dream_kzs':
    DREAMPar['N'] = 5
    DREAMPar['T'] = 5000
    DREAMPar['M'] = 24  	    # Size of initial archive Kalman jump    

# Dimension of the problem
DREAMPar['d'] = len(Par_info['min'])

# Define the model setup
plugin = {  'u0': [8, 2],               # initial count of each of the two species (prey, predator)
            't': np.arange(0, 61, 1)}   # specify the time step and end time

# Load synthetic data of predator prey data
data = np.loadtxt('data.txt')

# Add noise to data
Y = np.array([data]).reshape(-1,1)
Y[:61, 0] = Y[:61, 0] + np.random.normal(0, 0.1 * np.std(Y[:61, 0]), 61)
Y[61:122, 0] = Y[61:122, 0] + np.random.normal(0, 0.1 * np.std(Y[61:122, 0]), 61)

# Store measured data
Meas_info = {'Y': Y}

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (250) * np.ones(DREAMPar['d'])  # Number of discrete steps

# Optional settings
options = {	'modout': 'yes',  	# Return model simulations (yes/no)?
    		'parallel': 'no', 	# Run each chain on a different core
    		'save': 'yes',  	# Save workspace during run
    		'print': 'yes'}		# No figures printed to screen

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, [], [], plugin)

    # Postprocess: Create one matrix of Markov chains
    P = genparset(chain)
    # Get the number of samples (rows) in pars
    M = P.shape[0]
    # Burn-in: keep only the last 80% of the samples
    P = P[int(0.8 * M):, :DREAMPar['d']]
    # Burn-in to simulated abundances
    Y = FX[int(0.8 * M):, :]
