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
## Example 18: Lotka-Volterra model: Informal likelihood (GLUE)                       ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. and K.J. Beven (2018), Embracing equifinality with efficiency:       ##
##       Limits of Acceptability sampling using the DREAM_{(LOA)} algorithm,  Journal ##
##       of Hydrology,  559 , pp. 954-971, doi:10.1016/j.jhydrol.2018.02.026          ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys
import cProfile 

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

## Problem settings defined by user
DREAMPar = {
    'd': 4,              # Dimension of the problem
    'lik': 32,           # Informal likelihood function
    'GLUE': 10}          # Value (default) of likelihood shape parameter

# Define initial sampling and parameter ranges
Par_info = {
    'initial': 'latin',      # Latin hypercube sampling
    'boundhandling': 'reflect',  # Boundary handling: reflection
    'names': ['\\alpha', '\\beta', '\\gamma', '\\delta'],  # Parameter names
    'min': [0, 0, 0, 0],     # Lower bound parameters
    'max': [1, 10, 1, 10]}   # Upper bound parameters

# Load calibration data vector
Meas_info = {'Y': np.loadtxt('abundances.txt')}  # Load the food web dataset
Meas_info['Sigma'] = 1

# Define name of function for posterior exploration
Func_name = 'Lotka_Volterra.Lotka_Volterra'

# Optional settings
options = {
    'parallel': 'yes',	# Run chains in parallel?
    'modout': 'yes',    # Model returns model output [= 1st or 2nd output argument]
    'IO': 'no'}         # Master-worker communication using file writing?

# Define method to use {'dream', 'dream_zs', 'dream_d', 'dream_dzs', 'mtdream_zs'}
method = 'dream_zs'

# Set Markov chains and generations based on the method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10          # Number of Markov chains
    DREAMPar['T'] = 2500        # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3           # Number of Markov chains
    DREAMPar['T'] = 2000        # Number of generations
elif method in ['dream_kzs']:
    DREAMPar['N'] = 3           # Number of Markov chains
    DREAMPar['T'] = 2000        # Number of generations
    DREAMPar['M'] = 24          # Number of samples for Kalman jump    
    DREAMPar['a_1'] = 0         # To activate Kalman jump

# If method is dream_d or dream_dzs, set the number of discrete steps
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = 1000 * np.ones(DREAMPar['d'])  # Number of discrete steps

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options)
    # cProfile.run('chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options)', 'profile_report.prof')
