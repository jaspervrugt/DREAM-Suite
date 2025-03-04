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
## Example 12: Model-based geostatistics                                              ##
##                                                                                    ##
## Check the following paper                                                          ##
##  Minasny, B., J.A. Vrugt, and A.B. McBratney (2011), Confronting uncertainty in    ##
##      model-based geostatistics using Markov chain Monte Carlo simulation,          ##
##      Geoderma, 163, 150-162, doi:10.1016/j.geoderma.2011.03.011                    ##
##                                                                                    ##
## ---------------------------------------------------------------------------------- ##

import numpy as np
import os, sys
from scipy.stats import multivariate_normal

current_dir = os.getcwd()                                       # Get the current working directory
from blpmodel import postproc_variogram
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from DREAM_Suite_functions import genparset	                    # Import functions


# Define parameters
DREAMPar = {}
DREAMPar['d'] = 5       # Dimension of the problem
DREAMPar['lik'] = 2     # Model output is log-likelihood

# Provide information parameter space and initial sampling
Par_info = {}
Par_info['initial'] = 'latin'                   # Latin hypercube sampling
Par_info['boundhandling'] = 'reflect'           # Explicit boundary handling
Par_info['min'] = [0.0, 0.0, 0.00, 0.00, 0.00]  # If 'latin', min values
Par_info['max'] = [100, 100, 1000, 1000, 20.0]  # If 'latin', max values

# Define name of function (.py file) for posterior exploration
Func_name = 'blpmodel.blpmodel'

# Define method to use {'dream', 'dream_zs', 'dream_d', 'dream_dzs', 'mtdream_zs'}
method = 'dream_zs'

# Switch statement equivalent
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10          # Number of Markov chains
    DREAMPar['T'] = 2000        # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3           # Number of Markov chains
    DREAMPar['T'] = 2000        # Number of generations

# Handle the case where the method is 'dream_d' or 'dream_dzs'
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = [1000] * DREAMPar['d']  # Discrete steps

# Call the DREAM-Suite package
if __name__ == '__main__':
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info)

    P = genparset(chain)                           # Unpack the chain trajectories
    postproc_variogram(DREAMPar,P)                 # Postprocessor