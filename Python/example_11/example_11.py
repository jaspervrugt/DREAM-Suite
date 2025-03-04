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
## Example 11: Target is a d-variate t-distribution with df degrees of freedom        ##
##             and correlation matrix R                                               ##
##                                                                                    ##
## Check the following papers                                                         ##
##  Ter Braak, C.J.F., and J.A. Vrugt (2008), Differential Evolution Markov Chain     ##
##      with snooker updater and fewer chains, Statistics and Computing,              ##
##      10.1007/s11222-008-9104-9.                                                    ##
##                                                                                    ##
## ---------------------------------------------------------------------------------- ##

import numpy as np
import os, sys
# from scipy.stats import multivariate_normal

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite                             # Import DREAM_Suite

# Problem specific parameter settings
DREAMPar = {'d': 25,  # Dimension of the problem
            'lik': 2} # Model output is log-likelihood

# Provide information parameter space and initial sampling
Par_info = {'initial': 'latin',                     # Latin hypercube sampling
            'min': -5 * np.ones(DREAMPar['d']),     # Lower bound parameter values
            'max': 15 * np.ones(DREAMPar['d']),     # Upper bound parameter values
            'norm': 1}                              # Sample in normalized space [0-1]
 
# Define name of function (.m file) for posterior exploration
Func_name = 'mvt_lik.mvt_lik'

# Define method to use
method = 'dream_zs'

# Set parameters for method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 20  	# Number of Markov chains
    DREAMPar['T'] = 5000  	# Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3  		# Number of Markov chains
    DREAMPar['T'] = 30000  	# Number of generations

# Discrete steps for 'dream_d' or 'dream_dzs'
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = 1000 * np.ones(DREAMPar['d'])

# Call the DREAM-Suite package
if __name__ == '__main__':
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info)
