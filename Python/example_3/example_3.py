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
## Example 3: w1*N(µ1,Σ1) + w2*N(µ2,Σ2) with µ1 = -5, µ2 = 5, Σ1/Σ2 = identity matrix ##
##                   and w1 = 1/3 and w2 = 2/3                                        ##
##                                                                                    ##
## Check the following papers                                                         ##
##  Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M.    ##
##      Hyman (2009), Accelerating Markov chain Monte Carlo simulation by             ##
##      differential evolution with self-adaptive randomized subspace sampling,       ##
##      International Journal of Nonlinear Sciences and Numerical Simulation, 10(3),  ##
##      271-288.                                                                      ##
##  Ter Braak, C.J.F., and J.A. Vrugt (2008), Differential Evolution Markov Chain     ##
##      with snooker updater and fewer chains, Statistics and Computing,              ##
##      10.1007/s11222-008-9104-9.                                                    ##
##                                                                                    ##
## ---------------------------------------------------------------------------------- ##

import numpy as np
import os, sys
import matplotlib.pyplot as plt

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from DREAM_Suite_functions import genparset		                # Import functions

sys.path.append(os.path.join(parent_dir, 'gamesampling'))	    # Add gamesampling directory to Python path
from GAME_sampling import *					                    # Import functions

# Define the parameters (DREAMPar)
DREAMPar = {
    'd': 10,               # Dimension of the problem
    'thinning': 2,         # Only store every 2nd sample (adjusted to Python-style)
    'lik': 2               # Model output is log-likelihood
}

# Define the function name for posterior exploration
Func_name = 'mixture_lik.mixture_lik'

# Provide information about the parameter space and initial sampling
Par_info = {
    'initial': 'normal',                    # Nd(µ, Σ) d-variate normal
    'mu': np.zeros(DREAMPar['d']),  		# If 'normal', µ-mean of distribution
    'cov': 5 * np.eye(DREAMPar['d']),  		# If 'normal', Σ-covariance matrix
    'min': -20 * np.ones(DREAMPar['d']),  	# Min value for discrete sampling
    'max': 20 * np.ones(DREAMPar['d']),   	# Max value for discrete sampling
    'norm': 0					            # Normalized parameter space
}

# Define method to use {'dream', 'dream_zs', 'dream_d', 'dream_dzs', 'mtdream_zs'}
method = 'dream_zs'

# Set values for DREAMPar based on method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = DREAMPar['d']        # Number of Markov chains
    DREAMPar['T'] = 50000                # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 10                   # Number of Markov chains
    DREAMPar['T'] = 25000                # Number of generations

# If method is dream_d or dream_dzs, set the number of discrete steps
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = 400 * np.ones(DREAMPar['d'])  # Number of discrete steps

if __name__ == '__main__':

    # Call the DREAM-Suite function
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info)

    # Generate parset data
    parset = genparset(chain)

    # Check the mean and variance of MCMC samples
    P = parset[-5000:, :]

    # Marginal likelihood, Z
    Z, logZ, gmix = GAME_sampling(P, 'is', DREAMPar, Func_name, Par_info)

    print("gmix", gmix, "Z", Z)

# Note: need to check GAME_sampling one more time: Z should equal 1