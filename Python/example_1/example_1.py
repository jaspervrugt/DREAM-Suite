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
## Example 2: w1*N(µ1,Σ1) + w2*N(µ2,Σ2) with µ1 = -5, µ2 = 5, Σ1/Σ2 = identity matrix ##
##                   and w1 = 1/3 and w2 = 2/3                                        ##
##                                                                                    ##
## Check the following papers                                                         ##
##  Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M.    ##
##      Hyman (2009), Accelerating Markov chain Monte Carlo simulation by             ##
##      differential evolution with self-adaptive randomized subspace sampling,       ##
##      International Journal of Nonlinear Sciences and Numerical Simulation, 10(3),  ##
##      271-288.                                                                      ##
##   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian (2003), A Shuffled Complex  ##
##      Evolution Metropolis algorithm for optimization and uncertainty assessment of ##
##      hydrologic model parameters, Water Resour. Res., 39 (8), 1201,                ##
##      doi:10.1029/2002WR001642.                                                     ##
##                                                                                    ##
## ---------------------------------------------------------------------------------- ##

import numpy as np
import os, sys

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path

from DREAM_Suite import DREAM_Suite

## Problem settings defined by user
DREAMPar = {
    'd': 2,     # Dimension of the problem
    'lik': 2,   # Model output is log-likelihood
}

# Parameter info for the normal distribution (mean and covariance matrix)
Par_info = {
    'initial': 'Normal',                # Initial distribution is Normal
    'mu': np.zeros(DREAMPar['d']),      # Mean of the distribution (mu = 0)
    'cov': 10 * np.eye(DREAMPar['d']),  # Covariance matrix (10 * identity matrix)
}

# Define the function for posterior exploration (name of the function for likelihood)
Func_name = 'banana_lik.banana_lik'

# Method for DREAM algorithm
method = 'dream_zs'

# Switch-like structure for setting Markov chains and generations based on the method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10                  # Number of Markov chains
    DREAMPar['T'] = 100000              # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3                   # Number of Markov chains
    DREAMPar['T'] = 20000               # Number of generations

# Handle 'dream_d' or 'dream_dzs' methods for discrete sampling
if method in ['dream_d', 'dream_dzs']:
    Par_info['min'] = -100 * np.ones(DREAMPar['d'])     # Min value for discrete sampling
    Par_info['max'] = 100 * np.ones(DREAMPar['d'])      # Max value for discrete sampling
    Par_info['steps'] = 2001 * np.ones(DREAMPar['d'])   # Discrete steps

if __name__ == '__main__':
    # Call the DREAM-Suite function
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info)