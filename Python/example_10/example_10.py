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
## Example 10: w1*N(µ1,Σ1) + w2*N(µ2,Σ2) with µ1 = -5, µ2 = 5, Σ1/Σ2 =identity matrix ##
##                   and w1 = 1/3 and w2 = 2/3. Informative prior                     ##
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
from scipy.stats import multivariate_normal

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite


# Problem settings defined by user
DREAMPar = {'d': 2,  	# Dimension of the problem
            'lik': 2}  	# Model output is log-likelihood

# Provide information parameter space and initial sampling
Par_info = {'initial': 'normal',  		    # Multinormal initial sampling distribution
    'mu': np.zeros(DREAMPar['d']), 	        # If 'normal', define mean of distribution
    'cov': np.eye(DREAMPar['d']),  	        # If 'normal', define covariance matrix
    'a': np.array([-2, -2]),  		        # Mean of prior distribution
    'b': np.eye(2)}  			            # Covariance of prior distribution
# Bi-variate prior distribution
Par_info['prior'] = [multivariate_normal(mean = Par_info["a"], cov = Par_info["b"])]

# Define name of function (.py file) for posterior exploration
Func_name = 'mixture_lik.mixture_lik'

# Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs'

if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  	# Markov chains
    DREAMPar['T'] = 2000  	# generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3  		# Markov chains
    DREAMPar['T'] = 7000  	# generations

if method in ['dream_d', 'dream_dzs']:
    Par_info['min'] = np.array([-20, -20])  		    # Min value for discrete sampling
    Par_info['max'] = np.array([20, 20])  		        # Max value for discrete sampling
    Par_info['steps'] = 1000 * np.ones(DREAMPar['d'])   # discrete steps

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info)
