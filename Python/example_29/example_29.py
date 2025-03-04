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
## Example 29: w1*N(µ1,Σ1) + w2*N(µ2,Σ2) with µ1 = -5, µ2 = 5, Σ1/Σ2 = identity mtrix ##
##                   and w1 = 1/3 and w2 = 2/3. Informative prior on parameters       ##
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
from scipy.stats import norm, gamma, uniform, genpareto, multivariate_normal

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

# Problem settings defined by user
DREAMPar = {'d': 5,  	# Dimension of the problem
            'lik': 2}  	# Model output is log-likelihood

# Provide information parameter space and initial sampling
Par_info = {}
Par_info['initial'] = 'prior'  		    # Use a informative prior distribution

# Define univariate priors
Par_info['prior'] = [
    norm(loc=2, scale=1),  		        # Normal prior with mean 2 and std 1
    gamma(a=1, loc=0, scale=1),  	    # Gamma prior with shape 1 and scale 1
    uniform(loc=-10, scale=20),  	    # Uniform prior between -10 and 10
    uniform(loc=-10, scale=20),  	    # Uniform prior between -10 and 10
    genpareto(c=1, loc=-4, scale=2)] 	# Generalized Pareto distribution

# Alternatively with multivariate normal prior (commented out here)
# Par_info['mu'] = [-2, -2, -2, -2, -2]
# Par_info['Sigma'] = 5 * np.eye(5)
# Par_info['prior'] = [multivariate_normal(mean = Par_info["mu"], cov = Par_info["Sigma"])]

# Boundary handling
Par_info['boundhandling'] = 'fold'  		    # Boundary handling method
Par_info['min'] = -10 * np.ones(DREAMPar['d'])  # Min values
Par_info['max'] = 10 * np.ones(DREAMPar['d'])  	# Max values

# Define name of function for posterior exploration
Func_name = 'mixture_lik.mixture_lik'

# Define method to use {'dream', 'dream_zs', 'dream_d', 'dream_dzs', 'mtdream_zs'}
method = 'dream'

# Set number of Markov chains and generations based on method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  	# Number of Markov chains
    DREAMPar['T'] = 10000  	# Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3  		# Number of Markov chains
    DREAMPar['T'] = 30000  	# Number of generations

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (500) * np.ones(DREAMPar['d'])  # Discrete steps

# Call the DREAM-Suite package
if __name__ == '__main__':
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info)
