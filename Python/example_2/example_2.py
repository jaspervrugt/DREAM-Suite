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
## Example 2: N(µ,Σ) with µ = 0 and Σ = dxd matrix with 1:d on main diagonal and      ##
##                   0.5 correlation between dimensions                               ##
## Check: func_normal for details                                                     ##
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


# Define problem settings
DREAMPar = {}
DREAMPar['d'] = 10   		# Dimension of the problem
DREAMPar['thinning'] = 10  	# Only store every 2nd sample
DREAMPar['lik'] = 2  		# Model output is log-likelihood

# Provide information about the parameter space and initial sampling
Par_info = {}
Par_info['initial'] = 'latin'  			        # Latin hypercube sampling
Par_info['min'] = -15 * np.ones(DREAMPar['d'])  # Min values for latin sampling
Par_info['max'] = 15 * np.ones(DREAMPar['d'])  	# Max values for latin sampling
Par_info['boundhandling'] = 'none'

# Define function name for posterior exploration
Func_name = 'normal_lik.normal_lik'

# Define method to use {'dream', 'dream_zs', 'dream_d', 'dream_dzs', 'mtdream_zs'}
method = 'dream_zs'
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = DREAMPar['d'] 	# Number of Markov chains
    DREAMPar['T'] = 10000  		    # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 5     		    # Number of Markov chains
    DREAMPar['T'] = 15000  		    # Number of generations

if method in ['dream_d', 'dream_dzs']:
    Par_info['min'] = -50 * np.ones(DREAMPar['d'])  	# Min values for discrete sampling
    Par_info['max'] = 50 * np.ones(DREAMPar['d'])  	    # Max values for discrete sampling
    Par_info['steps'] = 1000 * np.ones(DREAMPar['d'])  	# Discrete steps

# Options for saving and printing
options = {'save': 'yes', 'print': 'no'}

if __name__ == '__main__':
    # Call the eDREAM package function
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, [], options)

    # Example of generating mock `chain` data (replace with actual data from DREAM-Suite)
    parset = genparset(chain)  				            # Simulating MCMC chain

    # Check the mean and variance of MCMC samples
    P = parset[-25000:, :DREAMPar['d']]  			    # Last 25000 samples, taking first 100 dimensions

    # Plot the mean of MCMC samples
    plt.figure(100)
    plt.plot(np.arange(1,DREAMPar['d'] + 1), np.mean(P, axis = 0), 'b+', markersize = 8, linewidth = 2)
    plt.axhline(0, color = 'r', linewidth = 2)
    plt.xlabel('Dimension', fontsize = 16)
    plt.ylabel(r'Mean, $\mu_{i}$', fontsize = 16)
    plt.legend(['Estimated', 'True'], loc = 'best', fontsize = 16)
    plt.grid(True)
    plt.axis('square')
    plt.title('Target mean against estimates from MCMC', fontsize = 16)

    # Plot the variance of MCMC samples
    plt.figure(101)
    plt.plot(np.arange(1,DREAMPar['d'] + 1),np.diag(np.cov(P.T)), 'b+', markersize = 8, linewidth = 2)      # Variance is diagonal of covariance matrix
    plt.plot(np.arange(1,DREAMPar['d'] + 1),np.arange(1,DREAMPar['d'] + 1), color = 'r', linewidth = 2)  	# True variance (set to DREAMPar['d'] for example)
    plt.xlabel('Dimension', fontsize = 16)
    plt.ylabel(r'Variance, $\sigma^{2}_{i}$', fontsize = 16)
    plt.legend(['Estimated', 'True'], loc = 'best', fontsize = 16)
    plt.grid(True)
    plt.axis('square')
    plt.title('Target variance against estimates from MCMC', fontsize = 16)

    # Show the plots
    plt.show()

    # Posterior samples with log-prior and log-likelihood
    P = parset[-10000:, :] 

    # Now we illustrate GAME sampling (marginal likelihood estimation = integral of posterior pdf)
    Z, logZ, gmix = GAME_sampling(P, 'is', DREAMPar, Func_name, Par_info)

    # Print marginal likelihood
    print("Marginal likelihood equal to:", Z)
    print("The expected integral should equal to one")
