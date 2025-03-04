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
## Example 15: Approximate Bayesian Computation: Bivariate normal benchmark test      ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation using Markov ##
##       chain Monte Carlo simulation: DREAM_(ABC), Water Resources Research,         ##
##       doi:10.1002/2014WR015386.                                                    ##
##   Sadegh, M., and J.A. Vrugt (2013), Bridging the gap between GLUE and formal      ##
##       statistical approaches: approximate Bayesian computation, Hydrology and      ##
##       Earth System Sciences, 17, 4831–4850                                         ##
##   Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model calibration and       ##
##       evaluation: Approximate Bayesian computation, Water Resources Research, 49,  ##
##       4335–4345, doi:10.1002/wrcr.20354                                            ##
##   Turner, B.M., and P.B. Sederberg (2013), Approximate Bayesian computation with   ##
##       differential evolution, Journal of Mathematical Psychology, In Press.        ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys
import matplotlib.pyplot as plt
from scipy import io

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from DREAM_Suite_functions import genparset, X_unnormalize      # Import functions


# Define how many bivariate normal distributions
n_pairs = 10

# Problem settings defined by user
DREAMPar = {'d': 2 * n_pairs,                       # Dimension of the problem
            'lik': 22}                              # ABC informal likelihood function

# Provide information parameter space and initial sampling
Par_info = {'initial': 'latin',         	        # Latin hypercube sampling
            'boundhandling': 'fold',    	        # Explicit boundary handling
            'names': ['\\mu^{1}_{x}','\\mu^{2}_{x}','\\mu^{3}_{x}','\\mu^{4}_{x}','\\mu^{5}_{x}','\\mu^{6}_{x}','\\mu^{7}_{x}','\\mu^{8}_{x}','\\mu^{9}_{x}','\\mu^{10}_{x}', \
                      '\\mu^{1}_{y}','\\mu^{2}_{y}','\\mu^{3}_{y}','\\mu^{4}_{y}','\\mu^{5}_{y}','\\mu^{6}_{y}','\\mu^{7}_{y}','\\mu^{8}_{y}','\\mu^{9}_{y}','\\mu^{10}_{y}'],
            'min': np.zeros(2 * n_pairs),  	        # If 'latin', min values
            'max': 10 * np.ones(2 * n_pairs),       # If 'latin', max values
            'norm': 1}                              # sample in normalized space [0-1]

# Define name of function for posterior exploration
Func_name = 'ABC_binormal.ABC_binormal'

# Create the observed summary metrics - the mean (mu) of ten bivariate normals
Meas_info = {'S': Par_info['min'] + np.random.rand(DREAMPar['d']) * (Par_info['max'] - Par_info['min'])}

# Optional settings for rho (distance function)
options = {'rho': 'lambda X, Y: np.sqrt(1 / 20 * np.sum((X - Y)**2))',
           'print': 'no'}

# Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs'

# Set up parameters for different methods
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10          # Number of Markov chains
    DREAMPar['T'] = 3000        # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3           # Number of Markov chains
    DREAMPar['T'] = 50000       # Number of generations

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (100) * np.ones(DREAMPar['d'])  # Number of discrete steps

# Call the DREAM-Suite package
if __name__ == '__main__':
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options)

    # Now plot results
    parset = genparset(chain)                   # Create matrix of 3D chain trajectories
    parset = parset[50000:-1, :]                # Apply burn-in
    P_un = X_unnormalize(parset[:, :DREAMPar['d']], Par_info)
    id_1 = np.arange(0, n_pairs)                # Indices for x-coordinates
    id_2 = np.arange(n_pairs, 2 * n_pairs)      # Indices for y-coordinates

    # Plot the results
    plt.figure(100)
    for i in range(0,n_pairs):
        if i < n_pairs - 1:
            plt.plot(P_un[:, id_1[i]], P_un[:, id_2[i]], 'r.')
            plt.plot(Meas_info['S'][id_1[i]], Meas_info['S'][id_2[i]], 'kx')
        else:
            plt.plot(P_un[:, id_1[i]], P_un[:, id_2[i]], 'r.', label = 'MCMC samples')
            plt.plot(Meas_info['S'][id_1[i]], Meas_info['S'][id_2[i]], 'kx', label = 'Observed data')

    plt.xlabel("X-coordinate [= means] of bivariate distributions")
    plt.ylabel("Y-coordinate [= means] of bivariate distributions")
    plt.title('DREAM-Suite inferred (red dots) and true (black cross) means of bivariate normals')
    plt.axis('square')
    plt.legend()
    plt.show()
