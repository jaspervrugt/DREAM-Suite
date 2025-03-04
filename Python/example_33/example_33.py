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
## Example 33: Two-dimensional target distribution: from Haario et al. 1999           ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, https://doi.org/10.1016/j.envsoft.2015.08.013 ##
##  Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M.    ##
##      Hyman (2009), Accelerating Markov chain Monte Carlo simulation by             ##
##      differential evolution with self-adaptive randomized subspace sampling,       ##
##      International Journal of Nonlinear Sciences and Numerical Simulation, 10(3),  ##
##      271-288                                                                       ##
##   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), ##
##       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       ##
##       backward with Markov chain Monte Carlo simulation, Water Resources Research, ##
##       44, W00B09, https://doi.org/10.1029/2007WR006720                             ##
##   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian (2003), A Shuffled Complex  ##
##      Evolution Metropolis algorithm for optimization and uncertainty assessment of ##
##      hydrologic model parameters, Water Resour. Res., 39 (8), 1201,                ##
##      https://doi.org/10.1029/2002WR001642                                          ##
##   Haario, H., E. Saksman, and J. Tamminen (1999), Adaptive proposal distribution   ##
##      for random walk Metropolis algorithm. Computational Statistics, 14, 375–395,  ##
##      https://doi.org/10.1007/s001800050022                                         ##
##                                                                                    ##
## ---------------------------------------------------------------------------------- ##

import numpy as np
import os, sys
import matplotlib.pyplot as plt

current_dir = os.getcwd()                                               # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))           # Go up one directory
sys.path.append(parent_dir)                                             # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	            # Add miscellaneous directory to Python path
from DREAM_Suite_functions import X_unnormalize, genparset              # Import functions

# Problem settings defined by user
DREAMPar = {'d': 2,  		# Dimension of the problem
			'lik': 1}  	    # Model output is a likelihood

# Determine the parameter ranges
Par_info = {'initial': 'latin',  			# Latin hypercube sampling
			'boundhandling': 'fold',  	    # Explicit boundary handling
			'min': [ -18, -3 ],  			# minimum of parameter values
			'max': [  18,  3 ],             # maximum of parameter values
            'norm': 1}    		            # sample in normalized space (= test)
Parinf = Par_info.copy()

## Define name of function (.m file) for posterior exploration
Func_name = 'rectangle_lik.rectangle_lik'

# Additional options for DREAM Package
options = {	'modout': 'no',  		# Return model simulations (yes/no)?
			'parallel': 'no',     	# Run each chain on a different core
			'save': 'yes'}  		# Save workspace DREAM during run

# Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs'

# Switch between methods
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  		# Number of Markov chains
    DREAMPar['T'] = 10000  		# Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3  			# Number of Markov chains
    DREAMPar['T'] = 20000  		# Number of generations

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = [7200, 600]  # Number of discrete steps

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info)

    # Make a 2D matrix from 3D chain trajectories
    P = genparset(chain)

    # Extract the last 30,000 samples from the samples
    Post = P[-30000:, :]

    # If sampled in normalized space, back-transform to original space
    if 'norm' not in Par_info:
        Par_info['norm'] = 0  # sampled in regular space

    if Par_info['norm'] == 1:
        # Back-transform the parameters
        Par_info['minun'] = np.array(Parinf['min'])
        Par_info['maxun'] = np.array(Parinf['max'])
        Par_info['min'] = np.array([0, 0])
        Par_info['max'] = np.array([1, 1])
        Post[:, 0:2] = X_unnormalize(Post[:, 0:2], Par_info)

    # Look only at the first parameter as the second is equal to the prior (inconsequential)
    edges = np.concatenate([
        np.arange(-18, -0.49, 0.25), 
        np.arange(-0.3, 0.31, 0.2), 
        np.arange(0.5, 18.01, 0.25)
    ])

    # Make a histogram of the posterior samples
    N, edges = np.histogram(Post[:, 0], bins = edges, density = True)

    # Midpoints of bins
    bin = 0.5 * (edges[:-1] + edges[1:])

    # Check if the empirical density integrates to one (should be sufficiently close)
    integral = np.trapezoid(N, bin)
    print(f"Integral of empirical density: {integral}")

    # Now plot a histogram of the marginal distribution according to MCMC samples
    plt.figure(100)
    #plt.bar(bin, N, width = np.diff(edges), color = 'b', align = 'edge', label = 'MCMC-estimate')
    plt.bar(bin, N, width = np.diff(edges), color = 'b', align = 'edge', label = 'MCMC-estimate')
    # JAV: Need to check histogram as it has moved a bit to the right - should match with MATLAB and true distribution

    # Now plot true target density
    x = np.linspace(-18, 18, 3601)
    f = lambda x: 35 * ((x >= -0.5) & (x <= 0.5)) + 1
    Z = np.trapezoid(f(x), x)

    plt.plot(x, f(x) / Z, 'r', linewidth = 2, label = 'True target')

    plt.legend(['MCMC-estimate', 'True target'], fontsize = 16)
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.show()

    # Check percentage of posterior samples outside S: [-0.5, 0.5]
    M = P.shape[0]
    out_MCMC = 100 * np.sum((P[:, 0] < -0.5) | (P[:, 0] > 0.5)) / M

    # True probability mass outside S
    out_true = 100 * 35 / (35 + 36)

    print(f"Percentage of MCMC samples outside S: {out_MCMC:.2f}%")
    print(f"Percentage of true target mass outside S: {out_true:.2f}%")

