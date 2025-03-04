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
## Example 34: Inverse modeling of Parlange's semi-analytic infiltration equation     ##
##             using HYDRUS simulated infiltration curves                             ##
##                                                                                    ##
## Check the following papers                                                         ##
##  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and               ##
##      H. Vereecken (2023), The time validity of Philip's two-term infiltration      ##
##      equation: An elusive theoretical quantity? Vadose Zone Journal, e20309,       ##
##      pp. 1-25, https://doi.org/10.1002/vzj2.20309                                  ##
##  Vrugt, J.A. and Y. Gao (2022), On the three-parameter infiltration equation of    ##
##      Parlange et al. (1982): Numerical solution, experimental design, and          ##
##      parameter estimation, Vadose Zone Journal, 21:e20167, pp. 1-25,               ##
##      https://doi.org/10.1002/vzj2.20167                                            ##
##                                                                                    ##
## ---------------------------------------------------------------------------------- ##

import numpy as np
import os, sys
import matplotlib.pyplot as plt
import scipy.io

current_dir = os.getcwd()                                               # Get the current working directory
from preprocess_data import preprocess_data

parent_dir = os.path.abspath(os.path.join(current_dir, '..'))           # Go up one directory
sys.path.append(parent_dir)                                             # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	            # Add miscellaneous directory to Python path
from DREAM_Suite_functions import X_unnormalize, genparset           # Import functions

# Set up plugin and model configuration
plugin = {'model_setup': 2}

# Problem settings defined by user
DREAMPar = {'lik': 11}  # Gaussian likelihood

# Provide information for parameter space and initial sampling
Par_info = {'initial': 'latin',  		    # Latin hypercube sampling
	 	    'boundhandling': 'reflect'}  	# Explicit boundary handling

# Define name of the function for posterior exploration
Func_name = 'Haverkamp_I_patch.Haverkamp_I_patch'

# Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs'

# Set number of chains and generations based on the method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  		# Markov chains
    DREAMPar['T'] = 10000  		# generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3  			# Markov chains
    DREAMPar['T'] = 5000  		# generations

# Set parameter bounds based on model setup
if plugin['model_setup'] == 1:
    Par_info['names'] = ['S','K_{\\rm s}','\\beta']
    Par_info['min'] = [0, 0, 0]  	    # Minimum parameter values
    Par_info['max'] = [20, 50, 2]  	    # Maximum parameter values
elif plugin['model_setup'] == 2:
    Par_info['names'] = ['S','K_{\\rm s}','K_{\\ rm i}','\\beta']
    Par_info['min'] = [0, 0, 0, 0]  	# Minimum parameter values
    Par_info['max'] = [10, 50, 2, 2]  	# Maximum parameter values

# Dimension of the problem
DREAMPar['d'] = len(Par_info['min'])

# Optional settings
options = {	'modout': 'yes',  	# Return model simulations (yes/no)?
    		'parallel': 'no', 	# Run each chain on a different core
    		'save': 'yes',  	# Save workspace during run
    		'print': 'no'} 		# No figures printed to screen

# Load measured infiltration data
data = scipy.io.loadmat('HYDRUS_1D_Data.mat')

# Preprocess the data
n_soil = 12  				    # How many soil types?
plugin['n'] = 100  			    # interpolated (t_meas, I_meas)
interpolation_method = 2  		# Interpolation method
I_max = 5  				        # Maximum infiltration in cm
interpolation_scheme = 'linear' # Interpolation scheme?
Parameters = data['Parameters'] # extract parameters from MATLAB mat file

# Assuming 'preprocess_data' is a function you've already defined
data_new, true_pars = preprocess_data(data['data'], Parameters, n_soil, plugin['n'], interpolation_method, interpolation_scheme, I_max)

# Initialize matrices to hold results
p95 = np.nan * np.zeros((n_soil, 5, DREAMPar['d']))
ML = np.nan * np.zeros((n_soil, DREAMPar['d']))

if __name__ == '__main__':
    # Loop through each soil type
    for soil_type in range(0, 12):
        dat = data_new[soil_type]    		        # Unpack data for current soil type
        Meas_info = {'Y': dat[:plugin['n'], 1]}  	# Measurement data (cm)
        plugin['t'] = dat[:plugin['n'], 0]  	    # Measurement times (hours)

        # Run DREAM-Suite package
        chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, [], [], plugin)

        # Unpack Markov chains
        parset = genparset(chain)  			        # Assuming genparset is a function defined elsewhere
        N = chain.shape[0]  			            # Number of samples per chain
        P = parset[N // 2:, :DREAMPar['d'] + 2]  	# Posterior samples (after burn-in)
        M = P.shape[0]

        # Find "best" solution based on maximum likelihood
        ii = np.argmax(np.sum(P[:, DREAMPar['d'] : DREAMPar['d'] + 2], axis=1))
        ML[soil_type, :] = P[ii, :DREAMPar['d']]

        # Calculate 95% confidence intervals and statistics
        for j in range(DREAMPar['d']):
            a = np.sort(P[:, j])
            p95[soil_type, :, j] = [
                a[int(0.025 * M)],
                np.mean(a),
                np.median(a),
                a[int(0.975 * M)],
                np.std(a) ]

        # Save the results to a file
        evalstr = f"DREAM_Suite_{soil_type}_{interpolation_scheme}_{plugin['n']}.npy"

        if os.name == 'nt':  			            # If running on Windows
            os.system(f"rename DREAM_Suite.npy {evalstr}")
        elif os.name == 'posix':  			        # If running on Unix-based systems (Linux/Mac)
            os.system(f"mv {evalstr}")
        else:
            raise OSError("Example_34: Unsupported operating system")

    # Save final results
    scipy.io.savemat('HYDRUS_results.mat', {'ML': ML, 'p95': p95, 'data_new': data_new, 'true_pars': true_pars,
        'Parameters': Parameters, 'data': data, 'I_max': I_max, 'interpolation_scheme': interpolation_scheme, 'interpolation_method': interpolation_method, 'n_soil': n_soil})
