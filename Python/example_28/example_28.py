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
## Example 28: Conceptual model parameter and rainfall estimation                     ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten (2005),     ##
##       Improved treatment of uncertainty in hydrologic modeling: Combining the      ##
##       strengths of global optimization and data assimilation, Water Resources      ##
##       Research, 41, W01017, doi:10.1029/2004WR003059                               ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys

current_dir = os.getcwd()                                       # Get the current working directory
from rainfall_runoff import load_data_dly, check_rainfall

parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from DREAM_Suite_functions import nuis_var_names	            # Import functions

# Define the problem settings
DREAMPar = {'nCR': 10,      # Increase subspace sampling (default is 3)
		    'delta': 2,     # Delta is 2
		    'lik': 16}      # Model output is simulation: Laplacian likelihood

# Parameter names of hmodel
par_names = ['I_{\\rm max}', 'S_{\\rm u,max}', 'Q_{\\rm s,max}', 
             '\\alpha_{\\rm e}', '\\alpha_{\\rm f}', 'K_{\\rm f}', 'K_{\\rm s}']

# Index:      1   2     3    4   5    6     7
# Parname:  Imax Smax Qsmax alE alF Kfast Kslow
fpar_mod = [1, 100, 10, 100, 0, 2, 70]
parmin_mod = [0.5, 10, 0, 1e-6, -10, 1e-2, 1e-2]    # Latin min values
parmax_mod = [10, 1000, 100, 100, 10, 10, 150]      # Latin max values

# Provide information on parameter space and initial sampling
Par_info = {'initial': 'latin',  			        # Latin hypercube sampling
            'boundhandling': 'reflect'}  		    # Explicit boundary handling

# Global variable for likelihood functions
LV = {}

# Case statement for likelihood function
if DREAMPar['lik'] == 13:       ## Normal distribution
    # index:      1  2  3   4
    # parname:   s0 s1 phi1 phi2
    fpar_nuis = [0.1, 0, 0, 0]
    parmin_nuis = [0, 0, 0, 0]
    parmax_nuis = [1, 1, 1, 1]
    LV['filename'] = 'Normal'
    id_nuis = [2]
elif DREAMPar['lik'] == 16:     ## Laplace distribution
    # index:     1   2    3
    # parname:  s0  s1  phi1
    fpar_nuis = [0.1, 0, 0]
    parmin_nuis = [0, 0, 0]
    parmax_nuis = [1, 1, 1]
    LV['filename'] = 'Laplace'
    id_nuis = [0, 1]
elif DREAMPar['lik'] == 17:     ## SST distribution
    # index:          1   2    3   4    5    6
    # parname:       s0  s1  nu   xi phi1  phi2
    fpar_nuis = [0.1, 0, 1e10, 1, 0, 0]
    parmin_nuis = [0, 0, 2, 0.1, 0, 0]
    parmax_nuis = [1, 1, 100, 10, 1, 1]
    LV['filename'] = 'SL'
    id_nuis = [2, 3, 4]
elif DREAMPar['lik'] == 44:     ## GL+
    # index:     1   2    3   4    5    6
    # parname:   s0  s1  beta xi phi1  phi2
    fpar_nuis = [0.1, 0, 0, 1, 0, 0]
    parmin_nuis = [0, 0, -1, 0.1, 0, 0]
    parmax_nuis = [1, 1, 1, 10, 1, 1]
    LV['filename'] = 'GL_plus'
    id_nuis = [2, 3, 4]
elif DREAMPar['lik'] == 45:     ## SGT (universal) likelihood
    # index:      1   2    3   4   5    6    7
    # parname    s0  s1  labda p   q  phi1 phi2
    fpar_nuis = [0.1, 0, 0, 2, 1e10, 0, 0]
    parmin_nuis = [0, 0, -1, 0.5, 2, 0, 0]
    parmax_nuis = [1, 1, 1, 100, 100, 1, 1]
    LV['filename'] = 'UL'
    id_nuis = [0, 2, 3, 4, 6]

# Define name of function for posterior exploration
Func_name = 'rainfall_runoff.rainfall_runoff'

data = load_data_dly('09497500')  		                        # Example data loading function
id = [0, 1, 2, 5, 4, 3]
data = data[:, id]                                              # Reorder data
wmp = 64                                                        # Warm-up period
idx_t1 = 64                                                     # Start of training period
idx_t2 = 795                                                    # Maximum time
idx_p = np.arange(idx_t1 - wmp, idx_t2)

# Rainfall storm detection and parameter adjustments
plugin = {'Tmax': idx_t2,                                       # Maximum time
          'nmod': len(parmin_mod), 		                        # Number of model parameters
          'P': data[idx_p, 5],  	                            # Daily rainfall data
          'Y': data[idx_p, 3],  	                            # Daily discharge data
          'Ep': data[idx_p, 4]}  	                            # Daily PET data

plugin, mult, mult_names, fpar_mult = check_rainfall(plugin)
parmin_mod = np.concatenate([parmin_mod, mult['min']])          # Add multiplier minimum
parmax_mod = np.concatenate([parmax_mod, mult['max']])          # Add multiplier maximum
par_names.extend(mult_names)  				                    # Add multiplier names
fpar_mod.extend(fpar_mult)  				                    # Add default multipliers

# Define nuisance variable names
nuis_names = nuis_var_names(DREAMPar)
names = par_names + nuis_names

# Model parameters and nuisance variables
n_modpar = len(fpar_mod)

LV['id_vpar'] = list(range(n_modpar)) + [x + n_modpar for x in id_nuis]
LV['fpar'] = fpar_mod + fpar_nuis

parmin = np.concatenate([parmin_mod, parmin_nuis])
parmax = np.concatenate([parmax_mod, parmax_nuis]) 

# Min/max values of parameters
Par_info['min'] = [parmin[i] for i in LV['id_vpar']]
Par_info['max'] = [parmax[i] for i in LV['id_vpar']]
Par_info['names'] = [names[i] for i in LV['id_vpar']]
DREAMPar['d'] = len(Par_info['min'])

# Define hmodel settings for numerical simulation
plugin['y0'] = 1e-5 * np.ones(5)  			            # Initial states (five reservoirs almost empty)
plugin['tout'] = np.arange(0, len(idx_p) + 1)           # Output simulation times of hmodel
plugin['idx'] = np.arange(idx_t1, idx_t2 + 1)           # Indices of simulated record
plugin['data'] = {'P': plugin['P'],
                  'Ep': plugin['Ep']}

# Integration settings for hmodel
plugin['options'] = {   'InitialStep': 0.2,		        # Initial time-step (days)
                        'MaxStep': 1,  			        # Maximum time-step (days)
                        'MinStep': 1e-5,  		        # Minimum time-step (days)
                        'RelTol': 1e-5, 		        # Relative tolerance
                        'AbsTol': 1e-5,                 # Absolute tolerances (mm)
                        'Order': 2}  			        # 2nd order accurate method (Heun)

# Specify discharge measurement error
Meas_info = {'Y': data[64:idx_t2, 3],   	            # Define measured discharge
            'Sigma': 1/10 * data[64:idx_t2, 3] + 1/100}

# Define optional settings
options = {'modout': 'yes',  	    # Return model simulations of samples?
            'parallel': 'yes', 		# Evolve each chain on a different node?
            'save': 'yes',  		# Save workspace during DREAM process
            'restart': 'no'}        # No restart run

# Define method to use
method = 'dream_zs'

# Set Markov chain parameters
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 50  		                                # Number of Markov chains
    DREAMPar['T'] = 20000  		                                # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 5  			                                # Number of Markov chains
    DREAMPar['T'] = 10000  		                                # Number of generations

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = np.ones(DREAMPar['d']) * (1000)  # Number of discrete steps

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, LV, [], plugin)
