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
## Example 23: SAC-SMA model: Limits of Acceptability                                 ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. and K.J. Beven (2018), Embracing equifinality with efficiency:       ##
##       Limits of Acceptability sampling using the DREAM_{(LOA)} algorithm,  Journal ##
##       of Hydrology,  559 , pp. 954-971, doi:10.1016/j.jhydrol.2018.02.026          ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Vrugt, J.A., Gupta H.V., Dekker, S.C., Sorooshian, S., Wagener, T. and W.        ##
##       Bouten (2006), Application of stochastic parameter optimization to the       ##
##       Sacramento Soil Moisture Accounting model, Journal of Hydrology, 325, pp.    ##
##       288-307                                                                      ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys
from scipy.stats import norm

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

# Problem settings defined by user
DREAMPar = {'d': 13,  			# Dimension of the problem
    		'thinning': 1,      # Only store each 5th sample
    		'lik': 23}          # Limits of acceptability likelihood

# Provide information parameter space and initial sampling
Par_info = {'initial': 'latin',                													# Initial sample from Latin hypercube sampling
    		'boundhandling': 'reflect',        													# Explicit boundary handling
               #    'uzfwm', 'uztwm', 'lzfpm', 'lzfsm', 'lztwm', 'zperc', 'rexp', 'uzk', 'pfree', 'lzpk', 'lzsk', 'acm', 'kf'     
               #       S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, nu_p, nu_s, a_cmax, kf
    		'names': ['S1_{\\rm Fmax}', 'S1_{\\rm Tmax}', 'S2_{\\rm FPmax}', 'S2_{\\rm FSmax}', 'S2_{\\rm Tmax}', '\\alpha', '\\psi', 'k_{\\rm i}', '\\kappa', '\\nu_{\\rm p}', '\\nu_{\\rm s}', 'a_{\\rm cmax}', 'k_{\\rm f}'],    # Parameter names
    		'min': [10, 10, 10, 10, 10, 1, 1, 0.01, 0.05, 1e-3, 1e-3, 0.05, 0.0],  	# Minimum values
    		'max': [500, 500, 1000, 1000, 500, 250, 5, 1000, 0.95, 0.25, 0.25, 0.95, 1.0]}  	# Maximum values

# Define model name
Func_name = 'sacsma.sacsma'

# Loading data (here we assume `bound.txt` is a file)
bound = np.loadtxt('bound.txt')

# Define warm-up period and indices for measured discharge record
plugin = {}
plugin['data'] = {}
plugin['data']['wmp'] = 64
idx_t1 = 64
idx_t2 = 795
idx_p = np.arange(idx_t1 - plugin['data']['wmp'], idx_t2)

# The 65 days are warmp-up
plugin = {'idx': np.arange(plugin['data']['wmp'], len(idx_p) + 1),        
          'data': {'P': np.sum(bound[idx_p, 5:9], axis = 1),    # Rainfall (mm/d)
                  'Ep': bound[idx_p, 4]},           # Potential evapotranspiration (mm/d)
          'y0': 1e-5 * np.ones(9),                  # Initial states
          'options': {'InitialStep': 0.2,           # Initial time step hmodel
                      'MaxStep': 1,                 # Maximum time step hmodel
                      'MinStep': 1e-5,              # Minimum time step hmodel
                      'RelTol': 1e-3,               # Relative tolerance 
                      'AbsTol': 1e-3,               # Absolute tolerance
                      'Order': 2},                  # Order of numerical solution
          'tout': np.arange(0, len(idx_p) + 1)}

# Measured streamflow data (converted from m3/s to mm/d)
F = 1944 * (1000 * 1000) / (1000 * 60 * 60 * 24)
Meas_info = {'S': 1 / F * bound[idx_t1:idx_t2, 3]}

# Optional settings
options = {	'modout': 'yes',                	# Return model simulations (yes/no)?
            'parallel': 'no',                  # Run chains in parallel 
    		'epsilon': 1.2 * Meas_info['S']}   	# Limits of acceptability

# Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'mtdream_zs'

# Set parameters for the chosen method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10               # Number of Markov chains
    DREAMPar['T'] = 5000             # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3                # Number of Markov chains
    DREAMPar['T'] = 25000            # Number of generations

# Discrete steps for methods with discrete sampling
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (500) * DREAMPar['d']  # Number of discrete steps for each parameter

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, [], [], plugin)
