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
## Example 14: Approx. Bayesian Comp. with watershed signatures & box-car likelihood  ##
##             Note: FDCFIT toolbox has a much larger built-in class of FDC functions ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Sadegh, M., J.A. Vrugt, H.V. Gupta, and C. Xu (2016), The soil water             ##
##       characteristic as new class of closed-form parametric expressions for the    ##
##       flow duration curve, J. Hydrol., 535, pp. 438–456                            ##
##   Lochbuhler, T., J.A. Vrugt, M. Sadegh, and N. Linde (2014), Summary statistics   ##
##       from training images as prior information in probabilistic inversion,        ##
##       Geophysical Journal International, 201, 157-171, doi:10.1093/gji/ggv008      ##
##   Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation using Markov ##
##       chain Monte Carlo simulation: DREAM_(ABC), Water Resources Research,         ##
##       doi:10.1002/2014WR015386.                                                    ##
##   Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model calibration and       ##
##       evaluation: Approximate Bayesian computation, Water Resources Research, 49,  ##
##       4335–4345, doi:10.1002/wrcr.20354.                                           ##
##   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), ##
##       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       ##
##       backward with Markov chain Monte Carlo simulation, Water Resources Research, ##
##       44, W00B09, doi:10.1029/2007WR006720                                         ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys
from scipy import io

current_dir = os.getcwd()                                       # Get the current working directory
from hmodel_summary_metrics import calc_signatures

parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite


## Problem settings defined by user
DREAMPar = {'d': 7, 	# Dimensionality of target
	        'lik': 22}  # ABC with summary metrics

# Provide information parameter space and initial sampling
Par_info = {'initial': 'latin',                                 # Latin hypercube sampling
            'boundhandling': 'reflect',                         # Explicit boundary handling
            'names': ['I_{\\rm max}', 'S_{\\rm u,max}', 'Q_{\\rm s,max}', '\\alpha_{\\rm e}', '\\alpha_{\\rm f}', 'K_{\\rm f}', 'K_{\\rm s}'],
            'min': [0.5, 10, 1e-4, 1e-6, -10, 1e-2, 1e-2],      # Minimum parameter values
            'max': [10, 1000, 100, 100, 10, 10, 150]}           # Maximum parameter values

# Define the name of the function for posterior exploration
Func_name = 'hmodel_summary_metrics.hmodel_summary_metrics'

# Load the French Broad data
daily_data = np.loadtxt('03451500.dly')

plugin = {}
plugin['data'] = {}
plugin['data']['wmp'] = 730
# idx_t1 = 731   # IN MATLAB
idx_t1 = 730
idx_t2 = daily_data.shape[0]
idx_p = np.arange(idx_t1 - plugin['data']['wmp'], idx_t2)

# The first two years are warm-up
plugin = {'idx': np.arange(plugin['data']['wmp'], len(idx_p) + 1),        
          'data': {'P': daily_data[idx_p, 3],       # Precipitation (mm/d)
                  'Ep': daily_data[idx_p, 4]},      # Potential evapotranspiration (mm/d)
          'y0': 1e-5 * np.ones(5),                  # Initial states
          'options': {'InitialStep': 1,             # Initial time step hmodel
                      'MaxStep': 1,                 # Maximum time step hmodel
                      'MinStep': 1e-3,              # Minimum time step hmodel
                      'RelTol': 1e-3,               # Relative tolerance 
                      'AbsTol': 1e-3,               # Absolute tolerance
                      'Order': 2},                  # Order of numerical solution
          'tout': np.arange(0, daily_data.shape[0] + 1)}

# Calculate summary metrics from the discharge data and define as prior distribution
Meas_info = {'S': calc_signatures(daily_data[idx_t1:idx_t2, 5]).T} 

# Optional settings
options = { 'DB': 'no',             # Diagnostic Bayes: ABC with summary metrics as prior
            'epsilon': [0.01],      # Epsilon of the noisy ABC implementation (illustration values)
            'parallel': 'yes',      # Run each chain on different core
            'modout': 'yes',  	    # Store model simulations
            'save': 'yes'}

# Define method to use {'dream', 'dream_zs', 'dream_d', 'dream_dzs', 'mtdream_zs'}
method = 'dream_zs'

if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  		    # Number of Markov chains
    DREAMPar['T'] = 10000  	        # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3  		        # Number of Markov chains
    DREAMPar['T'] = 5000  	        # Number of generations

# If the method is 'dream_d' or 'dream_dzs', set steps for discrete values
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (10000) * np.ones(DREAMPar['d'])

# Call the DREAM-Suite package
if __name__ == '__main__':
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, [], [], plugin)
