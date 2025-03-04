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
## Example 9: Application of spectral likelihood function to conceptual               ##
##            watershed models using measured discharge data                          ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
##       distribution-adaptive likelihood functions: Generalized and universal        ##
##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
##       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             ##
##       https://www.sciencedirect.com/science/article/pii/S002216942201112X          ##
##   Schoups, G., and J.A. Vrugt (2010), A formal likelihood function for parameter   ##
##       and predictive inference of hydrologic models with correlated,               ##
##       heteroscedastic and non-Gaussian errors, Water Resources Research, 46,       ##
##       W10531, doi:10.1029/2009WR008933                                             ##
##   Schoups, G., J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption   ##
##       of accuracy and efficiency of Markov Chain Monte Carlo simulation by         ##
##       inaccurate numerical implementation of conceptual hydrologic models, Water   ##
##       Resources Research, 46, W10530, doi:10.1029/2009WR008648                     ##
##   Montanari, A., and E. Toth (2007), Calibration of hydrological models in the     ##
##       spectral domain: An opportunity for scarcely gauged basins?, Water Resources ##
##       Research, 43, W05434, doi:10.1029/2006WR005184                               ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

## Problem settings defined by user
DREAMPar = {'d': 7,
            'lik': 15}    # 15: Spectral (Whittle) likelihood

# Define the name of the function for posterior exploration
Func_name = 'hmodel.hmodel'

# Parameter info
Par_info = {'min': [0.5, 10, 0, 1e-6, -10, 0, 0],       # Minimum values
            'max': [10, 1000, 100, 100, 10, 10, 150],   # Maximum values
            'names': ['I_{\\rm max}', 'S_{\\rm u,max}', 'Q_{\\rm s,max}', '\\alpha_{\\rm e}', '\\alpha_{\\rm f}', 'K_{\\rm f}', 'K_{\\rm s}'],
            'initial': 'latin',  		                # Latin hypercube sampling
            'boundhandling': 'reflect',  	            # Explicit boundary handling
            'norm': 0 }                                 # Normalized space?

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
          'tout': np.arange(0, len(idx_p) + 1)}

# Measurement information
Meas_info = {'Y': daily_data[idx_t1:idx_t2, 5]}  	    # Adjusted to match the MATLAB code (indexing)

# Define method to use
method = 'dream_zs'

# Set Markov chain parameters based on the method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10
    DREAMPar['T'] = 2500
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3
    DREAMPar['T'] = 20000
elif method == 'dream_kzs':
    DREAMPar['N'] = 3
    DREAMPar['T'] = 5000
    DREAMPar['M'] = 24

# Additional setup for 'dream_d' or 'dream_dzs'
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (1000) * np.ones(DREAMPar['d'])

# Optional settings
options = {'modout': 'no',
           'parallel': 'no',
           'save': 'yes',
           'restart': 'no'}

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, [], [], plugin)
