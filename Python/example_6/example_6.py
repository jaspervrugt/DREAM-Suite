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
## Example 6: Distribution-adaptive likelihood functions applied to conceptual        ##
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
##   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), ##
##       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       ##
##       backward with Markov chain Monte Carlo simulation, Water Resources Research, ##
##       44, W00B09, doi:10.1029/2007WR006720                                         ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys
import matplotlib.pyplot as plt

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from DREAM_Suite_functions import nuis_var_names	            # Import functions


## Problem settings defined by user
DREAMPar = {}
DREAMPar['lik'] = 44;  	# 13: Normal likelihood function with AR(2)-model of residuals
                    	# 14: Generalized likelihood function (Obsolete)
                    	# 16: Laplace likelihood function with AR(2)-model of residuals
                    	# 17: Skewed Student t likelihood function
                    	# 44: Generalized likelihood function PLUS
                    	# 45: Universal likelihood function
## Please refer to Vrugt et al. (2022) for theory on likelihood functions

# Define the name of the function for posterior exploration
Func_name = 'hmodel.hmodel'

# Parameter names of the hmodel
par_names = ['I_{\\rm max}', 'S_{\\rm u,max}', 'Q_{\\rm s,max}', '\\alpha_{\\rm e}', '\\alpha_{\\rm f}', 'K_{\\rm f}', 'K_{\\rm s}']
fpar_mod = [1, 100, 10, 100, 0, 2, 70]		    # Default values
parmin_mod = [0.5, 10, 0, 1e-6, -10, 0, 0]  	# Minimum values
parmax_mod = [10, 1000, 100, 100, 10, 10, 150]  # Maximum values

# Provide information for parameter space and initial sampling
Par_info = {}
Par_info['initial'] = 'latin'  		    # Latin hypercube sampling
Par_info['boundhandling'] = 'reflect'  	# Explicit boundary handling
#Par_info['norm'] = 1			        # Normalized space?

# Likelihood function - summary variable
LV = {}

# Handling likelihood function selection with a switch-like structure
if DREAMPar['lik'] == 13:   ## Normal distribution
    # index:      0   1    2   3
    # parname:   s0  s1  phi1 phi2    
    fpar_nuis = [0.1, 0, 0, 0]
    parmin_nuis = [0, 0, 0, 0]
    parmax_nuis = [1, 1, 1, 1]
    LV['filename'] = 'Normal'
    id_nuis = [2]
elif DREAMPar['lik'] == 14:  ## GL (OBSOLETE)
    # index:      0   1    2    3   4   5   6    7    8    9  10  
    # parname:  std0 std1 beta xi  mu1 phi1 phi2 phi3 phi4 K lambda
    fpar_nuis = [0.1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
    parmin_nuis = [0, 0, -1, 0.1, 0, 0, 0, 0, 0, 0, 0.1]
    parmax_nuis = [1, 1, 1, 10, 100, 1, 1, 1, 1, 1, 1]
    LV['filename'] = 'GL'
    id_nuis = [0, 1, 2, 3, 5]
elif DREAMPar['lik'] == 16:  ## Laplace distribution
    # index:      0   1   2
    # parname:   s0  s1  phi1
    fpar_nuis = [0.1, 0, 0]
    parmin_nuis = [0, 0, 0]
    parmax_nuis = [1, 1, 1]
    LV['filename'] = 'Laplace'
    id_nuis = [0]
elif DREAMPar['lik'] == 17:  ## SST (Skewed Student-t)
    # index:      0   1   2    3   4    5
    # parname:   s0  s1  nu   xi phi1  phi2
    fpar_nuis = [0.1, 0, 1e10, 1, 0, 0]
    parmin_nuis = [0, 0, 2, 0.1, 0, 0]
    parmax_nuis = [1, 1, 100, 10, 1, 1]
    LV['filename'] = 'SL'
    id_nuis = [2, 3, 4]
elif DREAMPar['lik'] == 44:  ## GL+
    # index:      0   1   2    3   4    5
    # parname:   s0  s1  beta xi phi1  phi2
    fpar_nuis = [0.1, 0, 0, 1, 0, 0]
    parmin_nuis = [0, 0, -1, 0.1, 0, 0]
    parmax_nuis = [1, 1, 1, 10, 1, 1]
    LV['filename'] = 'GL_plus'
    id_nuis = [0, 2, 3, 4]
elif DREAMPar['lik'] == 45:  ## Universal likelihood function
    # index:      0   1    2    3   4   5    6   
    # parname:    s0  s1  labda p   q  phi1 phi2
    fpar_nuis = [0.1, 0, 0, 2, 1e10, 0, 0]
    parmin_nuis = [0, 0, -1, 0.5, 2, 0, 0]
    parmax_nuis = [1, 1, 1, 100, 100, 1, 1]
    LV['filename'] = 'UL'
    id_nuis = [0, 2, 3, 4, 5]

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
Meas_info = {'Y': daily_data[idx_t1:idx_t2, 5],  	    # Adjusted to match the MATLAB code (indexing)
             'sigma2': 'nonconstant'}

# Define nuisance variable names
nuis_names = nuis_var_names(DREAMPar)
names = par_names + nuis_names

# Model parameters and nuisance variables
n_modpar = len(fpar_mod)
LV['id_vpar'] = list(range(n_modpar)) + [x + n_modpar for x in id_nuis]
LV['fpar'] = fpar_mod + fpar_nuis
parmin = parmin_mod + parmin_nuis
parmax = parmax_mod + parmax_nuis

# Min/max values of parameter selection
Par_info['min'] = [parmin[i] for i in LV['id_vpar']]
Par_info['max'] = [parmax[i] for i in LV['id_vpar']]
Par_info['names'] = [names[i] for i in LV['id_vpar']]
DREAMPar['d'] = len(Par_info['min'])

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
    Par_info['steps'] = (250) * np.ones(DREAMPar['d'])

# Optional settings
options = {'modout': 'no',
           'parallel': 'no',
           'save': 'yes',
           'restart': 'no'}

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, LV, [], plugin)
