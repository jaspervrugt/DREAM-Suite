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
## Example 4: Distribution-adaptive likelihood functions applied to conceptual        ##
##            watershed models using measured discharge data                          ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
##       distribution-adaptive likelihood functions: Generalized and universal        ##
##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
##       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             ##
##       https://www.sciencedirect.com/science/article/pii/S002216942201112X          ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Schoups, G., J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption   ##
##       of accuracy and efficiency of Markov Chain Monte Carlo simulation by         ##
##       inaccurate numerical implementation of conceptual hydrologic models, Water   ##
##       Resources Research, 46, W10530, doi:10.1029/2009WR008648                     ##
##   Schoups, G., and J.A. Vrugt (2010), A formal likelihood function for parameter   ##
##       and predictive inference of hydrologic models with correlated,               ##
##       heteroscedastic and non-Gaussian errors, Water Resources Research, 46,       ##
##       W10531, doi:10.1029/2009WR008933                                             ##
##   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson (2009),             ##
##       Equifinality of formal (DREAM) and informal (GLUE) Bayesian approaches in    ##
##       hydrologic modeling?, Stochastic Environmental Research and Risk Assessment, ##
##       23(7), 1011-1026, doi:10.1007/s00477-008-0274-y                              ##
##   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and        ##
##       J.M. Hyman (2009), Accelerating Markov chain Monte Carlo simulation by       ##
##       differential evolution with self-adaptive randomized subspace sampling,      ##
##       International Journal of Nonlinear Sciences and Numerical Simulation, 10(3), ##
##       271-288                                                                      ##
##   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), ##
##       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       ##
##       backward with Markov chain Monte Carlo simulation, Water Resources Research, ##
##       44, W00B09, doi:10.1029/2007WR006720                                         ##
##   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian (2003), A Shuffled Complex  ##
##       Evolution Metropolis algorithm for optimization and uncertainty assessment   ##
##       of hydrologic model parameters, Water Resour. Res., 39 (8), 1201,            ##
##       doi:10.1029/2002WR001642.                                                    ##
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
DREAMPar = LV = {}
DREAMPar['lik'] = 44;   # 13: Normal likelihood function with AR(2)-model of residuals
                        # 14: Generalized likelihood function (Obsolete)
                        # 16: Laplacian likelihood function with AR(2)-model of residuals
                        # 17: Skewed Student t likelihood function
                        # 44: Generalized likelihood function PLUS
                        # 45: Universal likelihood function
## Please refer to Vrugt et al. (2022) for theory on likelihood functions

# Function name for posterior exploration
Func_name = 'hymod.hymod'

# Parameter names for model
par_names = ['S_{\\rm u,max}', '\\beta', '\\alpha', 'K_{\\rm s}', 'K_{\\rm f}']

# Model parameters
fpar_mod = [200, 1.00, 0.5, 0.01, 0.50]
parmin_mod = [50, 0.10, 1e-4, 1e-4, 0.25]
parmax_mod = [1000, 10.0, 1.0, 0.25, 5.00]

# Parameter space and initial sampling information
Par_info = {}
Par_info['initial'] = 'latin'  		    	# Latin hypercube sampling
Par_info['boundhandling'] = 'reflect'  		# Boundary handling
Par_info['norm'] = 0  		                # Work in normalized space, if needed

# Handling likelihood function selection with a switch-like structure
if DREAMPar['lik'] == 13:                   # Normal distribution
    # index:      0   1    2   3
    # parname:   s0  s1  phi1 phi2    
    fpar_nuis = [0.1, 0, 0, 0]
    parmin_nuis = [0, 0, 0, 0]
    parmax_nuis = [1, 1, 1, 1]
    LV['filename'] = 'Normal'
    id_nuis = [2]
elif DREAMPar['lik'] == 14:  # GL (OBSOLETE)
    # index:      0   1    2    3   4   5   6    7    8    9  10  
    # parname:  std0 std1 beta xi  mu1 phi1 phi2 phi3 phi4 K lambda
    fpar_nuis = [0.1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
    parmin_nuis = [0, 0, -1, 0.1, 0, 0, 0, 0, 0, 0, 0.1]
    parmax_nuis = [1, 1, 1, 10, 100, 1, 1, 1, 1, 1, 1]
    LV['filename'] = 'GL'
    id_nuis = [0, 1, 2, 3, 5]
elif DREAMPar['lik'] == 16:  # Laplace distribution
    # index:      0   1   2
    # parname:   s0  s1  phi1
    fpar_nuis = [0.1, 0, 0]
    parmin_nuis = [0, 0, 0]
    parmax_nuis = [1, 1, 1]
    LV['filename'] = 'Laplace'
    id_nuis = [0]
elif DREAMPar['lik'] == 17:  # SST (Skewed Student-t)
    # index:      0   1   2    3   4    5
    # parname:   s0  s1  nu   xi phi1  phi2
    fpar_nuis = [0.1, 0, 1e10, 1, 0, 0]
    parmin_nuis = [0, 0, 2, 0.1, 0, 0]
    parmax_nuis = [1, 1, 100, 10, 1, 1]
    LV['filename'] = 'SL'
    id_nuis = [2, 3, 4]
elif DREAMPar['lik'] == 44:  # GL+
    # index:      0   1   2    3   4    5
    # parname:   s0  s1  beta xi phi1  phi2
    fpar_nuis = [0.1, 0, 0, 1, 0, 0]
    parmin_nuis = [0, 0, -1, 0.1, 0, 0]
    parmax_nuis = [1, 1, 1, 10, 1, 1]
    LV['filename'] = 'GL_plus'
    id_nuis = [2, 3, 4]
elif DREAMPar['lik'] == 45:  # Universal likelihood function
    # index:      0   1    2    3   4   5    6   
    # parname:    s0  s1  labda p   q  phi1 phi2
    fpar_nuis = [0.1, 0, 0, 2, 1e10, 0, 0]
    parmin_nuis = [0, 0, -1, 0.5, 2, 0, 0]
    parmax_nuis = [1, 1, 1, 100, 100, 1, 1]
    LV['filename'] = 'UL'
    id_nuis = [0, 2, 3, 4, 5]

# Loading data (here we assume `bound.txt` is a file)
bound = np.loadtxt('bound.txt')

# Define warm-up period and indices for measured discharge record
plugin = {}
plugin['data'] = {}
plugin['data']['wmp'] = 64
# idx_t1 = 65   # IN MATLAB
idx_t1 = 64
idx_t2 = bound.shape[0]
idx_p = np.arange(idx_t1 - plugin['data']['wmp'], idx_t2)

# The first two years are warm-up
plugin = {'idx': np.arange(plugin['data']['wmp'], len(idx_p) + 1),        
          'data': {'P': np.sum(bound[idx_p, 5:9], axis = 1),    # Rainfall (mm/d)
                  'Ep': bound[idx_p, 4]},           # Potential evapotranspiration (mm/d)
          'y0': 1e-5 * np.ones(6),                  # Initial states
          'options': {'InitialStep': 1,             # Initial time step hmodel
                      'MaxStep': 1,                 # Maximum time step hmodel
                      'MinStep': 1e-3,              # Minimum time step hmodel
                      'RelTol': 1e-3,               # Relative tolerance 
                      'AbsTol': 1e-3,               # Absolute tolerance
                      'Order': 2},                  # Order of numerical solution
          'tout': np.arange(0, len(idx_p) + 1)}

# Measured streamflow data (converted from m3/s to mm/d)
F = 1944 * (1000 * 1000) / (1000 * 60 * 60 * 24)
Meas_info = {}
Meas_info['Y'] = 1 / F * bound[idx_t1:idx_t2, 3]
Meas_info['sigma2'] = 'nonconstant'

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

# Define method for Markov chains
method = 'dream_kzs'

if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  	# Number of Markov chains
    DREAMPar['T'] = 2500  	# Generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 5
    DREAMPar['T'] = 15000
elif method == 'dream_kzs':
    DREAMPar['N'] = 5
    DREAMPar['T'] = 5000
    DREAMPar['M'] = 24  	# Archive samples Kalman jump

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (250) * DREAMPar['d']

# Optional settings
options = {}
options['modout'] = 'yes'  	    # Return model simulation samples?
options['parallel'] = 'yes'     # Run each chain on a different core
options['save'] = 'yes'  	    # Save workspace during run

if __name__ == '__main__':
    # Call the DREAM-Suite function
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, LV, [], plugin)
