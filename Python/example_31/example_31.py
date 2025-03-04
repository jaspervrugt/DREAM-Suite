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
## Example 31: AR(2)-model parameter estimation: Test of likelihood functions         ##
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
##                                                                                    ##
## ################################################################################## ##


import numpy as np
import os, sys
import scipy.io

current_dir = os.getcwd()                                               # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))           # Go up one directory
sys.path.append(parent_dir)                                             # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	            # Add miscellaneous directory to Python path
from DREAM_Suite_functions import nuis_var_names, SEPrnd, SGTrnd     # Import functions

# Provide information for parameter space and initial sampling
Par_info = {'initial': 'latin',  		    # Latin hypercube sampling
	 	    'boundhandling': 'reflect'}  	# Explicit boundary handling

par_names, fpar_mod, parmin_mod, parmax_mod = [], [], [], []

# Define name of function for posterior exploration
Func_name = 'ar2_model.ar2_model'

# Determine likelihood function
DREAMPar = {'lik': 44}

# Likelihood function - summary variable
LV = {}

# Handling likelihood function selection with a switch-like structure
if DREAMPar['lik'] == 44:  # GL+
    # index:      0   1   2    3   4    5
    # parname:   s0  s1  beta xi phi1  phi2
    fpar_nuis = [1.0, 0, 0, 1, 0, 0]
    parmin_nuis = [0, 0, -1, 0.1, 0, 0]
    parmax_nuis = [2, 1, 1, 10, 1, 1]
    LV['filename'] = 'GL_plus'
    id_nuis = [2, 3, 4, 5]
elif DREAMPar['lik'] == 45:  # Universal likelihood function
    # index:      0   1    2    3   4   5    6   
    # parname:    s0  s1  labda p   q  phi1 phi2
    fpar_nuis = [1.0, 0, 0, 2, 1e10, 0, 0]
    parmin_nuis = [0, 0, -1, 0.5, 2, 0, 0]
    parmax_nuis = [2, 1, 1, 100, 100, 1, 1]
    LV['filename'] = 'UL'
    id_nuis = [2, 3, 4, 5, 6]

# Load data
file_name = 'y.mat'

# Create autocorrelated data
if not os.path.isfile(file_name):
    M = 5000
    phi1 = 0.7
    phi2 = 0.2
    y = np.zeros((M, 2))
    
    for i in range(2):
        if i == 0:      # 'SEP' innovations
            beta = -0.5
            xi = 5
            rnd = SEPrnd(beta, xi, M, 1)
        elif i == 1:    # 'SGT' innovations
            lambda_ = -0.5
            p = 1.2
            q = 5.7
            rnd = SGTrnd(0, 1, lambda_, p, q, M, 1)
        
        for t in range(2, M):
            y[t, i] = phi1 * y[t-1, i] + phi2 * y[t-2, i] + rnd[t, 0]
    
    scipy.io.savemat(file_name, {'y': y, 'M': M})
else:
    data = scipy.io.loadmat(file_name)
    y = data['y']
    M = data['M']

plugin = {'N': int(M)}

if DREAMPar['lik'] == 44:
    Meas_info = {'Y': y[:, 0]}
elif DREAMPar['lik'] == 45:
    Meas_info = {'Y': y[:, 1]}

Meas_info['sigma2'] = 'constant'  # Constant measurement error variance

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

# Optional DREAM settings
options = {'parallel': 'no',    # This example runs in parallel
           'save': 'yes',       # Save DREAM workspace during run
           'modout': 'no'}      # Store model simulations

# Define method to use
method = 'dream_zs'
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10          # Number of Markov chains
    DREAMPar['T'] = 3000        # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3           # Number of Markov chains
    DREAMPar['T'] = 10000       # Number of generations

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, LV, [], plugin)
