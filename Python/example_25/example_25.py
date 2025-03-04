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
## Example 25: Bottom-up control hypothesis: Bedrock depth from high-resolution       ##
##             topographic data and a geomorphologic model of landscape evolution     ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Gomes, G.J.C., J.A. Vrugt, and E.A. Vargas Jr., Towards improved prediction of   ##
##       the bedrock depth underneath hillslopes: Bayesian Inference of the bottom-up ##
##       control hypothesis using high-resolution topographic data, Water Resources   ##
##       Research, 52, 3085–3112, doi:10.1002/2015WR018147, 2016                      ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from DREAM_Suite_functions import nuis_var_names	            # Import functions

# Define DREAMPar dictionary
DREAMPar = {'lik': 13}                      # Set the likelihood type (Log-likelihood with AR(1) residuals)

# Initial sampling and parameter ranges
Par_info = {'initial': 'latin',   		    # Latin hypercube sampling
            'boundhandling': 'reflect'}  	# Explicit boundary handling: reflection

# Provide information parameter space and initial sampling
par_names = ['\\Phi','\\lambda_{1}','S_{\\rm c}','\\lambda_{2}']
fpar_mod   = [  0.001,   1.0,  1.0,  10 ];     # Default values parameters
parmin_mod = [  0.0001,  0.1,  0.8,  1  ];     # Lower bound parameters
parmax_mod = [  0.1,     3.0,  1.5,  20 ];     # Upper bound parameters

# Define name of function (.m file) for posterior exploration
Func_name = 'DTB_model'

# Global variable LV for likelihood functions
LV = {}

# Handling likelihood function selection with a switch-like structure
if DREAMPar['lik'] == 13:  # Normal distribution
    # index:      0   1    2   3
    # parname:   s0  s1  phi1 phi2    
    fpar_nuis = [0.1, 0, 0, 0]
    parmin_nuis = [0, 0, 0, 0]
    parmax_nuis = [5, 1, 1, 1]   
    LV['filename'] = 'Normal'
    id_nuis = [0, 1]
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
    id_nuis = [2]
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

# Load depth to bedrock calibration data
DTB_data = np.loadtxt('DTB_data.txt')

# Store training data information
Meas_info = {   'Y': DTB_data[:, 0],    # Measured discharge data
 	            'sigma2': 'constant'}   # Constant measurement error variance}
nuis_names = nuis_var_names(DREAMPar)   # Extract names of nuisance variables
print(par_names,nuis_names)
names = par_names + nuis_names          # Combine with parameter names
n_modpar = len(fpar_mod)                # Number of model parameters

# Model parameter + nuisance variable selection
LV['id_vpar'] = list(range(n_modpar)) + [x + n_modpar for x in id_nuis]
# LV['id_vpar'] = list(range(1, n_modpar + 1)) + [i + n_modpar for i in id_nuis]
LV['fpar'] = fpar_mod + fpar_nuis       # Merge default values of parameters and nuisance variables

# Merging the parameter names and nuisance variables
parmin = parmin_mod + parmin_nuis
parmax = parmax_mod + parmax_nuis

Par_info['min'] = [parmin[i] for i in LV['id_vpar']]    # Min values of selected parameters
Par_info['max'] = [parmax[i] for i in LV['id_vpar']]    # Max values of selected parameters
Par_info['names'] = [names[i] for i in LV['id_vpar']]   # Names of parameters and nuisance variables

DREAMPar['d'] = len(Par_info['min'])    # Number of parameters

# Optional settings
options = {'parallel': 'no',  	# Model is so fast -> parallel slows down
    	   'modout': 'yes',     # Store the simulated output of the model
    	   'save': 'yes'}     	# Save options enabled

# Define method to use {'dream', 'dream_zs', 'dream_d', 'dream_dzs', 'mtdream_zs'}
method = 'dream_d'

if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  		# Number of Markov chains
    DREAMPar['T'] = 10000  		# Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3  			# Number of Markov chains
    DREAMPar['T'] = 15000  		# Number of generations

if method == 'dream_d' or method == 'dream_dzs':
    Par_info['steps'] = (1000) * np.ones(DREAMPar['d'])          # Number of discrete steps

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, LV)
