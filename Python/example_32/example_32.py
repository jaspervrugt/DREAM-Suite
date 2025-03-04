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
## Example 32: Distribution-adaptive likelihood functions applied to conceptual       ##
##             watershed models using measured discharge data                         ##
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
##   Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten (2005),     ##
##       Improved treatment of uncertainty in hydrologic modeling: Combining the      ##
##       strengths of global optimization and data assimilation, Water Resources      ##
##       Research, 41, W01017, doi:10.1029/2004WR003059                               ##
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

current_dir = os.getcwd()                                               # Get the current working directory
from hymod import load_data
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))           # Go up one directory
sys.path.append(parent_dir)                                             # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	            # Add miscellaneous directory to Python path
from DREAM_Suite_functions import define_lik, nuis_var_names         # Import functions

## Problem settings defined by user
DREAMPar = LV = {}
DREAMPar['lik'] = 17;   # 13: Normal likelihood function with AR(2)-model of residuals
                        # 14: Generalized likelihood function (Obsolete)
                        # 16: Laplacian likelihood function with AR(2)-model of residuals
                        # 17: Skewed Student t likelihood function
                        # 44: Generalized likelihood function PLUS
                        # 45: Universal likelihood function
## Please refer to Vrugt et al. (2022) for theory on likelihood functions

# Function name for posterior exploration
Func_name = 'hymod.hymod'

# Load the data
Meas_info, plugin = load_data()

# Model parameters
fpar_mod = [200, 1.00, 0.5, 0.01, 0.50]
parmin_mod = [50, 0.10, 0.0, 1e-5, 0.1]
parmax_mod = [1000, 10.0, 1.0, 0.25, 5.00]
par_names = ['S_{\\rm u,max}', '\\beta', '\\alpha', 'K_{\\rm s}', 'K_{\\rm f}']

# Parameter space and initial sampling information
Par_info = {'initial': 'latin',  		    	# Latin hypercube sampling
            'boundhandling': 'reflect',  		# Boundary handling
            'norm': 1}  		                # Work in normalized space, if needed

# Define likelihood properties: [user needs to select idx_nuis]
DREAMPar, Par_info, id_vpar, fpar, fname = define_lik(DREAMPar, Par_info, par_names, fpar_mod, parmin_mod, parmax_mod)

# Define likelihood variables
LV = {  'id_vpar': id_vpar,     # variable parameters 
        'fpar': fpar,           # fixed parameter values
        'filename': fname}      # name of likelihood function

# How many estimeable parameters?
DREAMPar['d'] = len(Par_info['min'])

# Optional settings
options = { 'modout': 'yes',  	    # Return model simulation samples?
            'parallel': 'no',    	# Run each chain on a different core
            'save': 'yes'}  	    # Save workspace during run

# Define method for Markov chains
method = 'dream_zs'

if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  	# Number of Markov chains
    DREAMPar['T'] = 5000  	# Generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 4
    DREAMPar['T'] = 30000

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (1000) * DREAMPar['d']

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, LV, [], plugin)
