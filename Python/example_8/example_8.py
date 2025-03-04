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
## Example 8: Approximate Bayesian Computation: Benchmark function                    ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation using Markov ##
##       chain Monte Carlo simulation: DREAM_(ABC), Water Resources Research,         ##
##       doi:10.1002/2014WR015386.                                                    ##
##   Sadegh, M., and J.A. Vrugt (2013), Bridging the gap between GLUE and formal      ##
##       statistical approaches: approximate Bayesian computation, Hydrology and      ##
##       Earth System Sciences, 17, 4831–4850                                         ##
##   Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model calibration and       ##
##       evaluation: Approximate Bayesian computation, Water Resources Research, 49,  ##
##       4335–4345, doi:10.1002/wrcr.20354                                            ##
##   Turner, B.M., and P.B. Sederberg (2013), Approximate Bayesian computation with   ##
##       differential evolution, Journal of Mathematical Psychology, In Press.        ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys

current_dir = os.getcwd()                                       # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

# Problem settings defined by user
DREAMPar = {}
DREAMPar['d'] = 1                         # Dimension of the problem
DREAMPar['lik'] = 22                      # ABC informal likelihood function
DREAMPar['delta'] = 1                     # Use only 1 pair of chains to create proposal

# Provide information about the parameter space and initial sampling
Par_info = {}
Par_info['initial'] = 'latin'             # Latin hypercube sampling
Par_info['boundhandling'] = 'fold'        # Explicit boundary handling
Par_info['min'] = [-10]                   # If 'latin', min values
Par_info['max'] = [10]                    # If 'latin', max values
Par_info['norm'] = 1                      # Sample in normalized space  

# Define name of function for posterior exploration
Func_name = 'ABC_func.ABC_func'

# Define Meas_info.S
Meas_info = {}
Meas_info['S'] = [0]

# Optional settings
options = {}
options['epsilon'] = 0.025               # Epsilon of the noisy ABC implementation
options['rho'] = 'lambda X, Y: X - Y'    # Define the distance function (Python equivalent of inline)

# Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs'

# Setting parameters based on the method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10                        # Markov chains
    DREAMPar['T'] = 10000                     # generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3                         # Markov chains
    DREAMPar['T'] = 25000                     # generations

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = 1000                  # discrete steps

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options)
