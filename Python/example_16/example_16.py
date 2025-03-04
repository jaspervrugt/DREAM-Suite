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
## Example 16: Geophysical inversion                                                  ##
##                                                                                    ##
## The problem is simplified compared to the paper cited below as it considers a      ##
## problem in which the true model represent a model with the same dimension as the   ##
## inverse parameterization and it uses straight rays. Results can be visualized by   ##
## the function visualize_results.m                                                   ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation using Markov ##
##       chain Monte Carlo simulation: DREAM_(ABC), Water Resources Research,         ##
##       doi:10.1002/2014WR015386.                                                    ##
##   Linde, N., and J.A. Vrugt (2013), Spatially distributed soil moisture from       ##
##       traveltime observations of crosshole ground penetrating radar using Markov   ##
##       chain Monte Carlo, Vadose Zone Journal, 12(1), 1-16                          ##
##   Laloy, E., N. Linde, and J.A. Vrugt (2012), Mass conservative three-dimensional  ##
##       water tracer distribution from Markov chain Monte Carlo inversion of         ##
##       time-lapse ground-penetrating radar data, Water Resour. Res., 48, W07510,    ##
##       doi:10.1029/2011WR011238                                                     ##
##   Carbajal, M.R., N. Linde, T. Kalscheuer, and J.A. Vrugt (2014), Two-dimensional  ##
##       probabilistic inversion of plane-wave electromagnetic data: Methodology,     ##
##       model constraints and joint inversion with electrical resistivity data,      ##
##       Geophysical Journal International}, 196(3), 1508-1524,                       ##
##       doi: 10.1093/gji/ggt482                                                      ##
##   Lochbühler, T., S.J. Breen, R.L. Detwiler, J.A. Vrugt, and N. Linde (2014),      ##
##       Probabilistic electrical resistivity tomography for a CO_2 sequestration     ##
##       analog, Journal of Applied Geophysics, 107, 80-92,                           ##
##       doi:10.1016/j.jappgeo.2014.05.013                                            ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import sys, os

current_dir = os.getcwd()                                       # Get the current working directory
from DCT_GPR import GPR_par_ranges
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite


# DCT order in x and z dimension
parx = 8
parz = 8

# Problem settings defined by user
DREAMPar = {}
DREAMPar['d'] = parx * parz       # Dimension of the problem
DREAMPar['thinning'] = 5          # Only store each 5th sample
DREAMPar['lik'] = 2               # Model output is log-likelihood

# Define name of function (.py file) for posterior exploration
Func_name = 'DCT_GPR.DCT_GPR'

# Provide information on parameter space and initial sampling
func = {'parx': parx, 'parz': parz}

# Define parameter ranges
Par_info, slowtrue = GPR_par_ranges(func)

# Initial sampling distribution
Par_info['initial'] = 'normal'                                                                      # N(µ, Σ) initial sampling distribution
Par_info['mu'] = Par_info['min'] + 0.5 * (np.array(Par_info['max']) - np.array(Par_info['min']))    # Mean µ of distribution
Par_info['cov'] = 0.001 * np.diag(np.array(Par_info['max']) - np.array(Par_info['min']))            # Covariance matrix Σ of distribution

# Define method to use (choices: 'dream', 'dream_zs', 'dream_d', 'dream_dzs', 'mtdream_zs')
method = 'dream_zs'

# Set method parameters
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 32                # Number of Markov chains
    DREAMPar['T'] = 10000             # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 5                 # Number of Markov chains
    DREAMPar['T'] = 60000             # Number of generations

# If method requires discrete steps, define the step size
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (1000) * np.ones(DREAMPar['d'])

# Call the DREAM-Suite package
if __name__ == '__main__':
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info)
