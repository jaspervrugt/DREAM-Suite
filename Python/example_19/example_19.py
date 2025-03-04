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
## Example 19: BMA mixture model estimation. I recommend using the MODELAVG toolbox!  ##
##             This provides a much more comprehensive treatment of the BMA model,    ##
##             including confidence/prediction intervals, scoring rules and other     ##
##             metrics of the BMA distribution forecast                               ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and Diagnostics:         ##
##       Elicitability, Propriety, and Scoring Rules for Hydrograph Functionals,      ##
##       Water Resources Research, 60, e2023WR036710,                                 ##
##       https://doi.org/10.1029/2023WR036710                                         ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Vrugt, J.A., C.G.H. Diks, and M.P. Clark (2008), Ensemble Bayesian model         ##
##       averaging using Markov chain Monte Carlo sampling, Environmental Fluid       ##
##       Mechanics, 8(5-6), 579-595, doi:10.1007/s10652-008-9106-3                    ##
##   Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty using ensemble   ##
##       methods: Comparison of sequential data assimilation and Bayesian model       ##
##       averaging, Water Resources Research, 43, W01411, doi:10.1029/2005WR004838    ##
##                                                                                    ##
## ################################################################################## ##


import numpy as np
import sys, os

current_dir = os.getcwd()                                       # Get the current working directory
from BMA_lik import setup_BMA
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	    # Add miscellaneous directory to Python path
from DREAM_Suite_functions import genparset	                    # Import functions

# Define problem settings
DREAMPar = {'lik': 2}  # Return argument is log-likelihood of BMA model

# Parameter space and sampling
Par_info = {'initial': 'latin',             # Latin hypercube sampling
            'boundhandling': 'reflect'}     # Boundary handling

## Define name of function
Func_name = 'BMA_lik.BMA_lik'

# Load the data (assuming 'discharge.txt' is a CSV file or similar)
data = np.loadtxt('discharge.txt')

# Define the training period and select the corresponding data
T_idx = np.arange(0, 3000)  # Indices for the training period (0-based)
D = data[T_idx, :8]         # Ensemble forecasts (assuming 8 models in columns)
y = data[T_idx, 8]          # Verifying observations (9th column in MATLAB, 8th in Python)

# Set up options for BMA method
options = {}
options['BMA'] = 'yes'      # Activate BMA method
PDF = 'normal'              # Forecast pdf ('normal'/'gamma')
VAR = '4'                   # Variance option ('1'/'2'/'3'/'4')

# Call the setup_BMA function to prepare the BMA model with bias correction
DREAMPar, Par_info, D_bc, A, B = setup_BMA(DREAMPar, Par_info, D, y, VAR)

# Structure for BMA info
plugin = {}
plugin['BMA'] = {'PDF': PDF,        # Choice of BMA pdf
                 'VAR': VAR,        # Choice of variance of BMA pdf
                 'D': D_bc,         # Bias-corrected forecasts of ensemble members: nxK matrix
                 'y': y,            # Measured discharge in mm/d: nx1 vector
                 'K': D.shape[1]}   # Number of ensemble members (columns)

# Set options for the optimization process
options['print'] = 'yes'    # Print output to screen (figures)
options['modout'] = 'yes'   # Return hmodel simulations
options['save'] = 'yes'     # Save memory during AMALGAM restart run
options['parallel'] = 'no'  # Do not run in parallel

## Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs'

if method in ['dream','dream_d']:
    DREAMPar['N'] = 10;             # Markov chains
    DREAMPar['T'] = 5000;           # generations
elif method in ['dream_zs','dream_dzs','mtdream_zs']:
    DREAMPar['N'] = 3;              # Markov chains
    DREAMPar['T'] = 25000;          # generations

if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (1000) * DREAMPar['d']

if __name__ == '__main__':

    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, [], options, [], [], plugin)

    # Create matrix of 3d chain trajectories
    parset = genparset(chain)
    # Apply burn in
    P = parset[-20000:, :DREAMPar['d']]
    # Weight normalization 
    P[:, :8] = P[:, :8] / np.sum(P[:, :8], axis=1, keepdims=True)
    # Print to screen
    print(P)
    ## Note: the MODElAVG toolbox is specifically developed for ensemble postprocessing using BMA methods
