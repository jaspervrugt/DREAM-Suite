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
## Example 21: Soil water flow modeling: Limits of Acceptability                      ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. and K.J. Beven (2018), Embracing equifinality with efficiency:       ##
##       Limits of Acceptability sampling using the DREAM_{(LOA)} algorithm,  Journal ##
##       of Hydrology,  559 , pp. 954-971, doi:10.1016/j.jhydrol.2018.02.026          ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Scharnagl, B., J.A. Vrugt, H. Vereecken, and M. Herbst (2011), Bayesian inverse  ##
##	     modeling of soil water dynamics at the field scale: using prior information  ##
##	     on soil hydraulic properties, Hydrology and Earth System Sciences, 15,       ##
##       3043–3059, doi:10.5194/hess-15-3043-2011                                     ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys
import pandas as pd
from scipy.stats import norm, uniform
from scipy import io

current_dir = os.getcwd()                                       # Get the current working directory
from HYDRUS import Load_data

parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

# Problem settings defined by user
DREAMPar = {'d': 7,                                 # Dimension of the problem
            'lik': 23}                              # Limits of acceptability likelihood

# Provide information on parameter space and initial sampling
Par_info = {'initial': 'latin',            			# Initial sample from Latin hypercube sampling
            'boundhandling': 'reflect',    			# Explicit boundary handling
            'names': ['\\theta_{\\rm r}$', '\\theta_{\rm s}', '\\alpha', 'n', 
                      'K_{\\rm s}', 'l', 'h_{\\rm bot}'],  		                # Parameter names
            'min': [0.0430, 0.4090, -2.5528, 0.1790, -2.2366, -5.49, -250],     # Lower bounds
            'max': [0.0910, 0.4810, -2.0706, 0.2670, -0.0800,  6.27, -50]}      # Upper bounds
# l is more or less fixed this way so we can use setup of case 7 

# Define name of function for posterior exploration
Func_name = 'HYDRUS.HYDRUS'

# Load the experimental data (assuming it's an Excel file)
data = pd.read_excel('all_observations.xlsx', sheet_name = 'TDR_data', usecols = "A, D:AG", header = 2).to_numpy()  # Adjust columns as needed

# Focus on the right columns
data = data[:, [1] + list(range(3, data.shape[1]))]

# Mean soil moisture at time t as observed summary metrics
Meas_info = {'S': np.nanmean(data, axis=0)}  # Mean values for each column, ignoring NaNs

# Determine limits of acceptability: epsilon for each entry in Meas_info.S
D = np.sort(data, axis = 0)		    # Sort each data column ascending order
m, n = D.shape  			        # Number of rows (observations) and columns (measurements)
alfa1, alfa2 = 0.025, 1 - 0.025  	# Confidence levels for significance level alfa = 0.05

# 95% Soil Moisture range at each measurement time (percentile calculation)
rng = np.nanpercentile(D, [100 * alfa1, 100 * alfa2], axis = 0)  	# 95% range
options = {'epsilon': 0.5 * np.diff(rng, axis = 0).flatten()}  	    # Limits of acceptability (LOA)

# Load other data that needs to be ported to HYDRUS
data_hydrus = Load_data()  		    # Assuming Load_data function is available

# Optional settings
options['IO'] = 'yes'      		    # Input-output writing of model files (for parallel processing)
options['parallel'] = 'yes'  	    # Run chains in parallel
options['save'] = 'yes'    		    # Save memory of DREAM during trial
options['modout'] = 'yes'  		    # Save model output

# Define method to use {'dream', 'dream_zs', 'dream_d', 'dream_dzs', 'mtdream_zs'}
method = 'dream_zs'

# Define the number of chains and generations based on the method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10   		    # Number of Markov chains
    DREAMPar['T'] = 2500 		    # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3    		    # Number of Markov chains
    DREAMPar['T'] = 8000 		    # Number of generations

# If the method involves discrete steps, define them
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = (500) * DREAMPar['d']  # Discrete steps for each parameter

if __name__ == '__main__':
    # Call the DREAM-Suite package
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, [], [], data_hydrus)
