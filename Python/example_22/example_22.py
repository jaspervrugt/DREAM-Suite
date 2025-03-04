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
## Example 22: Nash-Cascade model: Limits of Acceptability                            ##
##                                                                                    ##
## Check the following papers                                                         ##
##   Vrugt, J.A. and K.J. Beven (2018), Embracing equifinality with efficiency:       ##
##       Limits of Acceptability sampling using the DREAM_{(LOA)} algorithm,  Journal ##
##       of Hydrology,  559 , pp. 954-971, doi:10.1016/j.jhydrol.2018.02.026          ##
##   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software ##
##       package: Theory, concepts, and MATLAB implementation, Environmental Modeling ##
##       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             ##
##   Nash, J.E. (1960), A unit hydrograph study with particular reference to British  ##
##       catchments, Proceedings - Institution of Civil Engineers, 17, 249-282        ##
##   Nash, J.E., J.V. Sutcliffe (1970), River flow forecasting through conceptual     ##
##       models part I - A discussion of principles, Journal of Hydrology, 10(3),     ##
##       282-290                                                                      ##
##   Beven, K.J., and A.M. Binley (1992), The future of distributed models: Model     ##
##       calibration and uncertainty prediction. Hydrological Processes, 6, 279–298,  ##
##       doi: 10.1002/hyp.3360060305                                                  ##
##                                                                                    ##
## ################################################################################## ##

import numpy as np
import os, sys
from scipy.stats import norm

current_dir = os.getcwd()                                       # Get the current working directory
from Nash_Cascade_2 import Nash_Cascade_2

parent_dir = os.path.abspath(os.path.join(current_dir, '..'))   # Go up one directory
sys.path.append(parent_dir)                                     # add this to path
from DREAM_Suite import DREAM_Suite

# Define problem settings
DREAMPar = {'d': 2,     # Dimension of the problem
            'lik': 23   # Limits of acceptability
}
# Parameter space and sampling
Par_info = {'initial': 'latin',             # Latin hypercube sampling
            'boundhandling': 'reflect',     # Boundary handling
            'names': ['k', 'n'],            # Parameter names
	        'min': [1, 1],                  # Min values for Latin hypercube
            'max': [10, 10],                # Max values for Latin hypercube
}

# Function name for posterior exploration
Func_name = 'Nash_Cascade_2.Nash_Cascade_2'

# Create synthetic time series data
y = Nash_Cascade_2([4, 2])

# Define heteroscedastic measurement error
Meas_info = {'S': norm.rvs(loc = y, scale = 0.1*y)}

# Optional settings
options = {'modout': 'yes',  			    # Return model simulations
           'epsilon': [0.2*Meas_info['S']]}	# Epsilon values

# Define method and MCMC parameters
method = 'dream_dzs'

# Set MCMC parameters based on the method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10              # Number of Markov chains
    DREAMPar['T'] = 2000            # Number of generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3               # Number of Markov chains
    DREAMPar['T'] = 7000            # Number of generations

# If using discrete steps method
if method in ['dream_d', 'dream_dzs']:
    Par_info['steps'] = 900 * np.ones(DREAMPar['d'])

# Call the DREAM-Suite package
if __name__ == '__main__':
    chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options)
