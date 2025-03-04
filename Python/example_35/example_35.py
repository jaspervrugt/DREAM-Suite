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
## Example 35: Inverse modeling of Parlange's semi-analytic infiltration equation     ##
##             using the SWIG database of measured infiltration curves                ##
##                                                                                    ##
## Check the following papers                                                         ##
##  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and               ##
##      H. Vereecken (2023), The time validity of Philip's two-term infiltration      ##
##      equation: An elusive theoretical quantity? Vadose Zone Journal, e20309,       ##
##      pp. 1-25, https://doi.org/10.1002/vzj2.20309                                  ##
##  Vrugt, J.A. and Y. Gao (2022), On the three-parameter infiltration equation of    ##
##      Parlange et al. (1982): Numerical solution, experimental design, and          ##
##      parameter estimation, Vadose Zone Journal, 21:e20167, pp. 1-25,               ##
##      https://doi.org/10.1002/vzj2.20167                                            ##
##                                                                                    ##
## ---------------------------------------------------------------------------------- ##

import numpy as np
import os, sys
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
import scipy.io

current_dir = os.getcwd()                                               # Get the current working directory
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))           # Go up one directory
sys.path.append(parent_dir)                                             # add this to path
from DREAM_Suite import DREAM_Suite

sys.path.append(os.path.join(parent_dir, 'miscellaneous'))	            # Add miscellaneous directory to Python path
from DREAM_Suite_functions import X_unnormalize, genparset              # Import functions

# Set up plugin and model configuration
plugin = {'model_setup': 1}

# Provide information for parameter space and initial sampling
Par_info = {'initial': 'normal',  		    # Latin hypercube sampling
	 	    'boundhandling': 'fold',  	    # Explicit boundary handling
		    'mu': [10, 10, 1],              # µ of N(µ,Σ)
		    'cov': np.array([[2, 0, 0],
                   [0, 2, 0],
                   [0, 0, 0.3]])}           # Σ of N(µ,Σ)

## Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs'

# Set parameter bounds based on model setup
if plugin['model_setup'] == 1:
    Par_info['names'] = ['S','K_{\\rm s}','\\beta']
    Par_info['min'] = [0, 0, 0]  	    	# Minimum parameter values
    Par_info['max'] = [1e5, 1e5, 2]  	    # Maximum parameter values
elif plugin['model_setup'] == 2:
    Par_info['names'] = ['S','K_{\\rm s}','K_{\\ rm i}','\\beta']
    Par_info['min'] = [0, 0, 0, 0]  		# Minimum parameter values
    Par_info['max'] = [1e5, 1e5, 1e2, 2]  	# Maximum parameter values

# Dimension of the problem
DREAMPar = {'d': len(Par_info['min'])} 

# Set number of chains and generations based on the method
if method in ['dream', 'dream_d']:
    DREAMPar['N'] = 10  		            # Markov chains
    DREAMPar['T'] = 10000  		            # generations
elif method in ['dream_zs', 'dream_dzs', 'mtdream_zs']:
    DREAMPar['N'] = 3  			            # Markov chains
    DREAMPar['T'] = 15000  		            # generations

# Optional settings
options = {	'modout': 'yes',  	# Return model simulations (yes/no)?
    		'parallel': 'no', 	# Run each chain on a different core
    		'save': 'yes',  	# Save workspace during run
    		'print': 'no'} 		# No figures printed to screen

# Load data (from .mat file)
data = scipy.io.loadmat('SWIG_1D.mat')
data_SWIG = data['SWIG_1D']
n_soil = data_SWIG.shape[0]

method_to_use = 3
# [1] synthetic data: Paper 1                           [= example_34]
# [2] SWIG data - time & infiltration form: Paper 2     [= example 35]
# [3] SWIG data - time & infiltration form combined     [= example 35]

if __name__ == '__main__':
    # Loop through each soil type

    # Define the function for the implementation
    if method_to_use == 2:
        DREAMPar['lik'] = 11                                                    # Gaussian likelihood
        p95_SWIG = np.full((n_soil, 5, DREAMPar['d'], 2), np.nan)  			    # Initialize parameter ranges
        ML_SWIG = np.full((n_soil, DREAMPar['d'], 2), np.nan)  				    # Initialize max likelihood values
        alg_SWIG = P_SWIG = FX95_SWIG = None                                    # Initialize results
        conv_DREAM = MRstat_DREAM = np.full((n_soil, 2), np.nan)  			    # Initialize convergence diagnostic

        current_dir = os.getcwd()

        for app in range(2):  								                    # Loop over each approach
            for st in range(n_soil):  							                # Loop over each SWIG soil
                dat = data_SWIG[st]  							                # Unpack data
                dat = dat[0]                                                    # Python added by JAV
                n_data = dat.shape[0]  				        	                # Define plugin structure
                # if FX95_SWIG is None:
                #     FX95_SWIG = np.full((n_data, 4, n_soil, 2), np.nan)       # Initialize results

                if app == 0:  								                    # Minimize cumulative infiltration residuals
                    plugin['t'] = dat[:, 0]
                    Func_name = 'Haverkamp_I.Haverkamp_I'
                    Meas_info = {'Y': dat[:, 1]}  				                # Measurement data (in cm)
                else:                                                           # Minimize time residuals
                    plugin['I'] = dat[:, 1]
                    Func_name = 'Haverkamp_t.Haverkamp_t'
                    Meas_info = {'Y': dat[:, 0]}  				                # Measurement data (in hours)

                # Run DREAM-Suite package
                chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, [], [], plugin)
                if alg_SWIG is None:
                    alg_SWIG = np.full((output['MR_stat'].shape[0], 2, n_soil, 2), np.nan)
                else:
                    alg_SWIG[:, :, st, app] = output['MR_stat']                        

                parset = genparset(chain)  						                # Extract sample chains
                N = parset.shape[0]  							                # Number of samples
                P = parset[int(2 / 3 * N):, :DREAMPar['d'] + 2]                 # Burn-in posterior samples
                M = P.shape[0]  							                    # Number of posterior samples
                if P_SWIG is None:
                    P_SWIG = np.full((M, DREAMPar['d'] + 2, n_soil, 2), np.nan)
                else:
                    P_SWIG[:, :, st, app] = P                        

                # Index ML parameters
                ii = np.argmax(np.sum(P[:M, DREAMPar['d']:DREAMPar['d']+2], axis = 1))
                ML_SWIG[st, :, app] = P[ii, :DREAMPar['d']]  		            # Store ML parameters

                FX = FX[int(2 / 3 * N):, :n_data]  		                        # Get posterior simulations
                A = np.sort(FX, axis = 0)  						                # Sort simulations
                n1 = int(0.025 * M)  							                # 2.5% prediction limit
                n2 = int(0.5 * M)  							                    # Mean of prediction interval
                n3 = int(0.975 * M)  							                # 97.5% prediction limit

                fx95 = np.zeros((n_data, 4))  				                    # Initialize fx95 array
                fx95[:, 0] = A[n1, :]  							                # Store 2.5%
                fx95[:, 1] = A[n2, :]  							                # Store mean
                fx95[:, 2] = FX[ii, :]  						                # Store ML simulation
                fx95[:, 3] = A[n3, :]  							                # Store 97.5% prediction limit
                # FX95_SWIG[:, :, st, app] = fx95  						        # Store for soil

                for j in range(DREAMPar['d']):  				                # Do the same for parameters
                    a = np.sort(P[:, j])
                    p95_SWIG[st, :5, j, app] = [a[n1], np.mean(a), np.median(a), a[n3], np.std(a)]

                conv = 1 if output['MR_stat'][-1, 1] < 1.2 else 0  	            # Chains converged or not
                conv_DREAM[st, app] = conv  						            # Store convergence
                MRstat_DREAM[st, app] = output['MR_stat'][-1, 1]  	            # Store R_hat diagnostic

                # Save results
                mat_file = f"mat_files/{st+1}_{app+1}.mat"
                np.savez(mat_file, st=st, ML_SWIG=ML_SWIG, P=P, fx95=fx95, conv_DREAM=conv_DREAM, MRstat_DREAM=MRstat_DREAM)

        # Save optimal values
        np.savez("DREAM_SWIG.mat", ML_SWIG=ML_SWIG, p95_SWIG=p95_SWIG, P_SWIG=P_SWIG, data_SWIG=data_SWIG, alg_SWIG=alg_SWIG, n_soil=n_soil, MRstat_DREAM=MRstat_DREAM, conv_DREAM=conv_DREAM)

    elif method_to_use == 3:
        Func_name = 'Haverkamp_It.Haverkamp_It' 				                # Specify function name
        DREAMPar['lik'] = 12  								                    # Gaussian likelihood
        sigma_I = 0.1
        sigma_t = 0.05  								                        # Specify measurement error
        p95_IT_SWIG = np.full((n_soil, 5, DREAMPar['d']), np.nan)               # Initialize return matrices
        ML_IT_SWIG = np.full((n_soil, DREAMPar['d']), np.nan)  			        # Initialize max likelihood values
        alg_IT_SWIG = P_IT_SWIG = FX95_IT_SWIG = None                           # Initialize results
        conv_IT_DREAM = MRstat_IT_DREAM = np.full((n_soil,1), np.nan)  	        # Initialize convergence diagnostic
        
        current_dir = os.getcwd()

        for st in range(n_soil):  							                    # Loop over each SWIG soil
            dat = data_SWIG[st]  							                    # Unpack data
            dat = dat[0]                                                        # Python added by JAV
            n_data = 2 * dat.shape[0]  						                    # Define plugin structure
            n_d = dat.shape[0]  							                    # Number of data points
            Meas_info = {'Sigma': np.concatenate([sigma_t * np.ones((n_d, 1)), sigma_I * np.ones((n_d, 1))], axis=0)}  # Measurement error std.
            plugin['t'] = dat[:, 0]  							                # Time (hours)
            plugin['I'] = dat[:, 1]  							                # Cumulative infiltration (cm)
            Meas_info['Y'] = np.concatenate([dat[:, 0], dat[:, 1]], axis=0)  	# Measurement data
            Meas_info['Sigma'] = Meas_info['Sigma'].squeeze()
            # Run DREAM-Suite package
            chain, output, FX, Z, logL = DREAM_Suite(method, Func_name, DREAMPar, Par_info, Meas_info, options, [], [], plugin)
            if alg_IT_SWIG is None:
                alg_IT_SWIG = np.full((output['MR_stat'].shape[0], 2, n_soil), np.nan)
            else:
                alg_IT_SWIG[:, :, st] = output['MR_stat']  
            parset = genparset(chain)  							                # Extract sample chains
            N = parset.shape[0]  							                    # Number of samples
            P_IT = parset[int(2 / 3 * N):, :DREAMPar['d'] + 2]  			    # Burn-in posterior samples
            M = P_IT.shape[0]  								                    # Number of posterior samples
            if P_IT_SWIG is None:
                P_IT_SWIG = np.full((M, DREAMPar['d'] + 2, n_soil), np.nan)
            else:
                P_IT_SWIG[:, :, st] = P_IT                        

            # Index ML parameters
            ii = np.argmax(np.sum(P_IT[:M, DREAMPar['d']:DREAMPar['d'] + 2], axis=1))
            ML_IT_SWIG[st, :] = P_IT[ii, :DREAMPar['d']]  		                # Store ML parameters
            FX = FX[int(2 / 3 * N):, :n_data]            				        # Get posterior simulations
            A = np.sort(FX, axis=0)  							                # Sort simulations
            n1 = int(0.025 * M)  							                    # 2.5% prediction limit
            n2 = int(0.5 * M)  								                    # Mean of prediction interval
            n3 = int(0.975 * M)  							                    # 97.5% prediction limit

            fx95_IT = np.zeros((n_data, 4))  	    				            # Initialize fx95 array
            fx95_IT[:, 0] = A[n1, :]  							                # Store 2.5%
            fx95_IT[:, 1] = A[n2, :]  							                # Store mean
            fx95_IT[:, 2] = FX[ii, :]  							                # Store ML simulation
            fx95_IT[:, 3] = A[n3, :]  							                # Store 97.5% prediction limit
            # if FX95_IT_SWIG is None:
            #     FX95_IT_SWIG = np.full((200, 4, n_soil), np.nan)              # Initialize results
            # else:
            #     FX95_IT_SWIG[:, :, st] = fx95_IT         				        # Store for soil

            for j in range(DREAMPar['d']):  						            # Do the same for parameters
                a = np.sort(P_IT[:, j])
                p95_IT_SWIG[st, :5, j] = [a[n1], np.mean(a), np.median(a), a[n3], np.std(a)]

            conv = 1 if output['MR_stat'][-1, 1] < 1.2 else 0  				    # Chains converged or not
            conv_IT_DREAM[st] = conv                                            # Store convergence
            MRstat_IT_DREAM[st] = output['MR_stat'][-1, 1]  				    # Store R_hat diagnostic

            # Save results
            mat_file = f"mat_files/{st}_{3}.mat"
            np.savez(mat_file, st=st, ML_IT_SWIG=ML_IT_SWIG, P_IT=P_IT, fx95_IT=fx95_IT, conv_IT_DREAM=conv_IT_DREAM, MRstat_IT_DREAM=MRstat_IT_DREAM)

        # Save optimal values
        np.savez("DREAM_IT_SWIG.mat", ML_IT_SWIG=ML_IT_SWIG, p95_IT_SWIG=p95_IT_SWIG, P_IT_SWIG=P_IT_SWIG, data_SWIG=data_SWIG, alg_IT_SWIG=alg_IT_SWIG, n_soil=n_soil, MRstat_IT_DREAM=MRstat_IT_DREAM, conv_IT_DREAM=conv_IT_DREAM)
