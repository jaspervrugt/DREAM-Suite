# ####################################################################### #
#                                                                         #
#             DDDDD     RRRRR     EEEEEEE     AAA     MM    MM            #
#             DDDDDD    RRRRRR    EEEEEEE    AAAAA    MM    MM            #
#             DD  DD    RR   RR   EE        AA   AA   MMM  MMM            #
#             DD   DD   RR  RR    EEEE      AA   AA   MMMMMMMM            #
#             DD   DD   RRRRR     EEEE      AAAAAAA   MMM  MMM            #
#             DD  DD    RR RR     EE        AAAAAAA   MM    MM            #
#             DDDDDD    RR  RR    EEEEEEE   AA   AA   MM    MM            #
#             DDDDD     RR   RR   EEEEEEE   AA   AA   MM    MM            #
#                                                                         #
#              SSSSSSSS  UU    UU   II   TTTTTTTTTT   EEEEEEE             #
#              SSSSSSS   UU    UU   II   TTTTTTTTTT   EEEEEEE             #
#              SS        UU    UU   II       TT       EE                  #
#              SSSS      UU    UU   II       TT       EEEE                #
#                 SSSS   UU    UU   II       TT       EEEE                #
#                   SS   UU    UU   II       TT       EE                  #
#               SSSSSS   UUUUUUUU   II       TT       EEEEEEE             #
#              SSSSSSS   UUUUUUUU   II       TT       EEEEEEE             #
#                                                                         #
# ####################################################################### #
#                                                                         #
# DREAM-Suite: DiffeRential Evolution Adaptive Metropolis algorithm       #
# with continuous and/or discrete variables and single- or multi-try      #
# sampling from an archive of current or past states using parallel       #
# direction, snooker and/or Kalman candidate points. This toolbox         #
# implements in one package the DREAM, DREAM_D, DREAM_ZS, DREAM_ZS,       #
# MTDREAM_ZS and DREAM_KZS algorithms.                                    #
#                                                                         #
# ####################################################################### #
#                                                                         #
# SYNOPSIS:                                                               #
#  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          #
#      DREAMPar,Par_info)                                                 #
#  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          #
#      DREAMPar,Par_info,Meas_info)                                       #
#  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          #
#      DREAMPar,Par_info,Meas_info,options)                               #
#  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          #
#      DREAMPar,Par_info,Meas_info,options,MAP_info)                      #
#  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          #
#      DREAMPar,Par_info,Meas_info,options,MAP_info,plugin)               #
# WHERE                                                                   #
#  method      [input] Name (string) of MCMC method                       #
#   = 'dream'                                                             #
#   = 'dream_zs'                                                          #
#   = 'dream_d'                                                           #
#   = 'dream_dzs'                                                         #
#  Func_name   [input] Function (string) returns (log)lik or sim values   #
#   → Func_name must return a likelihood if DREAMPar.lik = 1              #
#   → Func_name must return a log-likelihood if DREAMPar.lik = 2          #
#   → Func_name returns vector of n simulated values: DREAMPar.lik > 2    #
#  DREAMPar    [input] Dictionary with algorithmic variables              #
#   .d             Dimensionality (# variables) target distribution       #
#   .N             # of Markov chains                                     #
#   .T             # of generations (= # samples of each Markov chain)    #
#   .lik           Choice of likelihood function                          #
#     = 1          User computes likelihood in Func_name                  #
#     = 2          User computes log-likelihood in Func_name              #
#     = 11         Normal likelihood (= post. density) [Box and Tiao]     #
#     = 12         Normal likelihood & (non)const err Meas_info.Sigma     #
#     = 13         Normal likelihood & (non)const var. + AR(2) process    #
#     = 14         Genrlzd likelihood Schoups & Vrugt (2010) [OBSOLETE]   #
#     = 15         Whittle likelihood function (see Whittle, 1953)        #
#     = 16         Laplace likelihood & (non)const var. + AR(1) proc.     #
#     = 17         Student t likelihood & (non)const var. + AR(2) proc.   #
#     = 21         Approximate Bayesian Computation: normal likelihood    #
#     = 22         Approximate Bayesian Computation: boxcar function      #
#     = 23         Limits of Acceptability: log-liklhood is # pts LOA     #
#     = 31         GLUE: log-likelihood: Table 1a Beven & Freer 2001      #
#     = 32         GLUE: log-likelihood: Table 1b Beven & Freer 2001      #
#     = 33         GLUE: log-likelihood: Table 1c Beven & Freer 2001      #
#     = 34         GLUE: log-likelihood: Page 284 Beven & Binley 1992     #
#     = 44         Generalized likelihood ++: Vrugt et al. 2022           #
#     = 45         Universal likelihood: Vrugt et al. 2022                #
#     = 52         Log-likelihood: Generalized Least Squares (GLS) form   #
#     = 61         Laplace power likelihood: unit intgrl, lambda estmtd   #
#     = 62         Normal power likelihood: unit intgrl, lambda estmtd    #
#   .nCR           # crossover values                    DEF: 3           #
#   .delta         # chain pairs for proposal            DEF: 3           #
#   .lambda        Random error for ergodicity           DEF: 0.05        #
#   .zeta          Randomization                         DEF: 0.05        #
#   .p_unit_gamma  Probability unit jumprate (gamma)     DEF: 0.2         #
#   .adapt_pCR     Adapt crossover probabilities?        DEF: 'yes'       #
#   .thinning      Each thinning(th) chain sample stored DEF: 1           #
#   .GLUE          GLUE likelihood parameter             DEF: 10          #
#   .beta0         Scaling factor built-in jump rate     DEF: 1           #
#   .outlier       Outlier chain detection test          DEF: 'iqr'       #
#                   → DREAM and DREAM_D                                   #
#   .psnooker      Selection probability snooker jump    DEF: 0.1         #
#                   → DREAM_ZS, DREAM_DZS and MTDREAM_ZS                  #
#   .m0            Initial size of external archive, Z   DEF: 10*d        #
#                   → DREAM_ZS, DREAM_DZS and MTDREAM_ZS                  #
#   .k             Growth rate of external archive       DEF: 10          #
#                   → DREAM_ZS, DREAM_DZS and MTDREAM_ZS                  #
#   .mt            Number of multi-try proposals         DEF: 5           #
#                   → MTDREAM_ZS                                          #
#   .M             # samples archive Z for Kalman jump   DEF: 20          #
#                   → DREAM_KZS                                           #
#   .a_1           # Kalman jump begins at a_1 *.T gens  DEF: 0.1         #
#                   → DREAM_KZS                                           #
#   .a_2           # Kalman jump ends at a_2 *.T gens    DEF: 0.25        #
#                   → DREAM_KZS                                           #
#  Par_info    [input] Parameter structure: Ranges, initial/prior & bnd   #
#   .names         1xd-cell array with parameter names   DEF: []          #
#   .min           1xd-vector of min parameter values    DEF: -inf(1,d)   #
#   .max           1xd-vector of max parameter values    DEF: inf(1,d)    #
#   .norm          Work in normlzed parameter space      DEF: 0           #
#     = 0          Work in unnormlzed parameter space    DEFault          #
#     = 1          Work in normalized [0-1] parameter space               #
#   .boundhandling Treat the parameter bounds or not?                     #
#     = 'reflect'  Reflection method                                      #
#     = 'bound'    Set to bound                                           #
#     = 'fold'     Folding [Vrugt&Braak: doi:10.5194/hess-15-3701-2011]   #
#     = 'reject'   Reject out of bound proposals                          #
#     = 'none'     No boundary handling                  DEFault          #
#   .initial       Method to draw initial chain states                    #
#     = 'uniform'  Uniform: U(Par_info.min,Par_info.max)                  #
#     = 'latin'    Latin hypercube: LH(Par_info.min,Par_info.max)         #
#     = 'normal'   Normal:  N(Par_info.mu,Par_info.cov)                   #
#     = 'prior'    User specified prior distribution                      #
#     = 'user'     Initial chain states taken from Par_info.x0            #
#   .mu            1xd-mean vector: µ if .initial = 'normal'              #
#   .cov           dxd-covariance matrix: Σ if .initial = 'normal'        #
#   .x0            N x d matrix of initial states if .initial = 'user'    #
#   .prior         Prior distribution (manual) if .initial = 'prior'      #
#                  Ex 1: Par_info.prior = @(x,a,b) mvnpdf(x,a,b);         #
#                            Par_info.a = [-2 -2];                        #
#                            Par_info.b = eye(2);                         #
#                  Ex 2: Par_info.prior = {'normpdf(x,0,1)',...           #
#                                          'unifpdf(x,-2,2)'}             #
#                  Note: prior handle can return log(pdf): Code checks    #
#   .steps         d-vector with # intervals for each parameter           #
#                   → DREAM_D/DREAM_DZS                                   #
#  Meas_info   [input] Dictionary with measurement information (fitting)  #
#   .Y             nx1 vector against which model output is compared      #
#   .Sigma         Measurement error standard deviation of Y              #
#     = scalar     Constant measurement error                             #
#     = vector     Nonconstant measurement error (nx1 vector)             #
#   .sigma2        Treatment meas. error var. likelihood 13/16/17/44/45   #
#     = 'constant      Homoscedastic measurement error variance           #
#     = 'nonconstant'  Heteroscedastic measurement error variance         #
#   .S             Scalar/vector with summary metrics                     #
#   .C             dxd meas. err. cov. matrix GLS likelihood [= 52]       #
#   .R             Measurement error covariance matrix                    #
#                   → DREAM_KZS                                           #
#  options     [input] Structure with computational settings/options      #
#   .parallel      Multi-core computation chains?        DEF: 'yes'       #
#   .IO            If parallel, IO writing model?        DEF: 'no'        #
#   .rho           ABC dist. func. (anonymous handle)  DEF: @(X,Y)|X-Y|   #
#   .epsilon       ABC epsilon value (scalar/vector)     DEF: 0.025       #
#   .DB            Diagnostic Bayes?                     DEF: 'no'        #
#   .modout        Return model simulations?             DEF: 'no'        #
#   .save          Save DREAM output during the run?     DEF: 'no'        #
#   .restart       Restart run? (only with "save")       DEF: 'no'        #
#   .diagnostics   Compute within-chain diagnostics?     DEF: 'yes'       #
#   .print         Output writing screen (tables/figs)   DEF: 'yes'       #
#   .burnin        Burn-in # chain for conv. diagnstcs   DEF: 50          #
#  MAP_info    [input] Dictionary with information about MAP solution     #
#   .map           1xd vector with MAP solution                           #
#   .An            dxd sensitivity matrix MAP solution                    #
#   .Betan         dxd variability matrix MAP solution = Bn if unbiased   #
#  plugin      [input] 2nd input argument Func_name. Class set by user    #
#                                                                         #
#  chain       [outpt] T x (d+2) x N array N chains: pars+logpr+loglik    #
#                      If thinning is used then length ~ T/thinning       #
#  output      [outpt] Structure summarizes algorithmic performance       #
#   .R_stat        Univariate \hat{R} convergence diagnostic              #
#   .MR_stat       Multivariate \hat{R} convergence diagnostic            #
#   .AR            Acceptance rate (#)                                    #
#   .CR            Crossover selection probabilities                      #
#   .outlier       Information about outlier chains                       #
#   .RunTime       CPU time in seconds                                    #
#  FX          [outpt] T/thinning x n matrix model simuls chain samples   #
#  Z           [outpt] External archive used by DREAM_ZS and DREAM_DZS    #
#  loglik      [outpt] Log-likelihood sampled chains [diagnostics only]   #
#                                                                         #
# #####################################################################   #
#                                                                         #
# The different components/algorithms of DREAM-Suite are described in     #
#   Vrugt, J.A., R. de Punder, and P. Grünwald, A sandwich with water:    #
#       Bayesian/Frequentist uncertainty quantification under model       #
#       misspecification, Submitted to Water Resources Research,          #
#       May 2024, https://essopenarchive.org/users/597576/articles/...    #
#           937008-a-sandwich-with-water-bayesian-frequentist-...         #
#           uncertainty-quantification-under-model-misspecification       #
#   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and           #
#       Diagnostics: Elicitability, Propriety, and Scoring Rules for      #
#       Hydrograph Functionals, Water Resources Research, 60,             #
#       e2023WR036710, https://doi.org/10.1029/2023WR036710               #
#   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022),    #
#       On the use of distribution-adaptive likelihood functions:         #
#       Generalized and universal likelihood functions, scoring rules     #
#       and multi-criteria ranking, Journal of Hydrology, 615, Part B,    #
#       2022, https://doi.org/10.1016/j.jhydrol.2022.128542               #
#   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the     #
#       DREAM software package: Theory, concepts, and MATLAB              #
#       implementation, Environmental Modeling and Software, 75,          #
#       pp. 273-316, https://doi.org/10.1016/j.envsoft.2015.08.013        #
#   Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation   #
#       using Markov chain Monte Carlo simulation: DREAM_(ABC), Water     #
#       Resources Research, https://doi.org/10.1002/2014WR015386          #
#   Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model            #
#       calibration and evaluation: Approximate Bayesian computation,     #
#       Water Resources Research, 49, pp. 4335–4345,                      #
#           https://doi.org/10.1002/wrcr.20354                            #
#   Laloy, E., and J.A. Vrugt (2012), High-dimensional posterior          #
#       exploration of hydrologic models using multiple-try DREAM_(ZS)    #
#       and high-performance computing, Water Resources Research, 48,     #
#       W01526, https://doi.org/10.1029/2011WR010608                      #
#   Vrugt, J.A., and C.J.F. ter Braak (2011), DREAM_(D): An adaptive      #
#       Markov chain Monte Carlo simulation algorithm to solve            #
#       discrete, noncontinuous, and combinatorial posterior parameter    #
#       estimation problems, Hydrology and Earth System Sciences, 15,     #
#       pp. 3701-3713, https://doi.org/10.5194/hess-15-3701-2011          #
#   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and                        #
#       B.A. Robinson (2009), Equifinality of formal (DREAM) and          #
#       informal (GLUE) Bayesian approaches in                            #
#       hydrologic modeling? Stochastic Environmental Research and Risk   #
#       Assessment, 23(7), pp. 1011-1026,                                 #
#           https://doi.org/10.1007/s00477-008-0274-y                     #
#   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon,                #
#       B.A. Robinson, and J.M. Hyman (2009), Accelerating Markov chain   #
#       Monte Carlo simulation by differential evolution with             #
#       self-adaptive randomized subspace sampling, International         #
#       Journal of Nonlinear Sciences and Numerical Simulation, 10(3),    #
#       pp. 271-288                                                       #
#   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and            #
#       B.A. Robinson (2008), Treatment of input uncertainty in           #
#       hydrologic modeling: Doing hydrology backward with Markov chain   #
#       Monte Carlo simulation, Water Resources Research, 44, W00B09,     #
#       https://doi.org/10.1029/2007WR006720                              #
#   Ter Braak, C.J.F., and J.A. Vrugt (2008), Differential Evolution      #
#       Markov Chain with snooker updater and fewer chains, Statistics    #
#       and Computing, https://doi.org/10.1007/s11222-008-9104-9          #
#   Ter Braak, C.J.F. (2006), A Markov Chain Monte Carlo version of the   #
#       genetic algorithm differential evolution: easy Bayesian           #
#       computing for real parameter spaces, Statistics and Computing,    #
#       16, pp. 239-249, doi:10.1007/s11222-006-8769-1                    #
#                                                                         #
# ####################################################################### #
#                                                                         #
# DIFFERENT TEST EXAMPLES                                                 #
#  example 1: d-dimensional banana shaped Gaussian distribution           #
#  example 2: d-dimensional Gaussian distribution                         #
#  example 3: d-dimensional multimodal normal mixture distribution        #
#  example 4: real-world example rainfall-runoff (hymod in C++/MATLAB)    #
#  example 5: rainfall-runoff (hymod as external executable)              #
#  example 6: hmodel with distribution-adaptive likelihood functions      #
#  example 7: HYDRUS-1D soil hydraulic model: multiplicative prior        #
#  example 8: Approximate Bayesian Computation: Benchmark function        #
#  example 9: Spectral likelihood function in watershed modeling          #
#  example 10: Gaussian mixture distibution: multivariate prior           #
#  example 11: d-variate t-distribution: df ° freedom & corr. matrix R    #
#  example 12: pedometrics problem involving variogram fitting            #
#  example 13: Nash-Cascade hydrograph                                    #
#  example 14: Approx. Bayesian Comp. watershed signatures                #
#  example 15: Approx. Bayesian Comp. bivariate normal benchmark test     #
#  example 16: Hydrogeophysical inversion                                 #
#  example 17: Watershed model, normal, AR(1) and heteroscedastic lik.    #
#  example 18: Lotka-Volterra model: informal likelihood (GLUE)           #
#  example 19: Bayesian Model Averaging: I recommend MODELAVG toolbox!    #
#  example 20: Limits of acceptability: Soil temperature modeling         #
#  example 21: Limits of acceptability: Soil moisture model HYDRUS-1D     #
#  example 22: Limits of acceptability: Nash-Cascade hydrograph           #
#  example 23: Limits of acceptability: SAC-SMA (old C-code Euler int.)   #
#  example 24: Flow duration curve fitting                                #
#  example 25: Bedrock depth from high-res topo data & geomorph model     #
#  example 26: Data assimilation Lorenz model (SODA: Vrugt et al. 2005)   #
#  example 27: Data assimilation interception model (Vrugt et al. 2003)   #
#  example 28: Rainfall and hmodel parameter estimation from streamflow   #
#  example 29: Gaussian mixture distribution & multiplicative prior       #
#  example 30: Predator prey interactions                                 #
#  example 31: AR(2)-parameter estimation: Test of likelihood functions   #
#  example 32: Distribution-adaptive likelihood functions                 #
#  example 33: 2-dimensional rectangular target distribution              #
#  example 34: Haverkamp infiltration equation using HYDRUS-1D data       #
#  example 35: Haverkamp infiltration equation using SWIG database        #
#  example 99: Bayesian inference & M-estimation: Sandwich correction     #
#                                                                         #
# ####################################################################### #
#                                                                         #
# COPYRIGHT (c) 2024  the author                                          #
#                                                                         #
#   This program is free software: you can modify it under the terms      #
#   of the GNU General Public License as published by the Free Software   #
#   Foundation, either version 3 of the License, or (at your option)      #
#   any later version                                                     #
#                                                                         #
#   This program is distributed in the hope that it will be useful, but   #
#   WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      #
#   General Public License for more details                               #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  PYTHON CODE:                                                           #
#  © Written by Jasper A. Vrugt using GPT-4 OpenAI's language model       # 
#    University of California Irvine                                      #
#  Version 2.1    Dec 2024                                                #
#                                                                         #
# ####################################################################### #

# FURTHER CHECKING                                                        
#  Website:  http://faculty.sites.uci.edu/jasper                          
#  Papers: http://faculty.sites.uci.edu/jasper/publications/              
#  Google Scholar: https://scholar.google.com/citations?user=zkNXecUAAAAJ&hl=nl                                 


import numpy as np                                      
import os
import random
import multiprocess as mp                               
import shutil, array, re, math
import importlib
from scipy.stats import t, norm
from itertools import combinations                      
import scipy.linalg as la

## for miscellaneous files 
from scipy.signal import lfilter                        
from scipy.optimize import fsolve                       
from scipy.optimize import minimize
from scipy import optimize                                  
import scipy.special as sp
from scipy import integrate                             
from scipy.interpolate import interp1d                  
from scipy.stats import gaussian_kde                    
import pandas as pd                                     
from scipy import stats                                     

# for prior
from scipy.stats import multivariate_normal
from scipy.stats import gamma, uniform, genextreme, genpareto

## for postprocessing/plotting
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import matplotlib.patches as patches
from statsmodels.tsa.stattools import acf               
from numpy import percentile                                
from screeninfo import get_monitors


def DREAM_Suite_check(method, Func_name, DREAMPar, Par_info, Meas_info, options):
    ## ########################################################################
    ## This function verifies the input arguments of DREAM-Suite defined by  ##
    ## the user                                                              ##
    ##                                                                       ##
    ## SYNOPSIS: [DREAMPar,Par_info,options] = DREAM_Suite_check( ...        ##
    ##               method,Func_name,DREAMPar,Par_info,Meas_info,options)   ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                ##
    ## Los Alamos National Laboratory 			        	                 ##
    ##                                                                       ##
    ## ##################################################################### ##

    # Create an output file for warnings
    with open('warning_file.txt', 'w+') as fid:
        fid.write('-------------- DREAM-Suite WARNING file --------------\n')
        
        # Check input data structures
        for param, name in zip([DREAMPar, Par_info, Meas_info, options],['DREAMPar', 'Par_info', 'Meas_info', 'options']):
            if not isinstance(param, dict):
                raise ValueError(f"DREAM-Suite ERROR: input argument {name} should be a dictionary with fields")

        # Convert all strings to lowercase except parameter names [= names do later]
        for structure in ['DREAMPar', 'Par_info', 'options']:
            vrbl = eval(structure)
            for key, value in vrbl.items():
                if isinstance(value, str):
                    vrbl[key] = value.lower()        

        # <><><><><><><><><><><><><><><> Func_name <><><><><><><><><><><><><><><><>
        if not Func_name:
            raise ValueError("DREAM-Suite ERROR: The variable Func_name has to be defined as string (between quotes)")

        if isinstance(Func_name, (int, float)):
            raise ValueError(f"DREAM-Suite WARNING: The variable Func_name is defined as numerical value: This should be a string (between quotes) with name of MATLAB model script (.m file)")

        # <><><><><><><><><><><><><><><> DREAMPar <><><><><><><><><><><><><><><><><>
        if method in ['dream', 'dream_d']:
            if 'N' not in DREAMPar:
                raise ValueError("DREAM/DREAM_{D} ERROR: Field 'N' of structure DREAMPar undefined")
        if 'N' not in DREAMPar:
            raise ValueError("DREAM-Suite ERROR: Field 'N' of structure DREAMPar undefined")
        elif not isinstance(DREAMPar['N'], (int)):
            DREAMPar['N'] = int(DREAMPar['N'])       
        if 'd' not in DREAMPar:
            raise ValueError("DREAM-Suite ERROR: Field 'd' of structure DREAMPar undefined")
        elif not isinstance(DREAMPar['d'], (int)):
            DREAMPar['d'] = int(DREAMPar['d'])
        if 'T' not in DREAMPar:
            raise ValueError("DREAM-Suite ERROR: Field 'T' of structure DREAMPar undefined")
        elif not isinstance(DREAMPar['T'], (int)):
            DREAMPar['T'] = int(DREAMPar['T'])
        if 'lik' not in DREAMPar:
            raise ValueError("DREAM-Suite ERROR: Field 'lik' of structure DREAMPar undefined")
        elif not isinstance(DREAMPar['lik'], (int)):
            DREAMPar['lik'] = int(DREAMPar['lik'])
        if DREAMPar['d'] <= 0:
            raise ValueError("DREAM-Suite ERROR: Number of parameters should be integer and larger than zero -> Set DREAMPar.d >= 1")

        if 'delta' in DREAMPar:
            if not isinstance(DREAMPar['delta'], (int)):
                DREAMPar['delta'] = int(DREAMPar['delta'])
            if DREAMPar['delta'] <= 0:
                raise ValueError("DREAM-Suite ERROR: Number of chains pairs used for offspring should be integer and larger than zero -> Use at least DREAMPar.delta ∈ [1,5] (default: 1)")
            elif DREAMPar['delta'] > 10:
                evalstr = f"DREAM-Suite WARNING: Field 'delta' of structure DREAMPar set rather large -> recommend to use values of delta in [1,10]\n"
                print(evalstr)
                fid.write(evalstr)

        if method in ['dream', 'dream_d']:
            delta = DREAMPar.get('delta', 3)
            N_min = max([8, 2 * delta + 1, (DREAMPar['d'] + 1) // 2])
            if DREAMPar['N'] < N_min:
                raise ValueError(f"DREAM/DREAM_(D) ERROR: Insufficient number of chains with delta = {delta} -> Use at least DREAMPar.N = {N_min} chains")

            N_max = max([15, 2 * delta + 1, DREAMPar['d']])
            if DREAMPar['N'] > N_max:
                evalstr = f"DREAM/DREAM_(D) WARNING: Field 'N' of structure DREAMPar set rather large --> DREAMPar.N = {N_max} more than sufficient\n"
                print(evalstr)
                fid.write(evalstr)

        # Other methods (dream_zs, dream_dzs, dream_kzs)
        if method in ['dream_zs', 'dream_dzs', 'dream_kzs']:
            delta = DREAMPar.get('delta', 1)
            if 'm0' in DREAMPar:
                m0_min = max(10 * DREAMPar['d'], max(3 * DREAMPar['N'], 2 * delta * DREAMPar['N']))
                if DREAMPar['m0'] < m0_min:
                    raise ValueError(f"DREAM_(ZS)/DREAM_(DZS)/DREAM_(KZS) ERROR: Field 'm0' of structure DREAMPar should be integer and set to at least {m0_min}")
                elif DREAMPar['m0'] > 50 * DREAMPar['d']:
                    evalstr = f"DREAM_(ZS)/DREAM_(DZS)/DREAM_(KZS) WARNING: Field 'm0' of structure DREAMPar set rather large --> recommend to use DREAMPar.m0 = 10 * DREAMPar.d = {10 * DREAMPar['d']}\n"
                    print(evalstr)
                    fid.write(evalstr)

            if 'k' in DREAMPar:
                if DREAMPar['k'] < 1:
                    raise ValueError("DREAM_(ZS)/DREAM_(DZS)/DREAM_(KZS) ERROR: Field 'k' of structure DREAMPar should be integer and set to at least 1 -> DREAMPar.k in [1,20] (default: DREAMPar.k = 10)")
                elif DREAMPar['k'] > 25:
                    evalstr = f"DREAM_(ZS)/DREAM_(DZS)/DREAM_(KZS) WARNING: Field 'k' of structure DREAMPar set rather large --> recommend to use DREAMPar.k = 10 (default)\n"
                    print(evalstr)
                    fid.write(evalstr)

        # Additional validation checks for other fields in DREAMPar
        if DREAMPar['T'] < 2:
            raise ValueError("DREAM-Suite ERROR: Number of generations smaller than one -> Set at least DREAMPar.T = 2")
        elif DREAMPar['T'] > 1e6:
            evalstr = "DREAM-Suite WARNING: Field 'T' of structure DREAMPar set rather large"
            print(evalstr)
            fid.write(evalstr)

        if 'nCR' in DREAMPar:
            # Python: Make sure nCR is an integer in case the user specified its value
            DREAMPar['nCR'] = int(DREAMPar['nCR'])
            if DREAMPar['nCR'] <= 0:
                raise ValueError("DREAM-Suite ERROR: Number of crossover values used for offspring should be integer and larger than zero -> Use at least DREAMPar.nCR = 1 (default: 3)")
            elif DREAMPar['nCR'] < max(DREAMPar['d'] // 5, 1):
                evalstr = f"DREAM_Suite SUGGESTION: Number of crossover values of {DREAMPar['nCR']} defined in field 'nCR' of structure DREAMpar is rather small for DREAMPar.d = {DREAMPar['d']}\n"
                print(evalstr)
                fid.write(evalstr)
                evalstr = f"DREAM_Suite SUGGESTION: Try another trial later using DREAMPar.nCR = max(floor(DREAMPar.d/5),1) = {max(DREAMPar['d'] // 5, 1)} (rule of thumb) as this should enhance the sampling and convergence speed of DREAM\n"
                print(evalstr)
                fid.write(evalstr)
            elif DREAMPar['nCR'] >= 1 / 4 * DREAMPar['d']:
                evalstr = f"DREAM-Suite WARNING: Number of crossover values of {DREAMPar['nCR']} defined in field 'nCR' of structure DREAMpar is rather large for DREAMPar.d = {DREAMPar['d']}\n"
                print(evalstr)
                fid.write(evalstr)
                evalstr = f"DREAM_Suite SUGGESTION: Start another trial using DREAMPar.nCR = max(floor(DREAMPar.d/5),1) = {max(DREAMPar['d'] // 5, 1)} (rule of thumb) as this should enhance the sampling and convergence speed of DREAM\n"
                print(evalstr)
                fid.write(evalstr)

        if 'thinning' in DREAMPar:
            if not isinstance(DREAMPar['thinning'], (int)):
                DREAMPar['thinning'] = int(DREAMPar['thinning'])    
            if DREAMPar['thinning'] < 1:
                raise ValueError("DREAM-Suite ERROR: Thinning parameter should be integer and larger than zero -> Set DREAMPar.thinning >= 1 (default: 1)")
            elif DREAMPar['thinning'] > 20:
                evalstr = f"DREAM-Suite WARNING: Field 'thinning' of structure DREAMPar set rather large -> recommend to use DREAMPar.thinning in [1,20]\n"
                print(evalstr)
                fid.write(evalstr)

        if 'beta0' in DREAMPar:
            if DREAMPar['beta0'] <= 0:
                raise ValueError("DREAM-Suite ERROR: Multiplier of jump rate should be larger than zero -> Set DREAMPar.beta0 > 0 (default: 1)")
            elif DREAMPar['beta0'] > 3:
                evalstr = f"DREAM-Suite WARNING: Field 'beta0' of structure DREAMPar set rather large -> recommend to use values of beta0 in [0.1,3]\n"
                print(evalstr)
                fid.write(evalstr)

        if 'GLUE' in DREAMPar:
            if DREAMPar['GLUE'] <= 0:
                # ERROR -- GLUE likelihood variable larger than zero
                raise ValueError("DREAM-Suite ERROR: Likelihood variable of 'GLUE' should be larger than zero -> Set DREAMPar['GLUE'] > 0")
            elif DREAMPar['GLUE'] > 1000:
                evalstr = f"DREAM-Suite WARNING: Field 'GLUE' of structure DREAMPar set rather large\n"
                # Print warning to screen and to file
                print(evalstr)
                fid.write(evalstr)

        if 'lambda' in DREAMPar:
            if DREAMPar['lambda'] <= 0:
                # ERROR -- lambda should be positive
                raise ValueError("DREAM-Suite ERROR: Value of 'lambda' should be larger than zero -> Set DREAMPar['lambda'] > 0 (default: 0.05)")
            elif DREAMPar['lambda'] > 1:
                evalstr = "DREAM-Suite WARNING: Field 'lambda' of structure DREAMPar set rather large -> recommend to use DREAMPar['lambda'] in [0.01, 0.25]\n"
                # Print warning to screen and to file
                print(evalstr)
                fid.write(evalstr)

        if 'zeta' in DREAMPar:
            if DREAMPar['zeta'] <= 0:
                # ERROR -- zeta should be positive
                raise ValueError("DREAM-Suite ERROR: Value of 'zeta' should be larger than zero -> Set DREAMPar['zeta'] > 0 (default: 1e-12)")
            elif DREAMPar['zeta'] > 0.1:
                evalstr = "DREAM-Suite WARNING: Field 'zeta' of structure DREAMPar set rather large -> recommend to use DREAMPar['zeta'] in [1e-20, 1e-3]\n"
                # Print warning to screen and to file
                print(evalstr)
                fid.write(evalstr)

        if DREAMPar['lik'] not in [1, 2, 11, 12, 13, 14, 15, 16, 17, 21, 22, 23, 31, 32, 33, 34, 44, 45]:
            # Unknown built-in likelihood function
            raise ValueError("DREAM-Suite ERROR: Unknown choice of likelihood function -> Select DREAMPar['lik'] = {1,2,11,12,13,14,15,16,17,21,22,23,31,32,33,34,44,45} or 99 (own likelihood)")

        if 'p_unit_gamma' in DREAMPar:
            if DREAMPar['p_unit_gamma'] < 0 or DREAMPar['p_unit_gamma'] > 1:
                # ERROR -- unit jump rate probability between 0 and 1
                raise ValueError("DREAM-Suite ERROR: Probability of unit jump rate should be between 0 and 1 -> Set DREAMPar['p_unit_gamma'] in [0,1] (default: 0.2)")

        if method in ['dream_zs', 'dream_dzs']:
            if 'psnooker' in DREAMPar:
                if DREAMPar['psnooker'] < 0 or DREAMPar['psnooker'] > 1:
                    # ERROR -- snooker jump probability must be between 0 and 1
                    raise ValueError("DREAM_(ZS)/DREAM_(DZS)/DREAM_(KZS) ERROR: Probability of snooker jump should be between 0 and 1 -> Set DREAMPar['psnooker'] in [0,1] (default: 0.1)")

        if method == 'dream_kzs':
            if DREAMPar['lik'] not in [11, 12, 13, 14, 15, 16, 17, 31, 32, 33, 34, 44, 45]:
                # ERROR -- DREAM_KZS requires simulated model output
                raise ValueError("DREAM_(KZS) ERROR: Kalman jump cannot be used: model output is a (log)-density and not a simulation")
            if 'R' not in Meas_info and 'Sigma' not in Meas_info:
                if DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:
                    # This implementation has not been done in MATLAB code yet [distribution-adaptive --> std_e --> Kalman jump]
                    warning_str = f"DREAM_(KZS) WARNING: Field 'R' and field 'Sigma' of structure Meas_info are not specified -> Data measurement errors in Kalman jump set equal to std_e from likelihood function {DREAMPar['lik']}"
                    print(warning_str)
                    fid.write(warning_str)
                else: # WARNING -- DREAM_KZS will not consider the measurement error in Kalman jump
                    warning_str = "DREAM_(KZS) WARNING: Field 'R' and field 'Sigma' of structure Meas_info are not specified -> Kalman jump will assume that data are observed without a measurement error"
                    print(warning_str)
                    fid.write(warning_str)

        if 'R' in Meas_info:
            if Meas_info['R'].shape[0] != len(Meas_info['Y']):
                # ERROR -- Measurement error covariance matrix has to have similar number of rows as length of data vector
                raise ValueError("DREAM_(KZS) ERROR: Measurement error covariance matrix, Meas_info['R'], has to have similar number of rows as length of data vector")
            if Meas_info['R'].shape[1] != len(Meas_info['Y']):
                # ERROR -- Measurement error covariance matrix has to have similar number of columns as length of data vector
                raise ValueError("DREAM_(KZS) ERROR: Measurement error covariance matrix, Meas_info['R'], has to have similar number of columns as length of data vector")

        if 'psnooker' in DREAMPar:
            if DREAMPar['psnooker'] < 0 or DREAMPar['psnooker'] > 0.4:
                # ERROR -- snooker jump probability must be between 0 and 0.4
                raise ValueError("DREAM_(KZS) ERROR: Probability of snooker jump should be between 0 and 0.4 -> Set DREAMPar['psnooker'] in [0,0.4] (default: 0.1)")

        if 'pkalman' in DREAMPar:
            if DREAMPar['pkalman'] < 0 or DREAMPar['pkalman'] > 0.5:
                # ERROR -- kalman jump probability must be between 0 and 0.5
                raise ValueError("DREAM_(KZS) ERROR: Probability of kalman jump should be between 0 and 0.5 -> Set DREAMPar['pkalman'] in [0,1] (default: 0.4)")

        if 'psnooker' in DREAMPar and 'pkalman' in DREAMPar:
            if DREAMPar['psnooker'] + DREAMPar['pkalman'] >= 1:
                # ERROR -- parallel direction jump probability will be zero or negative
                raise ValueError("DREAM_(KZS) ERROR: Sum of selection probabilities of parallel direction and snooker jump cannot be larger than unity")

        # Note: Kalman jump activated between a_1*T and a_2*T, where a_1 default = 0.1 and a_2 default = 0.25
        # Kalman jump is not reversible? Please check paper by Zhang and Vrugt, WRR, 2020
        if 'a_1' in DREAMPar:
            if DREAMPar['a_1'] < 0:
                # ERROR -- a_1 cannot be smaller than zero
                raise ValueError("DREAM_(KZS) ERROR: Value of a_1 cannot be smaller than zero")
            if DREAMPar['a_1'] >= 1:
                # ERROR -- a_1 cannot be larger than one
                raise ValueError("DREAM_(KZS) ERROR: Value of a_1 cannot exceed value of 0.5")

        if 'a_2' in DREAMPar:
            if DREAMPar['a_2'] < 0:
                # ERROR -- a_2 cannot be smaller than zero
                raise ValueError("DREAM_(KZS) ERROR: Value of a_2 cannot be smaller than zero")
            if DREAMPar['a_2'] >= 1:
                # ERROR -- a_2 cannot be larger than one
                raise ValueError("DREAM_(KZS) ERROR: Value of a_2 cannot exceed value of 0.5")

        if 'a_1' in DREAMPar and 'a_2' in DREAMPar:
            if DREAMPar['a_2'] <= DREAMPar['a_1']:
                # ERROR -- a_1 cannot be larger than a_2
                raise ValueError("DREAM_(KZS) ERROR: Value of a_1 cannot exceed or be equal to value of a_2")

        if DREAMPar['lik'] == 23:
            # Warning that multiple epsilon values can be used
            if 'epsilon' in options:
                if np.isscalar(options['epsilon']):
                    options['epsilon'] = np.array([options['epsilon']]) 
                elif isinstance(options['epsilon'], array.array) or not isinstance(options['epsilon'], (int, float)):
                    options['epsilon'] = np.array(options['epsilon']).reshape(1, -1)
                if isinstance(Meas_info['S'], array.array) or not isinstance(Meas_info['S'], (int, float)):
                    Meas_info['S'] = np.array(Meas_info['S'])
                    # Ensure Meas_info['S'] is a numpy array (if it's a list, convert it)
                if isinstance(Meas_info['S'], list):
                    Meas_info['S'] = np.array(Meas_info['S'])
                if len(Meas_info['S']) != len(options['epsilon']) and len(options['epsilon']) > 1:
                    # ERROR -- Meas_info.S incorrect length!!
                    raise ValueError("DREAM-Suite ERROR: Number of elements of 'epsilon' does not match that of Meas_info.S!!")
                if len(options['epsilon']) == 1:
                    evalstr = ("DREAM-Suite WARNING: Limits of Acceptability - All observations use same LOA - "
                            "You can define 'options.epsilon' as vector with a different value for each observation\n")
                    # Now print warning to screen and to file
                    print(evalstr)
                    with open('warning_file.txt', 'a') as fid:
                        fid.write(evalstr)
            else:
                evalstr = ("DREAM-Suite WARNING: Limits of Acceptability - Default value of 'epsilon' will be used (=0.025) "
                        "and applied to all observations stored in Meas_info.S\n")
                # Now print warning to screen and to file
                print(evalstr)
                with open('warning_file.txt', 'a') as fid:
                    fid.write(evalstr)

        # <><><><><><><><><><><><><><><><> Par_info <><><><><><><><><><><><><><><><>
        if 'initial' not in Par_info:
            raise ValueError("DREAM-Suite ERROR: Initial sampling distribution not defined -> Define Par_info.initial = 'latin' or 'uniform' or 'normal' or 'prior' or 'user' !!")

        if method in ['dream_d', 'dream_dzs']:
            if 'steps' not in Par_info:
                raise ValueError("DREAM_(D)/DREAM_(DZS) ERROR: The number of intervals of each parameter (integers) need to be provided -> Define Par_info.steps!!")
            else:
                if isinstance(Par_info['steps'], array.array) or not isinstance(Par_info['steps'], (int, float)):
                    Par_info['steps'] = np.array(Par_info['steps']).reshape(1,-1)
            if (Par_info['steps']).shape != (1, DREAMPar['d']):
                raise ValueError(f"DREAM_(D)/DREAM_(DZS) ERROR: The number of elements of Par_info.steps needs to be equal to {DREAMPar['d']}!!")
            if any(step < 3 for step in Par_info['steps'].reshape(-1,1)):
                raise ValueError("DREAM_(D)/DREAM_(DZS) ERROR: The values of Par_info.steps need to be integers and larger than 2!!")

        if Par_info['initial'] not in ['latin', 'uniform', 'normal', 'prior', 'user']:
            raise ValueError("DREAM-Suite ERROR: Initial sampling distribution unknown -> Set Par_info.initial = 'latin' or 'uniform' or 'normal' or 'prior' or 'user' !!")

        if Par_info['initial'] == 'latin':
            # ERROR -- if lhs is used -> requires explicit parameter ranges
            if 'min' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: Latin hypercube sampling selected but minimum parameter values not defined -> Set Par_info.min!!")
            if 'max' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: Latin hypercube sampling selected but maximum parameter values not defined -> Set Par_info.max!!")

        if Par_info['initial'] == 'uniform':
            # ERROR -- if uniform initial sampling is used -> requires explicit parameter ranges
            if 'min' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: Uniform initial sampling selected but minimum parameter values not defined -> Set Par_info.min!!")
            if 'max' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: Uniform initial sampling selected but maximum parameter values not defined -> Set Par_info.max!!")

        if Par_info['initial'] == 'normal':
            # ERROR -- if normal is used --> mean and covariance of this distribution need to be defined
            if 'mu' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: Normal distribution selected to sample from but unknown mean -> Define Par_info.mu!!")
            if 'cov' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: Normal distribution selected to sample from but unknown covariance -> Define Par_info.cov!!")

        if Par_info['initial'] == 'normal':
            # ERROR -- if normal is used --> mean and covariance of this distribution need to be defined
            if isinstance(Par_info['mu'], array.array) or not isinstance(Par_info['mu'], (int, float)):
                Par_info['mu'] = np.array(Par_info['mu']).reshape(1, -1)
            if Par_info['mu'].shape != (1, DREAMPar['d']):
                raise ValueError(f"DREAM-Suite ERROR: Mean of normal distribution (Par_info.mu) should be a row vector with DREAMPar.d = {DREAMPar['d']} values")
            if Par_info['cov'].shape != (DREAMPar['d'], DREAMPar['d']):
                raise ValueError(f"DREAM-Suite ERROR: Covariance of normal distribution ('Par_info.cov') should be a square matrix of size DREAMPar.d x DREAMPar.d = {DREAMPar['d']} x {DREAMPar['d']} values")

        if Par_info['initial'] == 'prior':
            # ERROR -- if explicit prior is used --> marginals need to be defined
            if 'prior' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: Prior distribution selected but unknown field 'prior' of structure Par_info -> Define Par_info.prior (see manual)!!")

        if Par_info['initial'] == 'user':
            # ERROR -- if user initial sampling is used -> requires explicit starting points
            if 'x0' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: User initial sampling selected but starting points chain not defined -> Set Par_info.x0!!")
            else:
                if Par_info['x0'].shape[0] != DREAMPar['N']:
                    raise ValueError("DREAM-Suite ERROR: Number of rows of matrix Par_info.x0 does not equal DREAMPar.N")
                if Par_info['x0'].shape[1] != DREAMPar['d']:
                    raise ValueError("DREAM-Suite ERROR: Number of columns of matrix Par_info.x0 does not equal DREAMPar.d")

        # If field boundhandling not defined --> define no use
        if 'boundhandling' not in Par_info:
            Par_info['boundhandling'] = 'none'

        # Check: do we sample in normalized [0-1] space, or not?
        if 'norm' in Par_info:
            Par_info['norm'] = int(Par_info['norm'])
            if not isinstance(Par_info['norm'], (int, float)):
                raise ValueError("DREAM-Suite ERROR: Par_info.norm should be a scalar -> Define Par_info.norm = 0 or 1")
            if Par_info['norm'] not in [0, 1]:
                raise ValueError("DREAM-Suite ERROR: Par_info.norm should be zero or one -> Define Par_info.norm = 0 or 1")
            if Par_info['norm'] == 1:
                if 'min' not in Par_info:
                    raise ValueError("DREAM-Suite ERROR: Parameter normalization is used but minimum parameter values not defined -> Set Par_info.min!!")
                if 'max' not in Par_info:
                    raise ValueError("DREAM-Suite ERROR: Parameter normalization is used but maximum parameter values not defined -> Set Par_info.max!!")
        else:
            Par_info['norm'] = (0)

        if 'min' in Par_info:
            # Ensure that 'min' and 'max' are consistent with AMALGAMPar['d']
            if isinstance(Par_info['min'], array.array) or not isinstance(Par_info['min'], (int, float)):
                Par_info['min'] = np.array(Par_info['min']).reshape(1,-1)
            if Par_info['min'].shape[1] != DREAMPar['d']:
                raise ValueError(f"DREAM ERROR: Number of elements of field 'min' of structure Par_info should be equal to {DREAMPar['d']}!!")
        if 'max' in Par_info:
            if isinstance(Par_info['max'], array.array) or not isinstance(Par_info['max'], (int, float)):
                Par_info['max'] = np.array(Par_info['max']).reshape(1,-1)
            if Par_info['max'].shape[1] != DREAMPar['d']:
                raise ValueError(f"DREAM ERROR: Number of elements of field 'max' of structure Par_info should be equal to {DREAMPar['d']}!!")

        if 'boundhandling' in Par_info and Par_info['boundhandling'] != 'none':
            if 'min' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: Boundary handling is used but minimum parameter values not defined -> Set Par_info.min!!")
            elif Par_info['min'].ndim == 2: # Added to python
                Par_info['min'] = Par_info['min'].flatten()
            if 'max' not in Par_info:
                raise ValueError("DREAM-Suite ERROR: Boundary handling is used but maximum parameter values not defined -> Set Par_info.max!!")
            elif Par_info['max'].ndim == 2: # Added to python
                Par_info['max'] = Par_info['max'].flatten()

        # Handle boundary methods and parameter ranges
#        if any(x in Par_info['boundhandling'] for x in ['fold', 'bound', 'reflect', 'reject']):
#            if len(Par_info['min']) != DREAMPar['d']:
#                raise ValueError(f"DREAM-Suite ERROR: Number of elements of field 'min' of structure Par_info should be equal to {DREAMPar['d']}!!")
#            if len(Par_info['max']) != DREAMPar['d']:
#                raise ValueError(f"DREAM-Suite ERROR: Number of elements of field 'max' of structure Par_info should be equal to {DREAMPar['d']}!!")

        # Remove dummy variable to make things easier
        if Par_info['boundhandling'] == 'none':
            del Par_info['boundhandling']

        # Check parameter names
        if 'names' in Par_info:
            if len(Par_info['names']) != DREAMPar['d']:
                raise ValueError(f"DREAM-Suite ERROR: Number of elements of field 'names' of structure Par_info should be equal to {DREAMPar['d']}!!")

        # <><><><><><><><><><><><><><><><> Meas_info <><><><><><><><><><><><><><><><>
        if 'Sigma' in Meas_info:
            # Check content of Sigma
            if isinstance(Meas_info['Sigma'], (object)):
                Meas_info['Sigma'] = np.array([Meas_info['Sigma']])
            if isinstance(Meas_info['Sigma'], array.array):
                Meas_info['Sigma'] = np.array(Meas_info['Sigma'])
            if isinstance(Meas_info['Sigma'], (list, np.ndarray)):  # Handle it as a regular array
                if len(Meas_info['Sigma']) != len(Meas_info['Y']) and len(Meas_info['Sigma']) > 1:
                    # ERROR -- Meas_info.Sigma incorrect length!!
                    raise ValueError("DREAM-Suite ERROR: Length of Meas_info.Sigma is not equal to that of the observations stored in Meas_info.Y!!")
                # Now check  this as well
                if any(sigma < 0 for sigma in Meas_info['Sigma']):
                    print("DREAM-Suite WARNING: One or more entries of Meas_info.Sigma is negative - we use absolute values")
                    Meas_info['Sigma'] = np.abs(Meas_info['Sigma'])
                if any(sigma == 0 for sigma in Meas_info['Sigma']):
                    print("DREAM-Suite WARNING: One or more entries of Meas_info.Sigma are zero - we use smallest positive value of entered values")
                    id_zero = [i for i, sigma in enumerate(Meas_info['Sigma']) if sigma == 0]
                    nonzero_sigma = [sigma for sigma in Meas_info['Sigma'] if sigma > 0]
                    if nonzero_sigma:
                        min_value = min(nonzero_sigma)
                    else:
                        min_value = 1e-3
                    for id in id_zero:
                        Meas_info['Sigma'][id] = min_value
    
        if 'Y' in Meas_info:
            if isinstance(Meas_info['Y'], array.array):
                Meas_info['Y'] = np.array(Meas_info['Y'])

        if 'S' in Meas_info:
            if isinstance(Meas_info['S'], array.array):
                Meas_info['S'] = np.array(Meas_info['S'])

        # Check if 'print' field exists in options
        if 'print' in options:
            if not isinstance(options['print'], str):
                print("DREAM-Suite WARNING: Field 'print' of structure options should be a string (content equal to 'yes' or 'no')")
                # Set to 'yes' by default
                options['print'] = 'yes'

        # Check if 'epsilon' field exists in options
        if 'epsilon' in options:
            if np.isscalar(options['epsilon']):
                options['epsilon'] = np.array(options['epsilon'])
            if isinstance(options['epsilon'], array.array) or not isinstance(options['epsilon'], (int, float)):
                options['epsilon'] = np.array(options['epsilon']).reshape(1,-1)
            # Fix next sentence for epsilon is a vector
            if (options['epsilon'] <= 0).any():
                raise ValueError("DREAM-Suite ERROR: Value of 'epsilon' of structure options should be larger than zero -> Set options.epsilon > 0 (default 0.025)")

        # Check if 'rho' field exists in options
        if 'rho' in options:
            if isinstance(options['rho'], float) or isinstance(options['rho'], int):  # check if rho is real
                raise ValueError("DREAM-Suite ERROR: Field 'rho' of structure options should be defined as inline function -> Default: options.rho = inline('abs(X-Y)');")
            # Must evaluate the rho expression [different in Python]
            # options['rho'] = eval(options['rho'])

        # Check if 'burnin' field exists in options
        if 'burnin' in options:
            if options['burnin'] <= 0:
                raise ValueError("DREAM-Suite ERROR: Value of 'burnin' of structure options should be larger than zero -> Set options.burnin > 0 (default 50 [%])")
            if options['burnin'] >= 100:
                raise ValueError("DREAM-Suite ERROR: Value of 'burnin' of structure options should be smaller than hundred -> Set options.burnin < 100 (default 50 [%])")
            print("DREAM-Suite WARNING: Burn in percentage defined in field 'burnin' of structure options -> Stay between 50 - 80% - unless you are an expert")

        # Check content of each field in the structure options
        for key, value in options.items():
            # Check if content is 'yes' or 'no' for specific fields
            if key not in ['epsilon', 'rho', 'burnin']:
                if value not in ['yes', 'no']:
                    raise ValueError(f"DREAM-Suite ERROR: Field '{key}' of structure options should be set equal to 'yes' or 'no'")

    return DREAMPar, Par_info, options


def DREAM_Suite_setup(method, DREAMPar, Func_name, Par_info, Meas_info, options, LV):
    ## ##################################################################### ##
    ## This function initializes the main variables used in DREAM-Suite      ##
    ##                                                                       ##
    ## SYNOPSIS: [DREAMPar,Par_info,Meas_info,Lik_info,options] = ...        ##
    ##               DREAM_Suite_setup(method,Func_name,DREAMPar, ...        ##
    ##               Par_info,Meas_info,options)                             ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                ##
    ## Los Alamos National Laboratory 			                             ##
    ##                                                                       ##
    ## ##################################################################### ##

    # Seed random number generator
    np.random.seed(1 + round(100 * np.random.rand()))
    
    # Method specific configuration
    if method in ['dream', 'dream_d']:
        # Name variable
        name = ['nCR', 'delta', 'steps', 'lambda', 'zeta', 'p_unit_gamma', 
                'adapt_pCR', 'thinning', 'beta0', 'GLUE', 'outlier', 
                'pparallel', 'psnooker', 'pkalman', 'mt']
        
        # Default values algorithmic variables DREAM - if not specified
        value = ['3', '3', str(max(max(int(DREAMPar['T'] / 50), 1), 50)), 
                 '0.05', '1e-12', '0.2', "'yes'", '1', '1', '10', "'iqr'", 
                 '1', '0', '0', '1']

    elif method in ['dream_zs', 'dream_dzs']:
        # Name variable
        name = ['nCR', 'delta', 'steps', 'lambda', 'zeta', 'p_unit_gamma', 
                'adapt_pCR', 'thinning', 'beta0', 'GLUE', 'N', 'k', 'psnooker', 
                'm0', 'pparallel', 'psnooker', 'pkalman', 'mt']
        
        # Default values algorithmic variables DREAM - if not specified
        value = ['3', '3', str(max(max(int(DREAMPar['T'] / 50), 1), 50)),
                 '0.05', '1e-12', '0.2', "'yes'", '1', '1', '10', '3', '10', '0.1', 
                 f"max(10*DREAMPar['d'], max(3*DREAMPar['N'], 2*DREAMPar['delta']*DREAMPar['N']))", 
                 '0.9', '0.1', '0', '1']

    elif method == 'dream_kzs':
        # Name variable
        name = ['nCR', 'delta', 'steps', 'lambda', 'zeta', 'p_unit_gamma', 
                'adapt_pCR', 'thinning', 'beta0', 'GLUE', 'N', 'k', 'psnooker', 
                'm0', 'pparallel', 'psnooker', 'pkalman', 'M', 'a_1', 'a_2', 'mt']
        
        # Default values algorithmic variables DREAM - if not specified
        value = ['3', '3', str(max(max(int(DREAMPar['T'] / 50), 1), 50)), 
                 '0.05', '1e-12', '0.2', "'yes'", '1', '1', '10', '3', '10', '0.1', 
                 f"max(10*DREAMPar['d'], max(3*DREAMPar['N'], 2*DREAMPar['delta']*DREAMPar['N']))", 
                 '0.5', '0.1', '0.4', 
                 f"min(20, max(10*DREAMPar['d'], max(3*DREAMPar['N'], 2*DREAMPar['delta']*DREAMPar['N'])))", 
                 '0.1', '0.25', '1']

    elif method == 'mtdream_zs':
        # Check if snooker probability is defined by the user
        if 'psnooker' in DREAMPar:
            DREAMPar['pparallel'] = 1 - DREAMPar['psnooker']       
        # Otherwise, parallel direction probability
        elif 'pparallel' in DREAMPar:
            DREAMPar['psnooker'] = 1 - DREAMPar['pparallel']
        
        # Name variable
        name = ['N', 'delta', 'k', 'mt', 'psnooker', 'pparallel', 'pkalman', 'nCR', 
                'steps', 'lambda', 'zeta', 'p_unit_gamma', 'adapt_pCR', 
                'thinning', 'beta0', 'GLUE']
        
        # Default values algorithmic variables DREAM - if not specified
        value = ['3', '1', '10', '5', '0.1', '0.9', '0.0', '3', 
                 str(max(max(int(DREAMPar['T'] / 50), 1), 50)), '0.05', 
                 '1e-12', '0.2', "'yes'", '1', '1', '10']

    # Set default values
    for i, var_name in enumerate(name):
        if var_name not in DREAMPar:
            DREAMPar[var_name] = eval(value[i])

    if 'm0' in DREAMPar:
        DREAMPar['m0'] = int(DREAMPar['m0'])
        DREAMPar['m'] = DREAMPar['m0']
    elif method == 'mtdream_zs':
        DREAMPar['m0'] = max(10 * DREAMPar['d'], max(3 * DREAMPar['N'], 2 * DREAMPar['delta'] * DREAMPar['N']))
        DREAMPar['m'] = DREAMPar['m0']

    # Additional settings for options
    default_options = {'parallel': 'no', 'IO': 'no', 'modout': 'no', 'save': 'no', 'restart': 'no', 'DB': 'no',
        'epsilon': 0.025, 'diagnostics': 'yes', 'print': 'yes', 'burnin': 50}
    
    # Make sure steps is an integer in case it was submitted 
    DREAMPar['steps'] = int(DREAMPar['steps'])
    # Same for integer
    if 'k' in DREAMPar:
        DREAMPar['k'] = int(DREAMPar['k'])
    if 'mt' in DREAMPar:
        DREAMPar['mt'] = int(DREAMPar['mt'])

    for key, value in default_options.items():
        if key not in options:
            options[key] = value
    
    # Approximate Bayesian computation (ABC)
    if 20 < DREAMPar['lik'] < 24:
        options['ABC'] = 'yes'
        if 'rho' not in options:
            options['rho'] = 'lambda X, Y: abs(X - Y)'  # a = lambda X, Y: X - Y
    else:
        options['ABC'] = 'no'

    if method in ['dream', 'dream_d']:
        # Matrix DREAMPar.R: Store for each chain (as row) the index of all other chains available for DE
        DREAMPar['R'] = np.empty((DREAMPar['N'],DREAMPar['N']-1)).astype(int)
        for i in range(DREAMPar['N']):
            z = list(range(0,DREAMPar['N']))
            z.remove(i)
            DREAMPar['R'][i, 0:DREAMPar['N'] - 1] = z

        # Define psnooker and pkalman as zero
        DREAMPar['psnooker'] = DREAMPar['pkalman'] = 0
        # Calculate selection probability of parallel direction jump
        DREAMPar['pparallel'] = 1 - DREAMPar['psnooker'] - DREAMPar['pkalman']

    elif method in ['dream_zs', 'dream_dzs', 'dream_kzs']:
        # Determine sample indices for proposal generation
        if DREAMPar['delta'] == 1:
            # Define DREAMPar select
            DREAMPar['select'] = 3 * DREAMPar['N']
            # Fixed indices randomly chosen samples of Z
            DREAMPar['R'] = np.reshape(np.arange(0, DREAMPar['select']), (3, DREAMPar['N'])).T
        elif DREAMPar['delta'] > 1:
            # Define DREAMPar select
            DREAMPar['select'] = 2 * DREAMPar['delta'] * DREAMPar['N']
            # Fixed indices randomly chosen samples of Z
            DREAMPar['R'] = np.reshape(np.arange(0, DREAMPar['select']), (2 * DREAMPar['delta'], DREAMPar['N'])).T
        
        # If likelihood function larger than 21 (22 or 23) or DB then no snooker update allowed!
        if DREAMPar['lik'] in [22, 23] or options.get('DB') == 'yes':
            DREAMPar['psnooker'] = 0
        
        # Set Kalman jump probability to zero in current case group for DREAM_ZS or DREAM_DZS
        if method in ['dream_zs', 'dream_dzs']:
            DREAMPar['pkalman'] = 0
        
        # Calculate selection probability of parallel direction jump
        DREAMPar['pparallel'] = 1 - DREAMPar['psnooker'] - DREAMPar['pkalman']

    elif method == 'mtdream_zs':
        # Determine sample indices for proposal generation
        if DREAMPar['delta'] == 1:
            # Define DREAMPar select
            DREAMPar['select'] = 3 * DREAMPar['mt']
            # Fixed indices randomly chosen samples of Z
            DREAMPar['R'] = np.reshape(np.arange(0, DREAMPar['select']), (3, DREAMPar['mt']))
        elif DREAMPar['delta'] > 1:
            # Define DREAMPar select
            DREAMPar['select'] = 2 * DREAMPar['delta'] * DREAMPar['mt']
            # Fixed indices randomly chosen samples of Z
            DREAMPar['R'] = np.reshape(np.arange(0, DREAMPar['select']), (2 * DREAMPar['delta'], DREAMPar['mt']))
        # In DREAM_Suite we use DREAMPar.R = DREAMPar.R' !!!
    
    # Python - make entries an integer
    DREAMPar['R'] = DREAMPar['R'].astype(int)
        
    if method in ['dream_d', 'dream_dzs']:
        Par_info['step_size'] = (Par_info['max'] - Par_info['min']) / Par_info['steps']
    
    # Training data observations
    Meas_info['n'] = len(Meas_info.get('Y', []))
    
    # Handling summary metrics and epsilon
    if 'S' in Meas_info:
        Meas_info['n_S'] = len(Meas_info['S'])
        options['epsilon'] = np.ones(Meas_info['n_S']) * options['epsilon']
    else:
        Meas_info['n_S'] = int(0)

    # Handling Sigma: Either a vector of nx1 or otherwise None
    if 'Sigma' in Meas_info:
        if np.size(Meas_info['Sigma']) == 1:
            Meas_info['Sigma'] = np.ones((Meas_info['n'], 1)) * Meas_info['Sigma']
    else:
        # We define Sigma but it is None [this definition is key in Calc_proposal, Calc_likelihood !!]
        Meas_info['Sigma'] = None

    # If likelihood function is 52, compute determinant and inverse
    if DREAMPar['lik'] == 52:
        if 'C' in Meas_info:
            Meas_info['detC'] = np.linalg.det(Meas_info['C'])
            Meas_info['invC'] = np.linalg.inv(Meas_info['C'])
    
    # Define prior handle (random samples & evaluation of pdf)
    if 'prior' in Par_info:
        if len(Par_info['prior']) > 1 or DREAMPar['d'] == 1:    ## Univariate case
            Par_info['u'] = 'yes'
        else:                                                   ## Multivariate case
            Par_info['u'] = 'no'
 
        # Now determine whether user returned the pdf or log(pdf)
        Par_info['pr'] = check_prior(Par_info, DREAMPar, M = 100)  
        print(f"DREAM_Suite: Code analyzed that prior handle returns a {Par_info['pr']}")

    # Check if Kalman jump is used between a_1*T and a_2*T generations
    if method == 'dream_kzs':
        # Save original value of DREAMPar['pkalman']
        DREAMPar['oldpkalman'] = DREAMPar.get('pkalman', None)
        # Check the value of DREAMPar['a_1'] and modify 'pkalman'
        if DREAMPar.get('a_1', 0) > 0:
            DREAMPar['pkalman'] = 0

    # Calculate selection probability of parallel direction jump
    DREAMPar['pparallel'] = 1 - DREAMPar['psnooker'] - DREAMPar['pkalman']

    # Remove unnecessary fields based on method
    remove_fields = []
    if method in ['dream', 'dream_d']:
        remove_fields = ['k', 'm0', 'm']
    elif method in ['dream_zs', 'dream_dzs', 'dream_kzs']:
        remove_fields = ['outlier']
    elif method == 'mtdream_zs':
        remove_fields = []

    for field in remove_fields:
        if field in DREAMPar:
            del DREAMPar[field]  # Remove field from dictionary

    # Order and print to screen DREAMPar fields based on method
    if method in ['dream', 'dream_d']:
        field_order = ['d', 'N', 'T', 'lik', 'delta', 'nCR', 'lambda', 'zeta', 'p_unit_gamma', 'thinning', 'beta0', 
                    'GLUE', 'adapt_pCR', 'outlier', 'steps', 'R', 'pparallel', 'psnooker', 'pkalman', 'mt']
        not_print = 5
    elif method in ['dream_zs', 'dream_dzs']:
        field_order = ['d', 'N', 'T', 'lik', 'delta', 'k', 'pparallel', 'psnooker', 'm0', 'nCR', 'lambda', 'zeta', 
                    'p_unit_gamma', 'thinning', 'beta0', 'GLUE', 'adapt_pCR', 'm', 'steps', 'select', 'R', 'pkalman', 'mt']
        not_print = 5
    elif method == 'dream_kzs':
        field_order = ['d', 'N', 'T', 'lik', 'delta', 'k', 'pparallel', 'psnooker', 'pkalman', 'M', 'a_1', 'a_2', 
                    'm0', 'nCR', 'lambda', 'zeta', 'p_unit_gamma', 'thinning', 'beta0', 'GLUE', 'adapt_pCR', 'm', 
                    'steps', 'select', 'R', 'mt', 'oldpkalman']
        not_print = 3
    elif method == 'mtdream_zs':
        field_order = ['d', 'N', 'T', 'lik', 'delta', 'k', 'mt', 'pparallel', 'psnooker', 'm0', 'nCR', 'lambda', 'zeta', 
                    'p_unit_gamma', 'thinning', 'beta0', 'GLUE', 'adapt_pCR', 'm', 'steps', 'select', 'R']
        not_print = 5

    # Reorder fields in DREAMPar - next Line added by GPT
    DREAMPar = {key: DREAMPar[key] for key in field_order if key in DREAMPar}

    # Setup likelihood function
    Lik_info, DREAMPar, Par_info = setup_lik(Func_name, DREAMPar, Par_info, Meas_info, LV)

    # Print summary of the settings if flag is 0
    flag = 0
    if flag == 0:
        print("----------------------- Summary of DREAMPar ------------------------------")
        for key, value in DREAMPar.items():
            if callable(value):
                print(f"{key:<12}: {value}")        
            elif isinstance(value, np.ndarray):
                flattened_values = value.flatten()  
                # Check if the array is `DREAMPar['R']` and format it as integers
                if key == 'R':
                    formatted_values = ', '.join([f"{int(v)}" for v in flattened_values])   # Cast to int
                else:
                    formatted_values = ', '.join([f"{v:.4f}" for v in flattened_values])    # Default float formatting
                print(f"{key:<12}: {formatted_values}")
            elif isinstance(value, list):
                pr_str = '{' + ', '.join(map(str, value)) + '}'
                print(f"{key:<12}: {pr_str}")
            elif isinstance(value, int):
                print(f"{key:<12}: {value}")
            else:
                print(f"{key:<12}: {value:<6}")
        print("----------------------- Summary of Par_info ------------------------------")
        for key, value in Par_info.items():
            if callable(value):
                print(f"{key:<13}: {value}")
            elif isinstance(value, np.ndarray):         # Handle numpy arrays
                flattened_values = value.flatten()
                formatted_values = ', '.join([f"{v:.2f}" for v in flattened_values if isinstance(v, (int, float))])  # Only format numeric values
                print(f"{key:<13}: {formatted_values}")
            elif isinstance(value, list):               # Handle lists
                pr_str = '{' + ', '.join(map(str, value)) + '}'
                print(f"{key:<13}: {pr_str}")
            elif isinstance(value, (int, float)):       # Handle numeric types (int, float)
                if key == 'norm':
                    print(f"{key:<13}: {int(value)}")   # Format as integer
                else:
                    print(f"{key:<13}: {value:.2f}")    # Apply numeric formatting only for numbers
            else:
                print(f"{key:<13}: {value}")            # Just print the value for other types                
        print("----------------------- Summary of Lik_info ------------------------------")
        for key, value in Lik_info.items():
            if callable(value):
                print(f"{key:<13}: {value}")
            elif isinstance(value, np.ndarray):         # Handle numpy arrays
                flattened_values = value.flatten()
                formatted_values = ', '.join([f"{v:.2f}" for v in flattened_values if isinstance(v, (int, float))])  # Only format numeric values
                print(f"{key:<13}: {formatted_values}")
            elif isinstance(value, list):  # Handle lists
                pr_str = '{' + ', '.join(map(str, value)) + '}'
                print(f"{key:<13}: {pr_str}")
            elif isinstance(value, (int, float)):       # Handle numeric types (int, float)
                if key == 't_start' or key == 'method' or key == 'filename':
                    print(f"{key:<13}: {int(value)}")   # Format as integer
                else:
                    print(f"{key:<13}: {value:.2f}")    # Apply numeric formatting only for numbers
            else:
                print(f"{key:<13}: {value}")            # Just print the value for other types                
        print("----------------------- Summary of options -------------------------------")
        for key, value in options.items():
            if callable(value):
                print(f"{key:<11}: {value}")        
            elif isinstance(value, np.ndarray):  
                flattened_values = value.flatten() 
                formatted_values = ', '.join([f"{v:.4f}" for v in flattened_values])  
                print(f"{key:<11}: {formatted_values}")
            elif isinstance(value, list):
                pr_str = '{' + ', '.join(map(str, value)) + '}'
                print(f"{key:<11}: {pr_str}")
            else:
                print(f"{key:<11}: {value:<6}")
        print("--------------------------------------------------------------------------")

    return DREAMPar, Par_info, Meas_info, Lik_info, options


def DREAM_Suite_calc_setup(DREAMPar, Func_name, options):
    ## ##################################################################### ##
    ## This function initializes computational environment of DREAM-Suite    ##
    ##                                                                       ##
    ## SYNOPSIS: [DREAMPar,f_handle,T] = DREAM_Suite_calc_setup(...          ##
    ##               DREAMPar,fname,options,plugin)                          ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                ##
    ## Los Alamos National Laboratory                                        ##
    ##                                                                       ##
    ## ##################################################################### ##

    # Set initial generation
    T = 1
    base_dir = None

    # Set up parallel execution based on options
    if options['parallel'] == 'no':
        DREAMPar['CPU'] = 1  # Use 1 CPU (processor)
        
    elif options['parallel'] == 'yes':
        # How many available workers? 
        workers = mp.cpu_count() 
        if workers > DREAMPar['N']:
            workers = DREAMPar['N']
        
        DREAMPar['CPU'] = workers
        
        # Write parallelization status to file
        with open('warning_file.txt', 'a+') as fid:
            evalstr = f"DREAM_Suite PARALLEL: Python pool opened with {DREAMPar['CPU']} workers for {DREAMPar['N']} chains\n"
            print(evalstr)
            fid.write(evalstr)
        
        # Handle I/O directories for parallel execution
        if options['IO'] == 'yes':
            if os.name in ['nt', 'posix']:  # For Windows or Unix systems

                example_dir = os.getcwd()                                   # Current working directory
                base_dir = os.path.join(example_dir, "worker_dirs")         # Base directory for worker directories
                # copy all files to model files
                temp_dir = os.path.join(os.getcwd(), 'model_files')         # Copy files to temporary directory
                if not os.path.exists(temp_dir):
                    os.makedirs(temp_dir)

                # Copy files to temporary directory
                for filename in os.listdir(example_dir):
                    file_path = os.path.join(example_dir, filename)
                    if os.path.isfile(file_path):  # Ensure only files are copied
                        shutil.copy(file_path, temp_dir)

                model_files_dir = os.path.join(example_dir, "model_files")  # Directory with all model files

            # Ensure the model files directory exists
            if not os.path.exists(model_files_dir):
                print(f"Model files directory '{model_files_dir}' does not exist!")
                
                return

            # Step 1: Initial setup, copy model files into each worker directory (done only once)
            if not os.path.exists(base_dir):
                os.makedirs(base_dir)

            for task_id in range(DREAMPar['CPU']):
                worker_dir = os.path.join(base_dir, f'worker_{task_id}')
                # Only copy the model files to the worker directories if not already done
                if not os.path.exists(worker_dir):
                    copy_model_files(model_files_dir, worker_dir)

    ## Function handle = dynamic function [not pickable, but ok with multiprocessing]
    func_handle = get_function_handle(Func_name)
 
    return DREAMPar, T, func_handle, base_dir


def DREAM_Suite_initialize(method, DREAMPar, func_handle, Par_info, Meas_info, Lik_info, options, MAP_info, base_dir, plugin, printed_warnings):
    ## ##################################################################### ##
    ## Initializes main variables and chains of DREAM-Suite                  ##
    ##                                                                       ##
    ## SYNOPSIS: [chain,output,X,FX,S,Table_gamma,CR,pCR,tCR,TsdX_Xp, ...    ##
    ##               sdX_Xp,loglik,EX,std_mX,iloc,iter,gen,LLX, ...          ##
    ##               MAP_info] = DREAM_Suite_initialize(method, ...          ##
    ##               DREAMPar,Par_info,Meas_info,Lik_info,options,           ##
    ##               f_handle,MAP_info)                                      ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                ##
    ## Los Alamos National Laboratory                                        ##
    ##                                                                       ##
    ## ##################################################################### ##

    # Initialize variables based on method
    if method in ['dream_zs', 'dream_dzs', 'dream_kzs', 'mtdream_zs']:
        N = DREAMPar['m0'] + DREAMPar['N']
    elif method in ['dream', 'dream_d']:
        N = DREAMPar['N']

    # Initialize X (parameter space)
    X = np.full((N, DREAMPar['d']), np.nan)
    
    # Generate initial states of the chains based on Par_info.initial
    if Par_info['initial'] == 'uniform':    # Uniform sampling
        X = np.random.uniform(Par_info['min'], Par_info['max'], (N, DREAMPar['d']))
    elif Par_info['initial'] == 'latin':    # Latin hypercube sampling
        X = LH_sampling(Par_info['min'], Par_info['max'], N)
    elif Par_info['initial'] == 'normal':   # Multivariate normal
        X = np.tile(Par_info['mu'], (N, 1)) + np.random.randn(N, DREAMPar['d']) @ np.linalg.cholesky(Par_info['cov'])
    elif Par_info['initial'] == 'prior':
        if Par_info['u'] == 'yes':          # Univariate prior
            for qq in range(DREAMPar['d']):
                fnc = Par_info['prior'][qq]               
                for zz in range(N):
                    X[zz, qq] = fnc.rvs()
        else:                               # Multivariate prior
            fnc = Par_info['prior'][0]
            for zz in range(N):
                X[zz, :] = fnc.rvs()
    elif Par_info['initial'] == 'user':     # User specified initial chain states
        for zz in range(N):
            X[zz, :] = Par_info['x0'][zz, :]

    # Normalize X if needed
    if Par_info['norm'] == 1:
        # Forward transformation with original ranges
        if 'min' in Par_info:
            X = (X - Par_info['min']) / (Par_info['max'] - Par_info['min'])
        else:
            X = (X - np.min(X, axis = 0)) / (np.max(X, axis = 0) - np.min(X, axis = 0))
        # Now, we must store the original bounds for back-transformation
        Par_info['minun'], Par_info['maxun'] = Par_info['min'].copy(), Par_info['max'].copy()
        # New min and max are the unit interval
        Par_info['min'], Par_info['max'] = np.zeros(DREAMPar['d']), np.ones(DREAMPar['d'])
    
    # Boundary handling if specified
    if 'boundhandling' in Par_info:
        X, v = Boundary_handling(X, Par_info)
    else:
        v = np.ones(N) < 0

    if method in ['dream_d', 'dream_dzs']:
        X = Discrete_space(X, Par_info)

    # Handle DREAM methods
    if method in ['dream_zs', 'dream_dzs', 'dream_kzs', 'mtdream_zs']:
        Z = np.hstack([X[:DREAMPar['m0'], :], np.full((DREAMPar['m0'], 2), np.nan)])
        with open('Z.bin', 'wb') as fid_Z:
            fid_Z.write(Z.tobytes())

        X = X[DREAMPar['m0']:, :]
        v = v[DREAMPar['m0']:]

        if method == 'dream_kzs':
            FXZ = Evaluate_target(Z[:, :DREAMPar['d']], DREAMPar, func_handle, Meas_info, options, base_dir, plugin, printed_warnings)[0]
        else:
            FXZ = np.full((DREAMPar['m0'], Meas_info['n'] + Meas_info['n_S']), np.nan)
        with open('FXZ.bin', 'wb') as fid_FXZ:
            fid_FXZ.write(FXZ.T.tobytes())
            # Z.bin and FXZ.bin are now always synchronized in row number [not in MATLAB code?]

    elif method in ['dream', 'dream_d']:
        pass        # No past samples used

    # Calculate the log-likelihood and other statistics
    if 'map' in MAP_info:
        map_un = X_unnormalize(MAP_info['map'], Par_info)
        Fmap, _ = Evaluate_target(map_un, DREAMPar, func_handle, Meas_info, options, base_dir, plugin, printed_warnings)
        # I modified Calc_likelihood to include nargout as a new input variable
        _, _, _, _, _, LLmax = Calc_likelihood(map_un, Fmap, DREAMPar, Par_info, Meas_info, Lik_info, options, 6)
        MAP_info['loglik'] = np.sum(LLmax)
        MAP_info['map_un'] = map_un
        MAP_info['Fisher'] = MAP_info['An']
        MAP_info['Godambe'] = (MAP_info['An'] @ np.linalg.inv(MAP_info['Betan'])) @ MAP_info['An']
        
    # Normalize parameter values if necessary
    X_un = X_unnormalize(X, Par_info)
    FX, SX = Evaluate_target(X_un, DREAMPar, func_handle, Meas_info, options, base_dir, plugin, printed_warnings)
    logPR_X = Eval_prior(X_un, SX, Par_info, Meas_info, options)
    # I modified Calc_likelihood to include nargout as a new input variable
    logL_X, EX, std_mX, _, _, LLX = Calc_likelihood(X_un, FX, DREAMPar, Par_info, Meas_info, Lik_info, options, 6, MAP_info)

    X = np.hstack([X, np.zeros((X.shape[0], 2))])
    X[:, DREAMPar['d']:DREAMPar['d'] + 2] = np.column_stack([logPR_X.reshape(-1,1), logL_X.reshape(-1,1)])
    X[np.array(v), DREAMPar['d'] + 1] = -np.inf

    # Store results
    if SX is None:
        prt_result = FX
    else:
        prt_result = np.vstack([FX, SX])

    DREAM_Suite_store_results(options, prt_result.T, Meas_info, 'wb')

    # Initialize log-likelihood matrix
    loglik = np.full((DREAMPar['T'] + 1, DREAMPar['N'] + 1), np.nan)

    loglik[0, :DREAMPar['N'] + 1] = np.hstack([DREAMPar['N'], X[:DREAMPar['N'], DREAMPar['d'] + 1]])

    # Prepare output structure
    nr = int(DREAMPar['T'] // DREAMPar['steps'] + 1)
    output = {'AR': np.full((nr, 2), np.nan), 'MR_stat': np.full((nr, 2), np.nan), 
              'R_stat': np.full((nr, DREAMPar['d'] + 1), np.nan), 'CR': np.full((nr, DREAMPar['nCR'] + 1), np.nan), 
              'outlier': None}
    output['AR'][0, 0] = DREAMPar['N']
    output['R_stat'][0, 0] = DREAMPar['N']
    output['MR_stat'][0, 0] = DREAMPar['N']

    pCR = (1 / DREAMPar['nCR']) * np.ones(DREAMPar['nCR'])
    CR = np.reshape(np.random.choice(np.arange(1, DREAMPar['nCR'] + 1), size=DREAMPar['N'] * DREAMPar['steps'], p=pCR), 
                    (DREAMPar['N'], DREAMPar['steps']))

    tCR = np.ones(DREAMPar['nCR'])
    TsdX_Xp = np.ones(DREAMPar['nCR'])
    sdX_Xp = np.zeros((DREAMPar['N'], DREAMPar['steps']))

    output['CR'][0, :DREAMPar['nCR'] + 1] = np.hstack([DREAMPar['N'], pCR])

    # Initialize chain trajectory array
    chain = np.full((int(DREAMPar['T'] // DREAMPar['thinning'] + 1), DREAMPar['d'] + 2, DREAMPar['N']), np.nan)
    chain[0, :DREAMPar['d'] + 2, :] = np.reshape(X.T, (1, DREAMPar['d'] + 2, DREAMPar['N']))

    # Set jump rates for proposal distributions
    Table_gamma = np.full((DREAMPar['d'], DREAMPar['delta']), np.nan)
    for zz in range(1,DREAMPar['delta']+1):
        Table_gamma[:, zz - 1] = 2.38 / np.sqrt(2 * zz * (np.arange(1, DREAMPar['d'] + 1)))
    Table_gamma *= DREAMPar['beta0']
    # Initialize counters
    #iloc, iter, gen = 1, 2, 2   # MATLAB code
    iloc, iter, gen = 0, 1, 0

    if method == 'mtdream_zs':
        Nmt = DREAMPar['mt'] * DREAMPar['N']
        id_r = np.arange(0, Nmt)
        id_c = np.arange(DREAMPar['mt'] - 1, Nmt, step = DREAMPar['mt'])
        id_r = np.setdiff1d(id_r, id_c)
    else:
        id_r = id_c = None

    if 'sigma2' in Meas_info or Meas_info['Sigma'] != None:
        method_mapping = {0: "s0/s1 not considered", 1: "CONSTANT measurement error: s0 fixed at default value", 
                          2: "CONSTANT measurement error: s0 treated as phantom variable", 3: "CONSTANT measurement error: s0 estimated",
                          4: "CONSTANT measurement error: s0 estimated", 5: "NONCONSTANT measurement error: s0 fixed at default value and s1 treated as phantom variable", 
                          6: "NONCONSTANT measurement error: s0 fixed at default value and s1 estimated", 
                          7: "NONCONSTANT measurement error: s0 estimated and s1 treated as phantom variable", 
                          8: "NONCONSTANT measurement error: s0 estimated and s1 estimated"}
        print(f"DREAM_Suite: {method_mapping.get(Lik_info['method'], 'Unknown method')}")
        
    print('\n')
    print('  -----------------------------------------------------------------------------------------------------------------------                ')
    print('  DDDDDDDD   RRRRRR     EEEEEEEEEE     AAA     MMM        MMM      SSSSSSSSSS UUU     UUU   III   TTTTTTTTTTT  EEEEEEEEEE                ')
    print('  DDDDDDDDD  RRREERRR   EEEEEEEEE     AAAAA    MMMM      MMMM      SSSSSSSSS  UUU     UUU   III   TTTTTTTTTTT  EEEEEEEEE                 ')
    print('  DDD    DDD RRR   RRR  EEE          AAA AAA   MMMMM    MMMMM      SSS        UUU     UUU   III       TTT      EEE                       ')
    print('  DDD    DDD RRR   RRR  EEE         AAA   AAA  MMMMMM  MMMMMM      SSS        UUU     UUU   III       TTT      EEE                       ')
    print('  DDD    DDD RRR  RRR   EEEEEE     AAA     AAA MMM MMMMMM MMM ---- SSSSSSSSSS UUU     UUU   III       TTT      EEEEEE         /^ ^\      ')
    print('  DDD    DDD RRR RRR    EEEEEE     AAAAAAAAAAA MMM  MMMM  MMM ---- SSSSSSSSSS UUU     UUU   III       TTT      EEEEEE        / 0 0 \     ')
    print('  DDD    DDD RRRRRR     EEE        AAA     AAA MMM   MM   MMM             SSS UUU     UUU   III       TTT      EEE           V\ Y /V     ')
    print('  DDD    DDD RRR  RRR   EEE        AAA     AAA MMM        MMM             SSS UUU     UUU   III       TTT      EEE            / - \      ')
    print('  DDDDDDDDD  RRR   RRR  EEEEEEEEE  AAA     AAA MMM        MMM       SSSSSSSSS UUUUUUUUUUU   III       TTT      EEEEEEEEE     /     |     ')
    print('  DDDDDDDD   RRR    RRR EEEEEEEEEE AAA     AAA MMM        MMM      SSSSSSSSSS UUUUUUUUUUU   III       TTT      EEEEEEEEEE    V__) ||     ')
    print('  -----------------------------------------------------------------------------------------------------------------------                ')
    print('  © Jasper A. Vrugt, University of California Irvine & GPT-4 OpenAI''s language model')
    print('    ________________________________________________________________________')
    print('    Version 2.1, Feb. 2025, Beta-release: MATLAB implementation is benchmark')
    print('\n')

    # print('\n')
    # print('  -------------------------------------------------------------------------------------------------------------------------------------------------------                ')
    # print('          DDDDDDDD   RRRRRR     EEEEEEEEEE     AAA     MMM        MMM      PPPPPPPPPP      AAA     CCCCCCCCC KKK    KKK      AAA     GGGGGGGGG EEEEEEEEEE                ')
    # print('          DDDDDDDDD  RRREERRR   EEEEEEEEE     AAAAA    MMMM      MMMM      PPPPPPPPPPP    AAAAA    CCCCCCCCC KKK   KKK      AAAAA    GGGGGGGGG EEEEEEEEE                 ')
    # print('          DDD    DDD RRR   RRR  EEE          AAA AAA   MMMMM    MMMMM      PPP     PPP   AAA AAA   CCC       KKK  KKK      AAA AAA   GGG   GGG EEE                       ')
    # print('  EEEEEEE DDD    DDD RRR   RRR  EEE         AAA   AAA  MMMMMM  MMMMMM      PPP     PPP  AAA   AAA  CCC       KKK KKK      AAA   AAA  GGG   GGG EEE                       ')
    # print('  EEEEEEE DDD    DDD RRR  RRR   EEEEEE     AAA     AAA MMM MMMMMM MMM ---- PPP     PPP AAA     AAA CCC       KKKKK       AAA     AAA GGG   GGG EEEEEE         /^ ^\      ')
    # print('  EEE     DDD    DDD RRR RRR    EEEEEE     AAAAAAAAAAA MMM  MMMM  MMM ---- PPPPPPPPPP  AAAAAAAAAAA CCC       KKK KKK     AAAAAAAAAAA  GGGGGGGG EEEEEE        / 0 0 \     ')
    # print('  EEEEE   DDD    DDD RRRRRR     EEE        AAA     AAA MMM   MM   MMM      PPPPPPP     AAA     AAA CCC       KKK  KKK    AAA     AAA       GGG EEE           V\ Y /V     ')
    # print('  EEE     DDD    DDD RRR  RRR   EEE        AAA     AAA MMM        MMM      PPP         AAA     AAA CCC       KKK   KKK   AAA     AAA       GGG EEE            / - \      ')
    # print('  EEEEEEE DDDDDDDDD  RRR   RRR  EEEEEEEEE  AAA     AAA MMM        MMM      PPP         AAA     AAA CCCCCCCCC KKK    KKK  AAA     AAA  GGGGGGGG EEEEEEEEE     /     |     ')
    # print('  EEEEEEE DDDDDDDD   RRR    RRR EEEEEEEEEE AAA     AAA MMM        MMM      PPP         AAA     AAA CCCCCCCCC KKK     KKK AAA     AAA GGGGGGGGG EEEEEEEEEE    V__) ||     ')
    # print('  -------------------------------------------------------------------------------------------------------------------------------------------------------                ')
    # print('  © Jasper A. Vrugt, University of California Irvine & GPT-4 OpenAI''s language model')
    # print('    ________________________________________________________________________')
    # print('    Version 2.1, Feb. 2025, Beta-release: MATLAB implementation is benchmark')
    # print('\n')

    if method == 'dream':
        print('%4s DDDDDD   RRRRRR   EEEEEEE   AAA   MM      MM         ' % '')
        print('%4s DD   DD  RR   RR  EE       AA AA  MM      MM         ' % '')
        print('%4s DD    DD RR   RR  EE      AA   AA MMM    MMM         ' % '')
        print('%4s DD    DD RR   RR  EEEEEE  AA   AA MMMM  MMMM         ' % '')
        print('%4s DD    DD RRRRRR   EE      AAAAAAA MM  MM  MM         ' % '')
        print('%4s DD   DD  RR   RR  EE      AA   AA MM      MM         ' % '')
        print('%4s DDDDDD   RR    RR EEEEEEE AA   AA MM      MM         ' % '')

    elif method == 'dream_d':
        print('%4s DDDDDD   RRRRRR   EEEEEEE   AAA   MM      MM         ' % '')
        print('%4s DD   DD  RR   RR  EE       AA AA  MM      MM         ' % '')
        print('%4s DD    DD RR   RR  EE      AA   AA MMM    MMM         ' % '')
        print('%4s DD    DD RR   RR  EEEEEE  AA   AA MMMM  MMMM         ' % '')
        print('%4s DD    DD RRRRRR   EE      AAAAAAA MM  MM  MM DDDDDD  ' % '')
        print('%4s DD   DD  RR   RR  EE      AA   AA MM      MM DD   DD ' % '')
        print('%4s DDDDDD   RR    RR EEEEEEE AA   AA MM      MM DD   DD ' % '')
        print('%4s                                              DD   DD ' % '')
        print('%4s                                              DDDDDD	 ' % '')

    elif method == 'dream_zs':
        print('%4s DDDDDD   RRRRRR   EEEEEEE   AAA   MM      MM                ' % '')
        print('%4s DD   DD  RR   RR  EE       AA AA  MM      MM                ' % '')
        print('%4s DD    DD RR   RR  EE      AA   AA MMM    MMM                ' % '')
        print('%4s DD    DD RR   RR  EEEEEE  AA   AA MMMM  MMMM                ' % '')
        print('%4s DD    DD RRRRRR   EE      AAAAAAA MM  MM  MM ZZZZZZZ SSSSSS ' % '')
        print('%4s DD   DD  RR   RR  EE      AA   AA MM      MM     ZZ  SS     ' % '')
        print('%4s DDDDDD   RR    RR EEEEEEE AA   AA MM      MM    ZZ   SSSSSS ' % '')
        print('%4s                                                ZZ        SS ' % '')
        print('%4s                                              ZZZZZZZ SSSSSS ' % '')
 
    elif method == 'dream_dzs':
        print('%4s DDDDDD   RRRRRR   EEEEEEE   AAA   MM      MM                        ' % '')
        print('%4s DD   DD  RR   RR  EE       AA AA  MM      MM                        ' % '')
        print('%4s DD    DD RR   RR  EE      AA   AA MMM    MMM                        ' % '')
        print('%4s DD    DD RR   RR  EEEEEE  AA   AA MMMM  MMMM                        ' % '')
        print('%4s DD    DD RRRRRR   EE      AAAAAAA MM  MM  MM DDDDDD  ZZZZZZZ SSSSSS ' % '')
        print('%4s DD   DD  RR   RR  EE      AA   AA MM      MM DD   DD     ZZ  SS     ' % '')
        print('%4s DDDDDD   RR    RR EEEEEEE AA   AA MM      MM DD   DD    ZZ   SSSSSS ' % '')
        print('%4s                                              DD   DD   ZZ        SS ' % '')
        print('%4s                                              DDDDDDD ZZZZZZZ SSSSSS ' % '')

    elif method == 'dream_kzs':
        print('%4s DDDDDD   RRRRRR   EEEEEEE   AAA   MM      MM                        ' % '')
        print('%4s DD   DD  RR   RR  EE       AA AA  MM      MM                        ' % '')
        print('%4s DD    DD RR   RR  EE      AA   AA MMM    MMM                        ' % '')
        print('%4s DD    DD RR   RR  EEEEEE  AA   AA MMMM  MMMM                        ' % '')
        print('%4s DD    DD RRRRRR   EE      AAAAAAA MM  MM  MM KK   KK ZZZZZZZ SSSSSS ' % '')
        print('%4s DD   DD  RR   RR  EE      AA   AA MM      MM KK  KK      ZZ  SS     ' % '')
        print('%4s DDDDDD   RR    RR EEEEEEE AA   AA MM      MM KKKK       ZZ   SSSSSS ' % '')
        print('%4s                                              KK  KK    ZZ        SS ' % '')
        print('%4s                                              KK   KK ZZZZZZZ SSSSSS ' % '')

    elif method == 'mtdream_zs':
        print('%4s MM      MM TTTTTTTT     DDDDDD   RRRRRR   EEEEEEE   AAA   MM      MM                ' % '')
        print('%4s MM      MM TTTTTTTT     DD   DD  RR   RR  EE       AA AA  MM      MM                ' % '')
        print('%4s MMM    MMM    TT        DD    DD RR   RR  EE      AA   AA MMM    MMM                ' % '')
        print('%4s MMMM  MMMM    TT   ---- DD    DD RR   RR  EEEEEE  AA   AA MMMM  MMMM                ' % '')
        print('%4s MM  MM  MM    TT   ---- DD    DD RRRRRR   EE      AAAAAAA MM  MM  MM ZZZZZZZ SSSSSS ' % '')
        print('%4s MM      MM    TT        DD   DD  RR   RR  EE      AA   AA MM      MM     ZZ  SS     ' % '')
        print('%4s MM      MM    TT        DDDDDD   RR    RR EEEEEEE AA   AA MM      MM    ZZ   SSSSSS ' % '')
        print('%4s                                                                        ZZ        SS ' % '')
        print('%4s                                                                      ZZZZZZZ SSSSSS ' % '')

    print('\n')
    # Final output
    return chain, output, X, FX, SX, Table_gamma, CR, pCR, tCR, TsdX_Xp, sdX_Xp, loglik, EX, std_mX, id_r, id_c, iloc, iter, gen, LLX, Par_info, MAP_info


def DREAM_Suite_end(DREAMPar, options, base_dir):
    # ####################################################################### #
    # This function closes pool of workers and removes IO files               #
    # ####################################################################### #

    # Close workers if there are multiple CPUs
    if DREAMPar['CPU'] > 1:
        # Close parallel pool: check in Python
        pass

        if options['IO'] == 'yes':  # If IO writing is enabled, remove directories
            # Step 3: Clean up (delete) the worker directories after all generations are done
            cleanup_worker_directories(base_dir, DREAMPar['CPU'])

    # Open the warning_file.txt file in append mode
    with open('warning_file.txt', 'a+') as fid:
        # Write final line of warning file
        fid.write('----------- End of DREAM-Suite WARNING file ----------\n')

    # Optionally open the warning file in the default editor for the platform
    # if platform.system() in ['Windows', 'Darwin']:  # Windows or macOS
    #     os.system('notepad warning_file.txt' if platform.system() == 'Windows' else 'open warning_file.txt')

    return


def DREAM_Suite_restart(Func_name, file_name):
    # ####################################################################### #
    # This function reloads variables for a restart run in Python adaptation  #
    # ####################################################################### #

    # Try-except block to handle loading failure
    try:
        # Load the restart data from the file `file_name`
        data = np.load(file_name + '.npy', allow_pickle=True)           # with shelve.open(file_name, 'r') as file:
        loaded_data = data.item()
        DREAMPar = loaded_data['DREAMPar']                              #   DREAMPar = file['DREAMPar']
        Par_info = loaded_data['Par_info']                              #   Par_info = file['Par_info']
        Meas_info = loaded_data['Meas_info']                            #   Meas_info = file['Meas_info'] 
        Lik_info = loaded_data['Lik_info']                              #   Lik_info = file['Lik_info']
        options = loaded_data['options']                                #   options = file['options']
        MAP_info = loaded_data['MAP_info']                              #   MAP_info = file['MAP_info']
        chain = loaded_data['chain']                                    #   chain = file['chain']
        output = loaded_data['output']                                  #   output = file['output']
        X = loaded_data['X']                                            #   X = file['X']
        FX = loaded_data['FX']                                          #   FX = file['FX']
        S = loaded_data['S']                                            #   S = file['S']
        Table_gamma = loaded_data['Table_gamma']                        #   Table_gamma = file['Table_gamma']
        CR = loaded_data['CR']                                          #   CR = file['CR']
        pCR = loaded_data['pCR']                                        #   pCR = file['pCR']
        cCR = loaded_data['cCR']                                        #   cCR = file['cCR']
        TdCR = loaded_data['TdCR']                                      #   TdCR = file['TdCR']
        sdX_Xp = loaded_data['sdX_Xp']                                  #   sdX_Xp = file['sdX_Xp']
        loglik = loaded_data['loglik']                                  #   loglik = file['loglik']
        EX = loaded_data['EX']                                          #   EX = file['EX']
        std_mX = loaded_data['std_mX']                                  #   std_mX = file['std_mX']
        id_r = loaded_data['id_r']                                      #   id_r = file['id_r']
        id_c = loaded_data['id_c']                                      #   id_c = file['id_c']
        iloc = loaded_data['iloc']                                      #   iloc = file['iloc']
        it = loaded_data['it']                                          #   it = file['it']
        g = loaded_data['g']                                            #   g = file['g']
        LLX = loaded_data['LLX']                                        #   LLX = file['LLX']
        t = loaded_data['t']                                            #   t = file['t']

        # Load the content from the pickle file
        #  with open("{}.pkl".format(file_name), 'rb') as file:
        #      plugin = pickle.load(file)

    except KeyError as e:
        # Handle missing key in the loaded data
        raise FileNotFoundError(f"DREAM-Suite ERROR: Cannot restart --> File {file_name} does not exist. "
                                 "Next run, to avoid this problem, set field 'save' of structure options to 'yes'")

    # Randomly sample the crossover values
    CR = np.random.choice(range(1, DREAMPar['nCR']+1), size=(DREAMPar['N'], DREAMPar['steps']), p=pCR)
    
    # Open a warning file
    with open('warning_file.txt', 'a+') as fid:
        # Check whether the previous run was aborted
        if t < DREAMPar['T']:
            warning_message = f"DREAM_Suite RESTART: Starting with t = {t}, but still using old budget of DREAMPar.T = {DREAMPar['T']}\n"
            print(warning_message)
            fid.write(warning_message)
        else:
            # Check for additional generations in T.txt
            fname_T = 'T.txt'
            if os.path.exists(fname_T):
                warning_message = f"DREAM_Suite RESTART: Located file 'T.txt' in the respective directory - now checking its content\n"
                print(warning_message), fid.write(warning_message)
                # Here, `checkfile_T` is assumed to be a function that reads 'T.txt' and returns additional generations
                T_new = checkfile_T(fname_T)
                # Assuming `T_new` is returned by checkfile_T
                warning_message = f"DREAM_Suite RESTART: User has requested additional generations in 'T.txt', completing {T_new} generations\n"
                print(warning_message), fid.write(warning_message)
                # Update total generations
                warning_message = f"DREAM_Suite RESTART: Initial t = {t} and completing {T_new} additional generations so that DREAMPar.T = {DREAMPar['T'] + T_new}\n"
                print(warning_message), fid.write(warning_message)
            else:
                # If T.txt doesn't exist, double the previous budget
                warning_message = "DREAM_Suite RESTART: Did not locate file 'T.txt'. Doubling the previous budget.\n"
                print(warning_message), fid.write(warning_message)
                T_new = DREAMPar['T']
                warning_message = f"DREAM_Suite RESTART: Initial t = {t} and completing {T_new} additional generations so that DREAMPar.T = {DREAMPar['T'] + T_new}\n"
                print(warning_message), fid.write(warning_message)
            
            # Append empty lines to various vectors, matrices and arrays
            loglik = np.pad(loglik, ((0, T_new), (0, 0)), mode='constant', constant_values=np.nan)
            nr = int(T_new // DREAMPar['steps'])    
            # Add nr lines to AR, pCR, R_stat and MR_stat
            output['AR'] = np.pad(output['AR'], ((0, nr), (0, 0)), mode='constant', constant_values=np.nan)
            output['CR'] = np.pad(output['CR'], ((0, nr), (0, 0)), mode='constant', constant_values=np.nan)
            output['R_stat'] = np.pad(output['R_stat'], ((0, nr), (0, 0)), mode='constant', constant_values=np.nan)
            output['MR_stat'] = np.pad(output['MR_stat'], ((0, nr), (0, 0)), mode='constant', constant_values=np.nan)
            nr = int(T_new // DREAMPar['thinning'])    
            # Add nr lines to chain
            chain = np.pad(chain, ((0, nr), (0, 0), (0, 0)), mode='constant', constant_values=np.nan)
            # Update DREAMPar.T with T_new
            DREAMPar['T'] += T_new

    # Define starting value of T
    T_start = t + 1

    # Initialize computational framework (assuming the function `DREAM_Suite_calc_setup` exists)
    DREAMPar, _, func_handle, base_dir = DREAM_Suite_calc_setup(DREAMPar, Func_name, options)

    # Return all the variables in a dictionary
    return DREAMPar, Par_info, Meas_info, Lik_info, func_handle, options, MAP_info, chain, output, X, FX, S, Table_gamma, CR, pCR, cCR, TdCR, sdX_Xp, loglik, EX, std_mX, id_r, id_c, iloc, it, g, LLX, base_dir, T_start


def DREAM_Suite_store_results(options, FX, Meas_info, id):
    # ####################################################################### #
    # Store the DREAM_Suite results to binary files                           #
    # ####################################################################### #
    """
    Parameters:
    - options: Dictionary containing options for the DREAM algorithm.
    - FX: Array of model simulations.
    - Meas_info: Dictionary containing measurement information, including n and n_S.
    - id: The identifier for the file to append to.

    Returns:
    - None. The function stores the results in a binary file called 'FX.bin'.
    """
    ### Consider writing this as a float 32
    # Check if we need to output the model simulations (modout == 'yes')
    if options.get('modout') == 'yes' and (Meas_info.get('n') > 0 or Meas_info.get('n_S') > 0):
        # Open the file 'FX.bin' in append mode
        if id == 'wb':      ## Open new file 
            with open('FX.bin', id) as fid_FX:  
                fid_FX.write(FX.tobytes())
        elif id == 'ab':    ## Append the new data
            with open('FX.bin', id) as fid_FX:  
                FX.tofile(fid_FX)            


def DREAM_Suite_postproc(method, DREAMPar, func_handle, Par_info, Meas_info, Lik_info, options, chain, output, iloc, MAP_info, base_dir, plugin, printed_warnings):
    ## ##################################################################### ##
    ## This function generates tables and figures of DREAM-Suite results     ##
    ##                                                                       ##
    ## SYNOPSIS: [chain,output,FX,Z] = DREAM_Suite_postproc(method,...       ##
    ##               DREAMPar,f_handle,Par_info,Meas_info,Lik_info,...       ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                ##
    ## Los Alamos National Laboratory 			        	                 ##
    ## University of California Irvine                                       ##
    ##                                                                       ##
    ## Release Version Jan. 2018     Publication of DREAM manual in EMS      ##
    ## Update          Nov. 2021     Changed output and DREAM_Package_end    ##
    ## Update          June 2023     Included new likelihood functions       ##
    ## Update          July 2023     Scoring rules and other additions       ##
    ## Version 2.0     July 2024     Major changes throughout code           ##
    ##                                                                       ##
    ## ##################################################################### ##

    # Print wait statement to the screen
    print('DREAM-Suite PREPARING OUTPUT: PLEASE WAIT ...')
    
    # Define name of program
    n_program = 'DREAM-Suite'
    
    # Define name of figures file
    file_name = f'{n_program}_figures.pdf'
    
    # Determine screen size (using matplotlib to get screen dimensions)
    monitor = get_monitors()[0]
    screen_width = monitor.width
    screen_height = monitor.height
    x_mult = screen_width / 1920
    y_mult = screen_height / 1080
    t_mult = min(x_mult, y_mult)

    # Define fontsize for figures
    fontsize_xylabel = 16 * t_mult
    fontsize_axis = 18 * t_mult
    fontsize_axis_numbers_small_graphs = 14 * t_mult
    fontsize_legend = 14 * t_mult
    fontsize_text = 14 * t_mult
    fontsize_title = 18 * t_mult
    fontsize_title_small_graphs = 15 * t_mult
    markersize_symbols = 5
    fontsize_titlepage = 20 * t_mult
    markersize_legend = 10 * t_mult
    fontsize_axis_numbers = 16 * t_mult
    length_minor_ticks = 2
    length_major_ticks = 4

    # Define color_order
    color_order = np.array([
        [0.00, 0.00, 1.00],
        [0.00, 0.50, 0.00],
        [1.00, 0.00, 0.00],
        [0.00, 0.75, 0.75],
        [0.75, 0.00, 0.75],
        [0.75, 0.75, 0.00],
        [0.25, 0.25, 0.25],
        [0.75, 0.25, 0.25],
        [0.95, 0.95, 0.00],
        [0.25, 0.25, 0.75],
        [0.75, 0.75, 0.75],
        [0.00, 1.00, 0.00],
        [0.76, 0.57, 0.17],
        [0.54, 0.63, 0.22],
        [0.34, 0.57, 0.92],
        [1.00, 0.10, 0.60],
        [0.88, 0.75, 0.73],
        [0.10, 0.49, 0.47],
        [0.66, 0.34, 0.65],
        [0.99, 0.41, 0.23]
    ])
    
    red_color = tuple(np.array([1.0, 0.0, 0.0]))
    black_color = tuple(np.array([0.0, 0.0, 0.0]))
    dark_gray = tuple(np.array([0.25, 0.25, 0.25]))
    medium_gray = tuple(np.array([0.5, 0.5, 0.5]))
    light_gray = tuple(np.array([0.75, 0.75, 0.75]))

    p_alfa = 0.05                               # significance level -> 95% confidence interval
    maxbins = 25;                               # Maximum number of bins for histograms?
    n_chains = min(DREAMPar['N'], 5)            # How many chains are plotted in relevant figures?

    # Extract min and max values from Par_info
    #Par_info['min'] = Par_info['min'][0, :]
    #Par_info['max'] = Par_info['max'][0, :]

    # Unpack the chains
    T, _, N = chain.shape
    # Create tables for posterior moments
    P = genparset(chain)                        # Assemble all chains in one matrix
    nP = P.shape[0]                             # Number of parameter samples
    P_post = P[int(3 / 4 * nP):nP, :DREAMPar['d'] + 2]
    if DREAMPar['lik'] != 23:                   # not limits of acceptability
        id_ML = np.argmax(P_post[:, DREAMPar['d'] + 1])
        ML = P_post[id_ML, :DREAMPar['d']]
        id_MAP = np.argmax(np.sum(P_post[:, DREAMPar['d']:DREAMPar['d'] + 2], axis = 1))
        MAP = P_post[id_MAP, :DREAMPar['d']]
        MAP_pr = P_post[id_MAP, DREAMPar['d'] ]
        MAP_lik = P_post[id_MAP, DREAMPar['d'] + 1]
    else:
        ML = np.nan * np.ones(DREAMPar['d'])
        MAP = np.nan * np.ones(DREAMPar['d'])
        MAP_pr = MAP_lik = np.nan

    MEAN = np.mean(P_post[:, :DREAMPar['d']], axis = 0)
    MED = np.median(P_post[:, :DREAMPar['d']], axis = 0)
    STD = np.std(P_post[:, :DREAMPar['d']], axis = 0)
    CORR = np.corrcoef(P_post[:, :DREAMPar['d']], rowvar=False)

    # Prepare variables for tables and figures
    str_par, str_table, str_CR, str_chain, str_S, method_table, method_fig, sndwch, sndwch_text = prepare_output(DREAMPar, Par_info, Meas_info, n_chains, MAP_info, method)

    # Print statement
    print("\n"), print('DREAM-Suite TABULATING RESULTS: PLEASE WAIT ...')

    # Write Tables with posterior statistics of parameters
    tabulate_output(method_table.upper(), DREAMPar, ML, MAP, MEAN, MED, STD, CORR, str_table, sndwch)
    
    # Write Tables with single chain convergence diagnostics
    tabulate_diagnostics(method_table.upper(), DREAMPar, options, chain, output, iloc, sndwch)
 
    # 2A. Read simulated output
    if options['modout'] == 'yes' and (Meas_info['n'] > 0 or Meas_info['n_S'] > 0):
        file_size = os.path.getsize('FX.bin')
        row_size = (Meas_info['n'] + Meas_info['n_S']) * np.dtype(np.float64).itemsize
        num_rows = file_size // row_size
        with open('FX.bin', 'rb') as fid_FX:
            raw_data = fid_FX.read()
            padding_size = (8 - len(raw_data) % 8) % 8
            padded_data = raw_data + b'\0' * padding_size
            FX = np.frombuffer(padded_data, dtype=np.float64)
        # FX = FX.reshape(num_rows,Meas_info['n'] + Meas_info['n_S'])
        FX = FX.reshape(nP,Meas_info['n'] + Meas_info['n_S'])
        # num_rows matches nP
    else:
        FX = None

    # 2B. Read past archive of samples
    if method in ['dream', 'dream_d']:
        Z = []
    elif method in ['dream_zs', 'dream_dzs', 'dream_kzs', 'mtdream_zs']:
        file_size = os.path.getsize('Z.bin')
        row_size = (DREAMPar['d'] + 2) * np.dtype(np.float64).itemsize
        num_rows = file_size // row_size
        with open('Z.bin', 'rb') as fid_Z:
            raw_data = fid_Z.read()
            padding_size = (8 - len(raw_data) % 8) % 8
            padded_data = raw_data + b'\0' * padding_size
            Z = np.frombuffer(padded_data, dtype=np.float64)
        # Z = Z.reshape(num_rows,DREAMPar['d']+2)
        Z = Z.reshape(DREAMPar['m'],DREAMPar['d']+2)
        # --> num_rows matches DREAMPar['m']

   # Now do we print figures or not?
    if options['print'] == 'no':
        return FX, Z
    else:
        pass  # Continue with figure output

    # Print wait statement to the screen
    print('DREAM-Suite PLOTTING RESULTS: PLEASE WAIT ...')

    # Check whether the model produces simulation or not
    sim_out = None  # Initialize simulation output
    S_post = None   # Initialize posterior samples
    FX_post = None  # Initialize posterior flux

    # Possibilities for model output:
    # 1. model output is likelihood, no Meas_info.Y, or Meas_info.S
    # 2. model output is log-likelihood, no Meas_info.Y or Meas_info.S
    # 3. model output are summary metrics, Meas_info.S (ABC)
    # 4. model output is simulation, Meas_info.Y, and/or Meas_info.S

    # First, create or unpack the model output simulations of the posterior samples
    if Meas_info['n'] > 0 or Meas_info['n_S'] > 0:
        if options['modout'] == 'yes':          # Check if the output was stored
            # If yes, take from FX (assuming FX is a numpy array)
            sim_out = FX[int(np.floor(3/4 * nP)):nP, :]
        else:
            # Generate posterior simulations
            if DREAMPar['d'] > 1:
                # Get unique posterior parameter vectorss
                UP_post, iiUP, idUP = np.unique(P_post[:, :DREAMPar['d']], axis = 0, return_index = True, return_inverse = True)
            else:
                UP_post, iiUP, idUP = np.unique(P_post[:, :DREAMPar['d']], return_index = True, return_inverse = True)
            # Now transform to unnormalized space, if necessary
            UP_postun = X_unnormalize(UP_post, Par_info)
            # Evaluate the model for posterior samples
            sim_out, S_out = Evaluate_target(UP_postun, DREAMPar, func_handle, Meas_info, options, base_dir, plugin, printed_warnings, verbose=1)
            # Duplicate simulations of unique parameter vectors
            if S_out is not None:
                sim_out = np.vstack([sim_out, S_out]).T
            else:
                sim_out = sim_out.T

            sim_out = sim_out[idUP, :]
            # Clear unused variable
            del UP_post  

    # Determine what sim_out actually is
    if options['DB'] == 'yes':  ## sim_out contains simulation and summary metrics
        FX_post = sim_out[:, :Meas_info.n]
        S_post = sim_out[:, Meas_info.n:Meas_info.n + Meas_info.n_S]
    else:                       ## Summary metrics or LOA
        if DREAMPar['lik'] in [21, 22, 23]:
            S_post = sim_out
        else:
            FX_post = sim_out

    sim_out = None 
    # Make sure that measured data have a single dimension only
    if 'n' in Meas_info:
        if Meas_info['n'] > 0:
            Y_meas = Meas_info['Y'].flatten()

    # Make sure that measured symmary metrics have a single dimension only
    if 'n_S' in Meas_info:
        if Meas_info['n_S'] > 0:
            if isinstance(Meas_info['S'], list):
                S_meas = np.array(Meas_info['S'])
            elif isinstance(Meas_info['S'], np.ndarray):
                if Meas_info['S'].ndim == 2:
                    S_meas = Meas_info['S'].flatten()
                else:
                    S_meas = Meas_info['S']

    # If FX_post is not empty, compute RMSE and PBIAS for MAP (Maximum A Posteriori) solution
    if FX_post is not None:
        FX_MAP = FX_post[id_MAP, :Meas_info['n']].T.copy()    # Extract MAP solution
        print(FX_MAP, FX_MAP.shape)
        RMSE_MAP = np.sqrt(np.sum((FX_MAP - Y_meas) ** 2) / Meas_info['n'])
        PBIAS_MAP = 100 * np.sum(FX_MAP - Y_meas) / np.sum(Y_meas)
        # Display RMSE and PBIAS
        print(RMSE_MAP, PBIAS_MAP)
    else:
        RMSE_MAP = PBIAS_MAP = FX_MAP = None

    # Extract t_start from Lik_info
    t_start = Lik_info['t_start']
    p_gamma = 1 - p_alfa

    # Now start plotting the figures
    with PdfPages(file_name) as pdf:

        ## --------------------------------------
        ## Plot Empty Figure for PDF
        ## --------------------------------------
        plt.figure(figsize=(12, 6))
        plt.plot([], [], 'ro')  # Empty plot
        plt.axis([0, 1, 0, 1])
        plt.gca().set_facecolor('w')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])        
        plt.gca().set_xticks([])
        plt.gca().set_yticks([])
        plt.text(0.25 * x_mult, 0.6 * y_mult, r'${\rm Visual \; results \; of \; DREAM-Suite \; toolbox}$', fontsize = fontsize_titlepage) #, ha = 'center', va = 'center')
        plt.text(0.27 * x_mult, 0.5 * y_mult, r'$\;\;\;\;\;{\rm Tables \; are \; not \; printed \; to \; PDF \; file}$', fontsize = fontsize_titlepage) #, ha = 'center', va = 'center') #, fontweight='bold')
        ax = plt.gca()  # Get current axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        pdf.savefig()    
        plt.show()

        ## --------------------------------------
        ## PLOT univariate scale reduction factor
        ## --------------------------------------        
        if 'R_stat' in output:
            fig, ax = plt.subplots(figsize=(15, 6))
            if np.any(np.sum(output['R_stat'][1:, 1:DREAMPar['d']] > 10, axis = 1) > 0):
                ha = ax.plot(output['R_stat'][:, 0]  / DREAMPar['N'], output['R_stat'][:, 1:], label = str_par, linewidth = 1.5)
            else:
                ha = ax.semilogy(output['R_stat'][:, 0]  / DREAMPar['N'], output['R_stat'][:, 1:], label = str_par, linewidth = 1.5)    
            # color_order = plt.cm.tab10.colors
            for i, line in enumerate(ha):
                line.set_color(color_order[i % len(color_order)])
            # Set x-axis to show integer ticks
            plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
            # Create custom legend handles with colored lines
            legend_handles = [Line2D([0], [0], color = color_order[i % len(color_order)], lw = 4, label = str_par[i]) for i in range(DREAMPar['d'])]
            legend = plt.legend(handles = legend_handles, fontsize = fontsize_legend, loc='best')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(color_order[z % len(color_order)])  
            ax.set_xlabel('Number of generations', fontsize = fontsize_xylabel, labelpad = 10)
            ax.set_ylabel(f'Scale reduction factor: $\\hat{{R}}$', fontsize = fontsize_xylabel, labelpad = 10)
            ax.set_title(f'{method_fig}: Univariate scale reduction diagnostic of GR (1992)', fontsize = fontsize_title)
            # ax.axhline(y=1.2, color = dark_gray, linestyle='--', linewidth = 1.5)
            x_max = output['R_stat'][-1, 0] / DREAMPar['N']
            y_max = np.ceil(np.max(np.minimum(output['R_stat'][1:, 1:],10)))
            ax.set_xlim([1, x_max])
            ax.set_ylim([0.95, y_max])
            ax.plot([1, x_max], [1.2, 1.2], '--', linewidth = 1.5, color = dark_gray)
            ax.set_ylim(bottom=0.95) 
            ax.minorticks_on()          
            # Check convergence
            check_conv = np.sum(output['R_stat'][1:, 1:DREAMPar['d'] + 1] > 1.2, axis = 1)
            if np.any(check_conv > 0):
                idx_conv = np.where(check_conv > 0)[0][-1] + 1
            else:
                idx_conv = 2
            if idx_conv:
                if idx_conv == (output['R_stat'].shape[0] - 1):
                    x_loc = 1.02 * output['R_stat'][-1, 0] / DREAMPar['N'] 
                    ax.text(x_loc, 1, 'not converged', fontsize = fontsize_legend, rotation=90, ha = 'center', color =  'red')
                else:
                    ax.axvline(output['R_stat'][idx_conv, 0] / DREAMPar['N'], color = 'gray', linestyle='--', linewidth = 1.5)
                    x_loc = ( output['R_stat'][idx_conv, 0] + 0.005 * output['R_stat'][-1, 0] ) / DREAMPar['N'] 
                    ax.text(x_loc, 1, 'converged', fontsize = fontsize_legend, rotation=90, ha = 'center', color =  'green')
            ax.tick_params(which='major', length = length_major_ticks, color = 'black')
            ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
            ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
            ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers)
            # Format the x and y axis with thousands separators
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            # ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            # ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
            pdf.savefig()
            plt.show() # (block=False)

        ## --------------------------------------
        ## PLOT multivariate scale reduction factor
        ## --------------------------------------        
        if 'MR_stat' in output:
            fig, ax = plt.subplots(figsize=(15, 6))
            if np.any(output['MR_stat'][1:, 1] > 10):
                ax.semilogy(output['MR_stat'][:, 0] / DREAMPar['N'], output['MR_stat'][:, 1], 'r', linewidth = 1.5)
            else:
                ax.plot(output['MR_stat'][:, 0] / DREAMPar['N'], output['MR_stat'][:, 1], 'r', linewidth = 1.5)
            ax.set_xlabel('Number of generations', fontsize = fontsize_xylabel, labelpad = 10)
            ax.set_ylabel(f"Scale reduction factor: $\\hat{{R}}^{{{{d}}}}$", fontsize = fontsize_xylabel, labelpad = 10)
            ax.set_title(f"{method_fig}: Multivariate scale reduction diagnostic of BG (1998)", fontsize = fontsize_title)            
            x_max = output['MR_stat'][-1, 0] / DREAMPar['N']
            y_max = np.ceil(np.max(np.minimum(output['MR_stat'][1:, 1],10)))
            ax.set_xlim([1, x_max])
            ax.set_ylim([0.95, y_max])
            ax.set_ylim(bottom=0.95)           
            ax.plot([1, x_max], [1.2, 1.2], '--', linewidth = 1.5, color = dark_gray)
            ax.minorticks_on()          
            ax.tick_params(which='major', length = length_major_ticks, color = 'black')
            ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
            ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
            ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers)
            # Format the x and y axis with thousands separators
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            # ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            # ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
            # save te figure
            pdf.savefig()
            plt.show() #(block=False)

        ## --------------------------------------
        ## PLOT Acceptance Rate 
        ## --------------------------------------        
        if 'AR' in output:
            fig, ax = plt.subplots(figsize=(15, 6))
            ax.plot(output['AR'][:, 0] / DREAMPar['N'], output['AR'][:, 1], 'r', linewidth = 1.5)
            ax.set_xlabel('Number of generations', fontsize = fontsize_xylabel, labelpad = 10)
            ax.set_ylabel('Acceptance rate', fontsize = fontsize_xylabel, labelpad = 10)
            ax.set_title(f"{method_fig}: Evolution of acceptance rate", fontsize = fontsize_title)
            x_max = output['AR'][-1, 0] / DREAMPar['N']
            ax.set_xlim([1, x_max])          
            ax.minorticks_on()          
            ax.tick_params(which='major', length = length_major_ticks, color = 'black')
            ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
            ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
            ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers)
            # Format the x and y axis with thousands separators
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            # ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            # ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
            # save te figure
            pdf.savefig()
            plt.show() # (block=False)

        ## --------------------------------------
        ## PLOT Selection probabilities
        ## --------------------------------------        
        if 'pCR' in output:
            fig, ax = plt.subplots(figsize=(15, 6))
            for i in range(DREAMPar['nCR']):
                ax.plot(output['pCR'][:, 0] / DREAMPar['N'], output['pCR'][:, i+1], label = str_CR[i], linewidth = 1.5)
            ax.set_xlabel('Number of generations', fontsize = fontsize_xylabel, labelpad = 10)
            ax.set_ylabel('Selection probability of crossover values', fontsize = fontsize_xylabel, labelpad = 10)
            ax.set_title(f"{method_fig}: Evolution of crossover selection probabilities", fontsize = fontsize_title)
            x_max = output['pCR'][-1, 0] / DREAMPar['N']
            ax.axis([1, x_max, 0, 1])
            ax.minorticks_on()          
            ax.tick_params(which='major', length = length_major_ticks, color = 'black')
            ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
            ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
            ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers)
            # Format the x and y axis with thousands separators
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            #ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            #ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
            # save te figure
            pdf.savefig()
            plt.show() # (block=False)

        ## --------------------------------------
        ## PLOT histograms of marginal densities of the parameters
        ## --------------------------------------
        row, col = 2, 4
        idx_y_label = np.arange(0, DREAMPar['d'] + 1 , col)
        Nbins = np.array([calcnbins(P_post[:, i]) for i in range(DREAMPar['d'])])
        nbins = min(np.min(Nbins), maxbins)
        nbins = max(5, round(nbins / 2))  # Adjust number of bins
        counter = 0
        while counter <= DREAMPar['d'] - 1:
            fig, axs = plt.subplots(row, col, figsize=(12, 8), zorder = 1)
            plt.subplots_adjust(wspace=0.4, hspace=0.4)  # Increase the space horizontally and vertically
            for ax in axs.flat:
                U, edges = np.histogram(P_post[:, counter], bins=nbins, density=True)
                X = 0.5 * (edges[:-1] + edges[1:])
                ax.bar(X, U / np.max(U), width=(edges[1] - edges[0]), color = 'gray', alpha=0.7,zorder=2)
                ax.set_xlabel(str_par[counter], fontsize = fontsize_xylabel, labelpad = 10)
                yticks = np.arange(0, 1.02, 0.2)
                ax.set_yticks(yticks)  # Set x-tick positions
                ax.set_yticklabels([str(round(tick,1)) for tick in yticks])
                ax.set_ylim(0, 1.02)                # Adjust y limits with padding
                if counter in idx_y_label:
                    ax.set_ylabel('Empirical density', fontsize = fontsize_xylabel, labelpad = 10)
                ax.plot(MAP[counter], 1.02, 'bx', markersize = 12, markeredgewidth=3, linewidth = 3, zorder=3, clip_on=False)
                label_plot = get_label(counter)  # This converts 0 -> 'A', 1 -> 'B', etc.
                ax.text(0.02,0.97, f'({label_plot})', transform = ax.transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left')
                ax.minorticks_on()          
                ax.tick_params(which='major', length = length_major_ticks, color = 'black')
                ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
                ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
                ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)
                # Format the x and y axis with thousands separators
                # ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
                # Optionally, you can set major ticks locator to improve the appearance
                # ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
                counter += 1
                if counter == DREAMPar['d']:
                    for ax in axs.flat:
                        if not ax.has_data():   # If the axis has no data, remove it
                            fig.delaxes(ax)     # delete empty subplots
                    break
            fig.suptitle(f"{method_fig}: Sampled marginal distribution of the parameters", fontsize = fontsize_title)
            pdf.savefig()
            plt.show() # (block=False)

        ## --------------------------------------
        ## PLOT the correlation matrix
        ## --------------------------------------
        if DREAMPar['d'] > 1:
            fig, ax = plt.subplots(figsize=(12, 8))
            cax = ax.imshow(CORR, cmap='viridis', aspect='auto')
            cbar = fig.colorbar(cax, ax=ax, orientation='vertical', pad = 0.01)
            cbar.set_label('r(x,y)', fontsize = fontsize_legend)
            ax.set_xticks(np.arange(DREAMPar['d']), str_par,fontsize = fontsize_xylabel)
            ax.set_yticks(np.arange(DREAMPar['d']), str_par,fontsize = fontsize_xylabel)
            ax.set_title(f"{method_fig}: Map of posterior correlation coefficients of parameters", fontsize = fontsize_title)
            pdf.savefig()
            plt.show() # (block=False)

        ## --------------------------------------       
        ## PLOT correlation plots of the posterior parameter samples
        ## --------------------------------------
        if DREAMPar['d'] in range(2, 15):
            fig, axs = plt.subplots(DREAMPar['d'], DREAMPar['d'], figsize=(15, 15))
            Nbins = np.array([calcnbins(P_post[:, i]) for i in range(DREAMPar['d'])])
            nbins = min(np.min(Nbins), maxbins)
            nbins = max(5, round(nbins / 2))  # Adjust number of bins
            for i in range(DREAMPar['d']):
                for j in range(DREAMPar['d']):
                    if i != j: # For scatter plots the unique posterior samples should do
                        axs[i, j].scatter(P_post[:, j], P_post[:, i], color = 'gray', s=12)
                        # Add a least-squares line for off-diagonal plots
                        fit = np.polyfit(P_post[:, j], P_post[:, i], 1)
                        axs[i, j].plot(P_post[:, j], np.polyval(fit, P_post[:, j]), 'b--', linewidth = 1)
                        axs[i, j].set_xlim([min(P_post[:, j]), max(P_post[:, j])])
                        axs[i, j].set_ylim([min(P_post[:, i]), max(P_post[:, i])])
                    if i == j: # For histograms we need all posterior samples
                        axs[i, j].hist(P_post[:, i], bins=nbins, density=True, alpha=0.6, color = 'gray', edgecolor = 'black')
                        axs[i, j].set_xlim([min(P_post[:, i]), max(P_post[:, i])])
                    x_min, x_max = axs[i, j].get_xlim()
                    dx = x_max - x_min
                    xticks = np.array([x_min + 1/12*dx, x_min + 6/12*dx, x_min + 11/12*dx])
                    axs[i, j].set_xticks(xticks)  # Set x-tick positions first, then labels, otherwise warning
                    axs[i, j].set_xticklabels([str(round(tick,2)) for tick in xticks])
                    y_min, y_max = axs[i, j].get_ylim()
                    dy = y_max - y_min
                    yticks = np.array([y_min + 1/12*dy, y_min + 6/12*dy, y_min + 11/12*dy])
                    axs[i, j].set_yticks(yticks)  # Set y-tick positions first, then labels, otherwise warning
                    axs[i, j].set_yticklabels([str(round(tick,2)) for tick in yticks])
                    axs[i, j].minorticks_on()
                    axs[i, j].tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
                    axs[i, j].tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)
                    axs[i, j].tick_params(which='major', length = length_major_ticks, color = 'black')
                    axs[i, j].tick_params(which='minor', length = length_minor_ticks, color = 'black')
                    # Format the x and y axis with thousands separators
                    # axs[i, j].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
                    # axs[i, j].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
                    if i == DREAMPar['d'] - 1:
                        axs[i, j].set_xlabel(str_par[j], fontsize = fontsize_xylabel)
                    else:
                        axs[i, j].set_xticklabels([])
                    if j == 0:
                        axs[i, j].set_ylabel(str_par[i], fontsize = fontsize_xylabel)
                    else:
                        axs[i, j].set_yticklabels([])
            fig.suptitle(f"{method_fig}: Marginal distribution and bivariate scatter plots of the posterior samples", fontsize = fontsize_title)
            pdf.savefig()
            plt.show() # (block=False)

        ## --------------------------------------
        ## PLOT traceplots for parameters
        ## --------------------------------------
        symbol = ['bs', 'rx', 'g+', 'yo', 'm<']
        symbol_color = ['b','r','g','y','m']
        symbol_marker = ['s','x','+','o','.']
        counter = 0
        while counter <= DREAMPar['d'] - 1:
            fig, axs = plt.subplots(2, 1, figsize=(15, 8), zorder=1)
            plt.subplots_adjust(wspace=0.4, hspace=0.4)  # Increase the space horizontally and vertically
            for ax in axs.flat:
                for i in range(n_chains):
                    ax.plot(np.arange(0, T), chain[:, counter, i], symbol[i], markersize = markersize_symbols, markeredgewidth=2, linewidth = 2, label = str_chain[i], zorder=2)
                ax.set_xlim(0, T-1)                                                          # Set x-axis limits
                y_min, y_max = ax.get_ylim()
                ax.set_ylim(min(0.9 * y_min, 1.1 * y_min), max(0.9 * y_max, 1.1 * y_max))  # Adjust y limits with padding
                xticks = np.concatenate([[0], np.arange((T-1)//5, T, (T-1)//5)])
                ax.set_xticks(xticks)  # Set x-tick positions
                ax.set_xticklabels([str(int(DREAMPar['thinning'] * tick)) for tick in xticks])  # Set x-tick labels with thinning included
                ax.plot(T-1,MAP[counter], 'x', markersize = 12, color = 'gray', markeredgewidth=3, linewidth = 5, zorder=3, label = str_chain[n_chains], clip_on=False)
                ax.set_xlabel('Sample number of chain', fontsize = fontsize_xylabel, labelpad = 10)
                ax.set_ylabel(f'Parameter {str_par[counter]}', fontsize = fontsize_xylabel, labelpad = 10)
                ax.set_title(f"{method_fig}: Sampled traceplot of parameter {str_par[counter]} ", fontsize = fontsize_title)
                legend_handles = [Line2D([0], [0], color = symbol_color[i], lw = 0, marker=symbol_marker[i], markersize = markersize_legend, markeredgewidth=3, label = str_chain[i]) for i in range(n_chains)]
                legend_handles[n_chains:] = [Line2D([0], [0], color = 'gray', lw = 0, marker='x', markersize = markersize_legend, markeredgewidth=3, label = str_chain[n_chains])]
                legend = ax.legend(handles = legend_handles, fontsize = fontsize_legend, loc='upper right')
                symbols_ext = symbol_color[0:n_chains] + ['gray']
                for z, text in enumerate(legend.get_texts()):
                    text.set_color(symbols_ext[z])
                # ax.legend(str_chain, loc='best', fontsize = fontsize_legend)
                ax.tick_params(axis = 'both', labelsize = fontsize_axis)
                label_plot = get_label(counter)  # This converts 0 -> 'A', 1 -> 'B', etc.
                ax.text(0.01,0.97, f'({label_plot})', transform = ax.transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left')
                ax.minorticks_on()
                ax.tick_params(which='major', length = length_major_ticks, color = 'black')
                ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
                ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
                ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers)
                # Format the x and y axis with thousands separators
                ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
                # ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
                # Optionally, you can set major ticks locator to improve the appearance
                ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
                # ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
                counter += 1                
                if counter == DREAMPar['d']:
                    for ax in axs.flat:
                        if not ax.has_data():   # If the axis has no data, remove it
                            fig.delaxes(ax)     # delete empty subplots
                    break
            #plt.tight_layout()
            pdf.savefig()                
            plt.show() # (block=False)

        ## --------------------------------------
        ## PLOT traceplot of prior and log-likelihood function
        ## --------------------------------------
        fig, axs = plt.subplots(2, 1, figsize=(15, 8))
        plt.subplots_adjust(wspace=0.4, hspace=0.4)  # Increase the space horizontally and vertically
        for i in range(2):
            for j in range(n_chains):
                if i == 0:
                    axs[0].plot(np.arange(0, T), chain[:, DREAMPar['d'], j], symbol[j], markersize = markersize_symbols, markeredgewidth=3, linewidth = 2, label = str_chain[j], zorder=2)
                else:
                    axs[1].plot(np.arange(0, T), chain[:, DREAMPar['d'] + 1, j], symbol[j], markersize = markersize_symbols, markeredgewidth=3, linewidth = 2, label = str_chain[j], zorder=2)
            axs[i].set_xlim(0, T-1)                # Set x-axis limits
            y_min, y_max = axs[i].get_ylim()
            axs[i].set_ylim(min(0.9 * y_min, 1.1 * y_min), max(0.9 * y_max, 1.1 * y_max))  # Adjust y limits with padding
            xticks = np.concatenate([[0], np.arange((T-1)//5, T, (T-1)//5)])
            axs[i].set_xticks(xticks)  # Set x-tick positions
            axs[i].set_xticklabels([str(int(tick * DREAMPar['thinning'])) for tick in xticks])  # Set x-tick labels with thinning included
            axs[i].set_xlabel('Number of samples', fontsize = fontsize_xylabel, labelpad = 10)
            if i == 0:
                axs[i].set_ylabel(r'$\mathcal{P}(\mathbf{x})$', fontsize = 18, loc='center', labelpad = 10)
                axs[i].set_title(f"{method_fig}: Sampled traceplot of log-prior", fontsize = fontsize_title)
            else:
                if DREAMPar['lik'] == 23:   ## Limits of acceptability
                    axs[i].set_ylabel(f'Number of limits satisfied', fontsize = 18, loc='center', labelpad = 10)
                    axs[i].set_title(f"{method_fig}: Sampled traceplot of # of limits of acceptability satisfied", fontsize = fontsize_title)
                else:                       ## Other likelihood functions          
                    axs[i].set_ylabel(r'$\mathcal{L}(\mathbf{x} \vert \widetilde{\mathbf{y}})$', fontsize = 18, loc='center', labelpad = 10)
                    axs[i].set_title(f"{method_fig}: Sampled traceplot of log-likelihood", fontsize = fontsize_title)
            legend_handles = [Line2D([0], [0], color = symbol_color[j], lw = 0, marker=symbol_marker[j], markersize = markersize_legend, markeredgewidth=3, label = str_chain[j]) for j in range(n_chains)]
            legend = axs[i].legend(handles = legend_handles, fontsize = fontsize_legend, loc='upper right')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(symbols_ext[z])
            label_plot = get_label(i)  # This converts 0 -> 'A', 1 -> 'B', etc.
            axs[i].text(0.01,0.97, f'({label_plot})', transform = axs[i].transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left')
            axs[i].minorticks_on()
            axs[i].tick_params(which='major', length = length_major_ticks, color = 'black')
            axs[i].tick_params(which='minor', length = length_minor_ticks, color = 'black')
            axs[i].tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
            axs[i].tick_params(axis = 'y', labelsize = fontsize_axis_numbers)
            # Format the x and y axis with thousands separators
            axs[i].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            # axs[i].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            axs[i].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            # axs[i].yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
            
        pdf.savefig()
        plt.show() # (block=False)
        # plt.tight_layout()

        ## ---------------------------------------------------
        ## PLOT sample autocorrelation functions of parameters
        ## ---------------------------------------------------
        maxlag = min(200,T-1) 
        symbol = ['bs', 'rx', 'g+', 'yo', 'm<']
        symbol_color = ['b','r','g','y','m']
        symbol_marker = ['s','x','+','o','.']
        counter = 0
        while counter <= DREAMPar['d'] - 1:
            fig, axs = plt.subplots(2, 1, figsize = (15, 8), zorder=1)
            plt.subplots_adjust(wspace = 0.4, hspace = 0.4)  # Increase the space horizontally and vertically
            acf_chain = np.zeros((maxlag + 1, n_chains))
            for ax in axs.flat:
                # Plotting the ACF of each parameter
                for z in range(n_chains):
                    # Compute ACF for each parameter using the first n_chains
                    acf_chain[:, z] = acf(chain[:, counter, z], nlags = maxlag)
                    ax.plot(np.arange(0, maxlag + 1), acf_chain[:, z], color = symbol_color[z], linewidth = 2, label = str_chain[i], zorder=2)
                y_min = np.min(acf_chain)
                if y_min < 0:
                    y_min = 1.1 * y_min
                else:
                    y_min = 0.9 * y_min
                ax.set_xlim([0, maxlag])
                ax.set_ylim([y_min, 1.02])
                ax.tick_params(axis='both', which='both', labelsize = fontsize_axis_numbers, direction = 'out')
                xticks = np.concatenate([[0], np.arange((maxlag)//5, maxlag+1, (maxlag)//5)])
                ax.set_xticks(xticks)  # Set x-tick positions
                ax.set_xticklabels([str(int(tick)) for tick in xticks])
                # ax.set_xticklabels([str(int(DREAMPar['thinning'] * tick)) for tick in xticks])  # Set x-tick labels with thinning included
                ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
                ax.set_xlabel('Lag', fontsize = fontsize_xylabel, labelpad = 10)
                ax.set_ylabel('Correlation coefficient', fontsize = fontsize_xylabel, labelpad = 10)
                ax.set_title(f"{method_fig}: Sample autocorrelation function of parameter {str_par[counter]} ", fontsize = fontsize_title)
                legend_handles = [Line2D([0], [0], color = symbol_color[i], lw = 4, label = str_chain[i]) for i in range(n_chains)]
                legend = ax.legend(handles = legend_handles, fontsize = fontsize_legend, loc='upper right')
                for z, text in enumerate(legend.get_texts()):
                    text.set_color(symbols_ext[z])
                ax.tick_params(axis = 'both', labelsize = fontsize_axis)
                label_plot = get_label(counter)  # This converts 0 -> 'A', 1 -> 'B', etc.
                ax.text(0.01,0.97, f'({label_plot})', transform = ax.transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left')
                ax.minorticks_on()
                ax.tick_params(which='major', length = length_major_ticks, color = 'black')
                ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
                ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
                ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers)                           
                counter += 1                
                if counter == DREAMPar['d']:
                    for ax in axs.flat:
                        if not ax.has_data():   # If the axis has no data, remove it
                            fig.delaxes(ax)     # delete empty subplots
                    break
            #plt.tight_layout()
            pdf.savefig()                
            plt.show() # (block=False)

        ## --------------------------------------
        ## PLOT histograms of summary statistics
        ## --------------------------------------
        if S_post is not None:
            row, col = 2, 4
            idx_y_label = np.arange(0, Meas_info['n_S'] + 1 , col)
            Nbins = np.array([calcnbins(S_post[:, i]) for i in range(Meas_info['n_S'])])
            nbins = min(np.min(Nbins), maxbins)
            nbins = max(5, round(nbins / 2))  # Adjust number of bins
            counter = 0
            if np.isscalar(options['epsilon']):
                # vector with a single index
                epsln = options['epsilon'] * np.ones((1,Meas_info['n_S']))
            else:
                epsln = options['epsilon']
            epsln = epsln.reshape(-1,1)

            while counter <= Meas_info['n_S'] - 1:
                fig, axs = plt.subplots(row, col, figsize=(12, 8), zorder = 1)
                plt.subplots_adjust(wspace=0.4, hspace=0.4)  # Increase the space horizontally and vertically
                for ax in axs.flat:
                    U, edges = np.histogram(S_post[:, counter], bins=nbins, density=True)
                    X = 0.5 * (edges[:-1] + edges[1:])
                    ax.bar(X, U / np.max(U), width=(edges[1] - edges[0]), color = 'gray', alpha=0.7,zorder=2)
                    ax.set_xlabel(str_S[counter], fontsize = fontsize_xylabel, labelpad = 10)
                    x_min = ax.get_xlim()[0]
                    x_max = ax.get_xlim()[1]
                    # Plot vertical lines for S +/- epsilon values
                    x_minpos = S_meas[counter] - epsln[counter] # - dx/10
                    x_maxpos = S_meas[counter] + epsln[counter] # + dx/10
                    ax.plot([x_minpos, x_minpos], [0, 1.02], color = 'blue', linestyle='--', linewidth = 1)
                    ax.plot([x_maxpos, x_maxpos], [0, 1.02], color = 'blue', linestyle='--', linewidth = 1)
                    # Now add a vertical label next to it
                    str_left = r'$\widetilde{{S}}_{{{idx}}} - \epsilon_{{{idx}}}$'.format(idx = counter+1)
                    str_right = r'$\widetilde{{S}}_{{{idx}}} + \epsilon_{{{idx}}}$'.format(idx = counter+1)
                    dx = x_max - x_min                    
                    if x_min > x_minpos - dx/10: ## Vertical lines will not print well
                        x_min = x_minpos - dx
                    if x_max < x_maxpos + dx/10: ## Vertical lines will not print well    
                        x_max = x_maxpos + dx
                    # adjust axis
                    ax.set_xlim([x_min, x_max])
                    ax.text(x_minpos - dx/4, 0.3, str_left, rotation = 90, color = 'blue', ha = 'center', fontsize = fontsize_legend)
                    ax.text(x_maxpos + dx/3.2, 0.3, str_right, rotation = 90, color = 'blue', ha = 'center', fontsize = fontsize_legend)
                    yticks = np.arange(0, 1.02, 0.2)
                    ax.set_yticks(yticks)               # Set x-tick positions
                    ax.set_yticklabels([str(round(tick,1)) for tick in yticks])
                    ax.set_ylim(0, 1.02)                # Adjust y limits with padding
                    if counter in idx_y_label:
                        ax.set_ylabel('Empirical density', fontsize = fontsize_xylabel, labelpad = 10)
                    ax.plot(S_meas[counter], 1.02, 'bx', markersize = 12, markeredgewidth=3, linewidth = 3, zorder=3, clip_on=False)
                    label_plot = get_label(counter)  # This converts 0 -> 'A', 1 -> 'B', etc.
                    ax.text(0.02,0.97, f'({label_plot})', transform = ax.transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left')
                    ax.minorticks_on()
                    ax.tick_params(which='major', length = length_major_ticks, color = 'black')
                    ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
                    ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
                    ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)
                    # Format the x and y axis with thousands separators
                    # ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
                    # Optionally, you can set major ticks locator to improve the appearance
                    # ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
                    counter += 1
                    if counter == Meas_info['n_S']:
                        for ax in axs.flat:
                            if not ax.has_data():   # If the axis has no data, remove it
                                fig.delaxes(ax)     # delete empty subplots
                        break
                fig.suptitle(f"{method_fig}: Sampled marginal distribution of the summary statistics", fontsize = fontsize_title)
                #fig.tight_layout(rect=[0, 0, 1, 0.95])
                pdf.savefig()
                plt.show() # (block=False)        

        ## ------------------------------------------
        ## 5. Time series of summary statistics: LOA
        ## ------------------------------------------
        if S_post is not None and DREAMPar['lik'] == 23:
            # Limits of acceptability
            fig, ax = plt.subplots(figsize=(15, 6))
            # Find solutions that satisfy limits of acceptability
            idx_LOA = np.where(P_post[:, DREAMPar['d'] + 1] == Meas_info['n_S'])[0]       
            ax.set_facecolor('white')                 # Default white background color
            ax.set_prop_cycle(color = color_order)    # Color cycle
            if len(idx_LOA) == 0:       ## No "posterior" solutions that satisfy all limits of acceptability
                ax.plot(np.arange(1, Meas_info['n_S'] + 1), S_meas, 's', color = 'red', markerfacecolor = 'white', markersize = markersize_symbols, markeredgewidth=2, label = "Measured summary statistcs")
                ax.set_xlim([0, Meas_info['n_S']])
                #ax.set_ylim([min(0.9 * a[2], 1.1 * a[2]) for a in ax.axis()])
                ax.text(0.05, 0.90, 'No solution satisfies all limits of acceptability', fontsize = fontsize_xylabel, transform=ax.transAxes)
            else:                       ## solution that satisfies all limits
                # Derive minimum and maximum simulated values
                tot_S_unc = np.column_stack((np.max(S_post[idx_LOA, :Meas_info['n_S']], axis=0),
                                             np.min(S_post[idx_LOA, :Meas_info['n_S']], axis=0)))               
                # Fill the ranges (replace 'Fill_Ranges' with fill_between)
                y1 = tot_S_unc[:, 0]; y2 = tot_S_unc[:, 1]
                ax.fill_between(np.arange(1, Meas_info['n_S'] + 1), y1, y2, where=(y1 > y2), color = [0.7, 0.7, 0.7], alpha = 0.7, interpolate=True)
                # Now add the observations
                ax.plot(np.arange(1, Meas_info['n_S'] + 1), S_meas, 's', color = 'red', markerfacecolor = 'white', markersize = markersize_symbols, markeredgewidth=2, label = "Measured summary statistcs")
                ax.set_xlim([0, Meas_info['n_S']])
                # Add labels
                ax.set_xlabel('Observation number', fontsize = fontsize_xylabel)
                ax.set_ylabel('Data/simulated values', fontsize = fontsize_xylabel)
                ax.tick_params(axis = 'both', which = 'major', labelsize = fontsize_xylabel)
                # Calculate percentage inside the bounds of the behavioral simulation/prediction space
                coverage_LOA = 100 * np.sum((y2 < S_meas) & (S_meas < y1)) / Meas_info['n_S']
                ax.text(0.05, 0.90, f'Limits of acceptability envelope {coverage_LOA:5.2f}% of training data', fontsize = fontsize_xylabel, transform = ax.transAxes, ha = 'left')
            ax.set_title(f'{method_fig}: Limits of acceptability: Posterior simulations that satisfy all limits' + sndwch_text, fontsize = fontsize_title)
            label_plot = get_label(0)  # This converts 0 -> 'A', 1 -> 'B', etc.
            ax.text(0.01,0.98, f'({label_plot})', transform = ax.transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left')
            ax.minorticks_on()
            ax.tick_params(which='major', length = length_major_ticks, color = 'black')
            ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
            ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
            ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)
            # Format the x and y axis with thousands separators
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            # ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
            # ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
            pdf.savefig()
            plt.show() # (block=False)        

        ## --------------------------------------
        ## PLOT confidence and prediction intervals
        ## --------------------------------------
        if FX_post is not None:
            alfa1 = 100 * p_alfa / 2                    # Lower percentile
            alfa2 = 100 * (1 - p_alfa/2)                # Upper percentile
            # Set random seed for reproducibility
            np.random.seed(1 + round(100 * np.random.rand()))  # Equivalent to rng(1+round(100*rand),'twister')
            # Compute the p_alfa/2 & (1-p_alfa/2) confidence limits: parameter uncertainty
            par_unc = np.percentile(FX_post, [alfa1, alfa2], axis = 0).T
            if Meas_info['n'] == 1:
                par_unc = par_unc.T  # If only 1 measurement, transpose the result
                # Method for generating total uncertainty
                # [1] Take MAP parameter values and add structural uncertainty 
                #     --> parameter uncertainty explains part of total uncertainty
                # [2] Sample posterior parameters and add structural uncertainty on top
                #     --> does not guarantee that 95% contains 95%
            unc_method = 1
            # Method = 2: Does not work yet if modout = 'yes' -> need to add 
            # How: get "unique" samples of P_post
            if unc_method == 1:
                # Method 1: Using only MAP value
                MAP_un = X_unnormalize(MAP, Par_info)       # Transform to unnormalized space
                # Nr = FX_post.shape[0]                     # Number of replicates desired
                Nr = min(FX_post.shape[0],(1000))           # Maximum number of 1000 replicates      
                tot_unc, fx_mod = Bayes_pdf(MAP_un, FX_MAP.T, np.ones(Nr), Nr, RMSE_MAP, DREAMPar, Meas_info, Lik_info, p_gamma)
            elif unc_method == 2:
                # Method 2: Using posterior parameter samples
                nUP = len(iiUP)                             # Number of unique samples
                Nr = np.zeros(nUP, dtype=int)               # Vector with number of replicates
                for z in range(nUP):
                    Nr[z] = np.sum(idUP == z)               # How many replicates each row?
                # Get parameter and total uncertainty
                tot_unc, fx_mod = Bayes_pdf(UP_postun,FX_post[iiUP, :Meas_info['n']], idUP, Nr, RMSE_MAP, DREAMPar, Meas_info, Lik_info, p_gamma)
            # Compute spread of 100*p_gamma% confidence intervals for parameters
            width_par_unc = np.mean(par_unc[:, 1] - par_unc[:, 0])
            # Compute spread of 100*p_gamma% prediction intervals
            width_tot_unc = np.mean(tot_unc[:, 1] - tot_unc[:, 0])
            # Open new figure with normalized size and outer position
            fig, ax = plt.subplots(1,1, figsize=(16, 9), dpi=80)
            # Fill the ranges (replace 'Fill_Ranges' with fill_between)
            y1 = tot_unc[:, 1]; y2 = tot_unc[:, 0]
            ax.fill_between(np.arange(1, Meas_info['n'] + 1), y1, y2, where=(y1 > y2), color = [0.7, 0.7, 0.7], alpha = 0.7, interpolate=True)
            # Add gamma% uncertainty due to parameter uncertainty (dark gray)
            y1 = par_unc[:, 1]; y2 = par_unc[:, 0]
            ax.fill_between(np.arange(1, Meas_info['n'] + 1), y1, y2, where=(y1 > y2), color = [0.2, 0.2, 0.2], alpha = 0.7, interpolate=True)
            # Plot the verifying observations (red circles)
            ax.plot(np.arange(1, Meas_info['n'] + 1), Y_meas, 'ro', markersize = 6, linewidth = 1, markerfacecolor = 'w')
            # Plot the MAP simulation (black line)
            ax.plot(np.arange(1, Meas_info['n'] + 1), FX_MAP, 'k', linewidth = 1.5)
            ax.set_xlabel('Row (index) of training data set', fontsize = fontsize_xylabel, labelpad = 10)
            ax.set_ylabel('Simulated/measured values', fontsize = fontsize_xylabel, labelpad = 10)
            ax.minorticks_on()
            ax.tick_params(which='major', length = length_major_ticks, color = 'black')
            ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
            ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
            ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers)
            # Format the x and y axis with thousands separators
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            # ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            # ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
            # Create title
            fig.suptitle(f"{method_fig}: Snapshot of confidence (dark gray) and prediction (light gray) intervals", fontsize = fontsize_title)
            # Define legend location: get the y-axis limits
            min_X = 0
            max_X = min(500, Meas_info['n'])
            ax.set_xlim(0, max_X)       
            y_limits = ax.get_ylim()
            min_Y = y_limits[0]
            max_Y = y_limits[1]        
            x_loc = min_X + 0.05 * (max_X - min_X)
            y_loc = min_Y + 1 * (max_Y - min_Y)
            dx = 0.03 * (max_X - min_X)
            dy = 0.07 * (max_Y - min_Y)
            # Add axes to the figure manually (with [left, bottom, width, height] in normalized figure coordinates)
            values = int(100 * (1 - p_alfa))
            #pred_valstr = '/'.join(map(str, values)) + '%'
            pred_valstr = str(values) + '%'
            pred_str = 'Prediction interval'
            conf_valstr = str(values) + '%'
            conf_str = 'Confindence interval'
            # First entry (Prediction interval)
            x1 = x_loc
            x2 = x1 + dx
            # Add a rectangle patch to highlight a specific region
            rect = patches.Rectangle((x1, y_loc), dx, dy / 2, linewidth = 2, edgecolor = 'k', facecolor=[0.7, 0.7, 0.7])
            ax.add_patch(rect)
            # Add text for the prediction interval
            ax.text(x_loc + 1.4 * dx, y_loc + dy / 5, f"{pred_valstr} {pred_str}", fontsize = fontsize_text, color = 'k', verticalalignment = 'center', horizontalalignment = 'left')
            # Update the location for the next legend entry
            y_loc -= dy / 1.3
            # Second entry (Confidence interval)
            rect = patches.Rectangle((x_loc, y_loc), dx, dy/2 , linewidth = 2, edgecolor = 'k', facecolor=[0.2, 0.2, 0.2])
            ax.add_patch(rect)
            ax.text(x_loc + 1.4 * dx, y_loc + dy / 5, f"{conf_valstr} {conf_str}", fontsize = fontsize_text, color = [0.2, 0.2, 0.2], verticalalignment = 'center', horizontalalignment = 'left')
            # Third entry (MAP simulation)
            y_loc -= dy / 1.7
            ax.plot([x_loc, x_loc + dx], [y_loc, y_loc], color = 'k', linewidth = 4)
            ax.text(x_loc + 1.4 * dx, y_loc, f'MAP simulation', fontsize = fontsize_text, verticalalignment = 'center', horizontalalignment = 'left')
            # Last entry (Verifying observations)
            y_loc -= dy / 1.3
            ax.plot(x_loc + dx / 2, y_loc, 'ro', markersize = 8, linewidth = 3, markerfacecolor = 'w')
            ax.text(x_loc + 1.4 * dx, y_loc, 'Verifying observations', fontsize = fontsize_text, color = 'r', verticalalignment = 'center', horizontalalignment = 'left')
            # Add label with % contained (Coverage of the prediction intervals)
            x_loc = min_X + 0.6 * (max_X - min_X)
            y_loc = min_Y + 1.05 * (max_Y - min_Y)
            dx = 0.02 * (max_X - min_X)
            dy = 0.04 * (max_Y - min_Y)
            # Make a white patch for the coverage information
            x1 = x_loc - dx
            x2 = x_loc + 16 * dx
            rect = patches.Rectangle((x_loc, y_loc - 2 * dy), 15*dx, 2 * dy , linewidth = 1, edgecolor = 'k', facecolor='w')
            ax.add_patch(rect)
            # Update location for the coverage information
            y_loc = min_Y + 1.01 * (max_Y - min_Y)
            # Display coverage information for each prediction interval
            coverage_value = int(100 * p_gamma )  # Get the g value and convert it to a percentage
            coverage_percent = 100 * round( 100 * np.sum((tot_unc[:, 0].squeeze() < Y_meas) & (Y_meas < tot_unc[:, 1].squeeze())) / Meas_info['n']) / 100
            ax.text(x_loc + dx/2, y_loc, f'Coverage of {coverage_value}% pred interval is {coverage_percent}%', fontsize = fontsize_text, color = [0.7, 0.7, 0.7], verticalalignment = 'center', horizontalalignment = 'left')
            ax.tick_params(axis = 'both', labelsize = fontsize_axis)
            plt.tight_layout()
            pdf.savefig()
            plt.show() #  (block=False)


        ## --------------------------------------
        ## DO RESIDUAL ANALYSIS 
        ## --------------------------------------
        if FX_MAP is not None and DREAMPar['lik'] > 10 and Meas_info['n'] > 1:
            MAP_un = X_unnormalize(MAP, Par_info)
            # Compute normalized (partial) residuals and corresponding density --> MAP and FX_MAP need to be 2-dimensional
            _, e, __, X_n, XX_n = Calc_likelihood(np.array([MAP_un]), FX_MAP.reshape(-1,1), DREAMPar, Par_info, Meas_info, Lik_info, options, 5, MAP_info)

            # Residual computation based on the likelihood
            if DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:
                if Lik_info['t_start'] > 0:     ## Residual correlation treated
                    eps_n = X_n
                    f_eps_n = XX_n
                    par = Lik_info['fpar']                  # all parameters
                    par[Lik_info['id_vpar']] = MAP_un       # variable parameters
                    nuisvar = par[Lik_info['id_nuisvar']]   # isolate nuisance vars
                else:                           ## Residual correlation not treated
                    e_n = X_n
                    f_e_n = XX_n
            else: ## Likelihoods do not treat serial correlation
                e_n = X_n
                f_e_n = XX_n

            ## --------------------------------------
            # Plot A: Time series of observed and simulated data
            ## --------------------------------------
            fig, axs = plt.subplots(2, 1, figsize=(15, 10), zorder = 1)
            plt.subplots_adjust(wspace=0.3, hspace=0.3)  # Increase the space horizontally and vertically
            axs[0].plot(np.arange(1, Meas_info['n'] + 1), Y_meas, 's', color = 'red', markerfacecolor = 'white', markersize = markersize_symbols, markeredgewidth=2, label = "Measurement data")
            axs[0].plot(np.arange(1, Meas_info['n'] + 1), FX_MAP, color = 'black', linewidth = 1.5, label = "MAP simulation")
            axs[0].set_xlim([1, Meas_info['n']])
            axs[0].set_xlabel(f'Excerpt of training data set', fontsize = fontsize_xylabel, verticalalignment = 'top', horizontalalignment = 'center', labelpad = 10)
            axs[0].set_ylabel(f'Model output', fontsize = fontsize_xylabel, verticalalignment = 'top', horizontalalignment = 'center', labelpad = 20)
            legend_handles = [Line2D([0], [0], color = 'red', lw = 0, marker='s', markersize=markersize_legend, markerfacecolor = 'white', label = f'Measurement data')]                   
            legend_handles[1:] = [Line2D([0], [0], color = 'black', lw = 3, label = f'MAP simulation')]
            symbol_ext = ['red','black']
            legend = axs[0].legend(handles = legend_handles, fontsize = fontsize_legend, loc='upper right')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(symbol_ext[z])  
            # Add letter
            label_plot = get_label(0)  # This converts 0 -> 'A', 1 -> 'B', etc.
            axs[0].text(0.002,0.98, f'({label_plot})', transform = axs[0].transAxes, fontsize = fontsize_text, horizontalalignment ='left', va = 'top', ha = 'left', color = 'black')
            axs[0].minorticks_on()
            axs[0].tick_params(which='major', length = length_major_ticks, color = 'black')
            axs[0].tick_params(which='minor', length = length_minor_ticks, color = 'black')
            axs[0].tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
            axs[0].tick_params(axis = 'y', labelsize = fontsize_axis_numbers)
            # Format the x and y axis with thousands separators
            axs[0].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            # axs[0].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            axs[0].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            # axs[0].yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            

            ## ---------------------------------------------------
            # B: Now plot raw and/or studentized partial residuals
            ## ---------------------------------------------------
            axs[1].plot(np.arange(1, Meas_info['n'] + 1), e, color = dark_gray, linestyle = '-', linewidth = 1.5, label = 'Raw residuals')
            axs[1].set_ylabel('Raw residuals', fontsize = fontsize_xylabel, color = dark_gray)
            axs[1].tick_params(axis='y', labelcolor = dark_gray)
            axs[1].set_xlim([1, Meas_info['n']])
            axs[1].set_xticks(np.arange(0, Meas_info['n'], Meas_info['n'] // 5))
            axs[1].xaxis.set_minor_locator(plt.MaxNLocator(integer=True))
            axs[1].set_xlabel("Observation number", fontsize = fontsize_xylabel)
            legend_handles = [Line2D([0], [0], color = dark_gray, lw = 3, label = f'Raw residuals')]
            axs12 = axs[1].twinx()
            # Now check whether we usse autocorrelation treatment
            if Lik_info['t_start'] > 0:
                # Plot studentized partial residuals on right axis
                axs12.plot(np.arange(t_start+1, Meas_info['n'] + 1), eps_n[t_start:Meas_info['n']], 's', color = light_gray, markersize = 3, linewidth = 1.5, markerfacecolor = 'white', label = 'Studentized partial residuals')
                axs12.set_ylabel('Studentized partial residuals', fontsize = fontsize_xylabel, color = light_gray)
                legend_handles[1:] = [Line2D([0], [0], color = light_gray, lw = 0, marker='s', markersize=markersize_legend, markerfacecolor = 'white', label = f'Studentized partial residuals')]                   
            else:
                # Plot studentized raw residuals on right axis
                axs12.plot(np.arange(1, Meas_info['n'] + 1), e_n[0:Meas_info['n']], 's', color = light_gray, markersize = 3, linewidth = 1.5, markerfacecolor = 'white', label = 'Studentized raw residuals')
                axs12.set_ylabel('Studentized raw residuals', fontsize = fontsize_xylabel, color = light_gray)
                legend_handles[1:] = [Line2D([0], [0], color = light_gray, lw = 0, marker = 's', markersize = markersize_legend, markerfacecolor = 'white', label = f'Studentized raw residuals')]                   
            # Right axis (for studentized residuals)
            axs12.tick_params(axis='y', labelcolor = light_gray)
            axs12.set_ylim([-3 ,3])
            axs[1].legend(fontsize = fontsize_legend, loc='upper right')
            axs[1].tick_params(axis = 'both', direction = 'out', length = 5, width=1, colors = 'black')
            symbol_ext = [dark_gray,light_gray]
            legend = axs[1].legend(handles = legend_handles, fontsize = fontsize_legend, loc='upper right')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(symbol_ext[z])  
            # Add letter
            label_plot = get_label(1)  # This converts 0 -> 'A', 1 -> 'B', etc.
            axs[1].text(0.002,0.98, f'({label_plot})', transform = axs[1].transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left', color = 'black')
            axs[1].minorticks_on()
            axs[1].tick_params(which='major', length = length_major_ticks, color = 'black')
            axs[1].tick_params(which='minor', length = length_minor_ticks, color = 'black')
            axs[1].tick_params(axis = 'x', labelsize = fontsize_axis_numbers)
            axs[1].tick_params(axis = 'y', labelsize = fontsize_axis_numbers)
            # Format the x and y axis with thousands separators
            axs[1].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            # axs[1].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f'{int(y):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            axs[1].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            # axs[1].yaxis.set_major_locator(ticker.MaxNLocator(integer=True))            
            pdf.savefig()
            plt.show() #  (block=False)
            
            ## ------------------------------------------------------------------------------------------------------------
            ## Plot C: Variance of standardized partial [with AR treatment] or standardized raw residuals [no AR treatment]
            ## ------------------------------------------------------------------------------------------------------------
            fig, axs = plt.subplots(2, 2, figsize=(12, 10), zorder = 1)
            plt.subplots_adjust(wspace=0.35, hspace=0.35)  # Increase the space horizontally and vertically
            # For the first case
            if Lik_info['t_start'] > 0: # DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:
                axs[0, 0].plot(FX_MAP[t_start:Meas_info['n']], eps_n[t_start:Meas_info['n']], 's', color = light_gray, markersize = 5, linewidth = 1.5, markerfacecolor = 'white', label = 'Studentized partial residuals')
                # check if there is variation in FX_MAP, otherwise do not do regression
                if np.std(FX_MAP[t_start:Meas_info['n']]) > 0:
                    c = np.polyfit(FX_MAP[t_start:Meas_info['n']], eps_n[t_start:Meas_info['n']], 1)
                axs[0, 0].set_ylabel('Studentized partial residuals', fontsize = fontsize_xylabel)
                legend_handles = [Line2D([0], [0], color = light_gray, lw = 0, marker='s', markersize = markersize_legend, markerfacecolor = 'white', label = 'Studentized partial residuals')]
                axs[0, 0].set_ylim([np.min(eps_n[t_start:Meas_info['n']]), np.max(eps_n[t_start:Meas_info['n']])])

            else:  # For the second case
                axs[0, 0].plot(FX_MAP[t_start:Meas_info['n']], e_n[t_start:Meas_info['n']], 's', color = light_gray, markersize = 5, linewidth = 1.5, markerfacecolor = light_gray, label = 'Studentized raw residuals')
                if np.std(FX_MAP[t_start:Meas_info['n']]) > 0:
                    c = np.polyfit(FX_MAP[t_start:Meas_info['n']], e_n[t_start:Meas_info['n']], 1)
                axs[0, 0].set_ylabel('Studentized raw residuals', fontsize = fontsize_xylabel)
                legend_handles = [Line2D([0], [0], color = light_gray, lw = 0, marker='s', markersize =markersize_legend, markerfacecolor = 'white', label = 'Studentized raw residuals')]
                axs[0, 0].set_ylim([np.min(e_n[t_start:Meas_info['n']]), np.max(e_n[t_start:Meas_info['n']])])

            axs[0, 0].tick_params(axis = 'both', direction = 'out', length = 6)
            axs[0, 0].set_title(f'Stable variance? {sndwch_text}', fontsize = fontsize_title_small_graphs)
            axs[0, 0].set_xlabel('Model output', fontsize = fontsize_xylabel)
            axs[0, 0].set_xlim([np.min(FX_MAP[t_start:Meas_info['n']]), np.max(FX_MAP[t_start:Meas_info['n']])])
            if np.std(FX_MAP[t_start:Meas_info['n']]) > 0:
                x_vals = np.array(axs[0, 0].get_xlim())
                y_vals = c[0] * x_vals + c[1]
                axs[0, 0].plot(x_vals, y_vals, 'r', linewidth = 1.5, linestyle='--', label = 'Linear regression')
                legend_handles[1:] = [Line2D([0], [0], color = 'red', lw = 3, label = 'Linear regression')]
                symbol_ext = [light_gray, red_color]
            else:
                symbol_ext = [light_gray]

            legend = axs[0, 0].legend(handles = legend_handles, fontsize = fontsize_legend, loc='upper right')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(symbol_ext[z])  
            # Add letter
            label_plot = get_label(2)  # This converts 0 -> 'A', 1 -> 'B', etc.
            axs[0, 0].text(0.01,0.98, f'({label_plot})', transform = axs[0, 0].transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left', color = 'black')
            axs[0, 0].minorticks_on()
            axs[0, 0].tick_params(which='major', length = length_major_ticks, color = 'black')
            axs[0, 0].tick_params(which='minor', length = length_minor_ticks, color = 'black')
            axs[0, 0].tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
            axs[0, 0].tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)

            ## --------------------------------------
            # D: Histogram of residuals
            ## --------------------------------------
            nl = min(25, Meas_info['n'] // 20)
            if Lik_info['t_start'] > 0: # DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:
                res_plot = eps_n[t_start:Meas_info['n']]
            else:
                res_plot = e_n[t_start:Meas_info['n']]
            # Determine number of bins
            if nl >= 1:
                nbins = int(min(calcnbins(res_plot), 2 * maxbins))
                n_res, edges = np.histogram(res_plot, bins = nbins, density = True)
                # Adjust edges and n_res for plotting
                n_res = np.insert(n_res, 0, 0)
                n_res = np.append(n_res, 0)
                edges = np.insert(edges, 0, -np.inf)
                edges = np.append(edges, np.inf)
                # Calculate midpoints for the bins
                x_res = 0.5 * (edges[:-1] + edges[1:])
                # Plot marginal PDF of residuals using stem plot
                # print(x_res,n_res)
                if Lik_info['t_start'] > 0:
                    axs[0, 1].stem(x_res, n_res, linefmt = light_gray, markerfmt = 's', basefmt = " ", label = 'Studentized partial residuals')
                    axs[0, 1].set_xlabel('Studentized partial residuals', fontsize = fontsize_xylabel)
                    legend_handles = [Line2D([0], [0], color = light_gray, lw = 0, marker = 's', markersize = markersize_legend, markerfacecolor = 'white', label = 'Studentized partial residuals')]
                else:
                    axs[0, 1].stem(x_res, n_res, linefmt = light_gray, markerfmt = 's', basefmt = " ", label = 'Studentized raw residuals')
                    axs[0, 1].set_xlabel('Studentized raw residuals', fontsize = fontsize_xylabel)                   
                    legend_handles = [Line2D([0], [0], color = light_gray, lw = 0, marker = 's', markersize = markersize_legend, markerfacecolor = 'white', label = 'Studentized raw residuals')]
                axs[0, 1].set_ylabel('Density', fontsize = fontsize_xylabel)
                axs[0, 1].set_title(f'Histograms {sndwch_text}', fontsize = fontsize_title_small_graphs)
                axs[0, 1].tick_params(axis = 'both', direction = 'out', length = 5, width=1, colors = 'black')
                r_99 = percentile(res_plot, [0.5, 99.5])
                axs[0, 1].set_xlim([r_99[0], r_99[1]])
                axs[0, 1].set_xticks(np.arange(int(r_99[0]), int(r_99[1])+1, 1))
            # Add likelihood plot (line of likelihood)
            if Lik_info['t_start'] > 0:
                data_res = np.column_stack([eps_n[t_start:Meas_info['n']], f_eps_n[t_start:Meas_info['n']]])
                # Sort the residual data and plot
                data_res = data_res[data_res[:, 0].argsort()]
                axs[0, 1].plot(data_res[:, 0], data_res[:, 1], color = 'r', linewidth = 2, label = 'Studentized partial residuals')
            else:
                data_res = np.column_stack([e_n[t_start:Meas_info['n']], f_e_n[t_start:Meas_info['n']]])
                # Sort the residual data and plot
                data_res = data_res[data_res[:, 0].argsort()]
                axs[0, 1].plot(data_res[:, 0], data_res[:, 1], color = 'r', linewidth = 2, label = 'Studentized raw residuals')
            # determine a suitable maximum maximum y-value
            y_max = 1.1 * max(np.max(data_res[:,1]), np.max(n_res))
            axs[0, 1].set_ylim([0,y_max * 1.1])
            axs[0, 1].yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1f}'))
            lik_number = DREAMPar['lik']
            legend_handles[1:] = [Line2D([0], [0], color = 'red', lw = 3, label = f'Likelihood function: {lik_number}')]
            legend = axs[0, 1].legend(handles = legend_handles, fontsize = fontsize_legend, loc = 'upper right')
            symbol_ext = [light_gray, red_color]
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(symbol_ext[z])  
            # Add letter
            label_plot = get_label(3)  # This converts 0 -> 'A', 1 -> 'B', etc.
            axs[0, 1].text(0.01,0.98, f'({label_plot})', transform = axs[0, 1].transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left', color = 'black')
            axs[0, 1].minorticks_on()
            axs[0, 1].tick_params(which='major', length = length_major_ticks, color = 'black')
            axs[0, 1].tick_params(which='minor', length = length_minor_ticks, color = 'black')
            axs[0, 1].tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
            axs[0, 1].tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)
            # Format the x and y axis with thousands separators
            # axs[0, 1].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'{int(x):,}'))
            # Optionally, you can set major ticks locator to improve the appearance
            # axs[0, 1].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        
            ## -------------------------------------------
            # E: Autocorrelation function of the residuals
            ## -------------------------------------------
            numLags = min(50, Meas_info['n'] - 1)
            lags = np.arange(numLags + 1)
            # Calculate ACF based on likelihood model
            if Lik_info['t_start'] > 0: # DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:
                acfunc, bnds = acf(eps_n[t_start:Meas_info['n']], nlags = numLags, alpha = 0.05)
                legend_handles = [Line2D([0], [0], color = light_gray, lw = 0, marker = 's', markersize = markersize_legend, markerfacecolor = 'white', label = 'Studentized partial residuals')]
            else:
                acfunc, bnds = acf(e_n[t_start:Meas_info['n']], nlags = numLags, alpha = 0.05)
                legend_handles = [Line2D([0], [0], color = light_gray, lw = 0, marker = 's', markersize = markersize_legend, markerfacecolor = 'white', label = 'Studentized raw residuals')]
            # We do not use bnds [= exact confidence intervals at each lag] but rather compute the theoretical confidence interval for large n
            bnds[0] = - 2/np.sqrt(Meas_info['n'])
            bnds[1] = - bnds[0]
            axs[1, 0].stem(lags, acfunc, markerfmt = 's', basefmt = " ", linefmt = light_gray)
            axs[1, 0].set_xlabel('Lag', fontsize = fontsize_xylabel)
            axs[1, 0].set_ylabel('Autocorrelation', fontsize = fontsize_xylabel)
            axs[1, 0].set_title(f'Sample ACF {sndwch_text}', fontsize = fontsize_title_small_graphs)
            axs[1, 0].tick_params(axis = 'both', direction = 'out', length = 5, width = 1, colors = 'black')
            axs[1, 0].set_xlim([0, numLags])
            # Plot bounds
            axs[1, 0].plot([0.5, numLags], [bnds[0], bnds[0]], 'r--', linewidth = 1)
            axs[1, 0].plot([0.5, numLags], [bnds[1], bnds[1]], 'r--', linewidth = 1)
            axs[1, 0].set_yticks(np.arange(-1, 1.1, 0.2))
            axs[1, 0].set_yticklabels(['-1.0', '-0.8', '-0.6', '-0.4', '-0.2', '0.0', '0.2', '0.4', '0.6', '0.8', '1.0'])
            axs[1, 0].set_ylim([min(-0.2, 1.1 * np.min(acfunc)), max(1, 1.1 * np.max(acfunc))])
            legend_handles[1:] = [Line2D([0], [0], color = 'red', lw = 3, label = '95% confidence limits')]
            symbol_ext = [light_gray, red_color]
            legend = axs[1, 0].legend(handles = legend_handles, fontsize = fontsize_legend, loc='upper right')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(symbol_ext[z])
            # Add letter
            label_plot = get_label(4)  # This converts 0 -> 'A', 1 -> 'B', etc.
            axs[1, 0].text(0.01,0.98, f'({label_plot})', transform = axs[1, 0].transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left', color = 'black')
            axs[1, 0].minorticks_on()
            axs[1, 0].tick_params(which='major', length = length_major_ticks, color = 'black')
            axs[1, 0].tick_params(which='minor', length = length_minor_ticks, color = 'black')
            axs[1, 0].tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
            axs[1, 0].tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)

            ## -------------------------------------
            # G: Quantile-Quantile plot of residuals
            ## -------------------------------------
            if Lik_info['t_start'] > 0: # DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:
                y = np.sort(eps_n[t_start:Meas_info['n']])
            else:
                y = np.sort(e_n[t_start:Meas_info['n']])
            # Compute plotting position
            pp = np.linspace(1 / (len(y) + 1), len(y) / (len(y) + 1), len(y))
            # Determine theoretical quantiles based on the distribution
            if DREAMPar['lik'] in [14, 44]:
                x = SEPinv(pp, nuisvar[2], nuisvar[3])
                str_leg = 'SEP'
                x_lab = 'SEP quantiles'
            elif DREAMPar['lik'] == 16:
                x = LAPinv(pp, 0, 1)
                str_leg = 'Laplace'
                x_lab = 'Laplace quantiles'
            elif DREAMPar['lik'] == 17:
                x = SSTinv(pp, nuisvar[2], nuisvar[3])
                str_leg = 'SST'
                x_lab = 'SST quantiles'
            elif DREAMPar['lik'] == 45:
                x = SGTinv(pp, 0, 1, nuisvar[2], nuisvar[3], nuisvar[4])
                str_leg = 'SGT'
                x_lab = 'SGT quantiles'
            else:
                x = norm.ppf(pp)
                str_leg = 'Normal'
                x_lab = 'Normal quantiles'
            # Plot sorted residuals vs. theoretical quantiles
            axs[1, 1].scatter(x, y, marker = 's', color = light_gray, s = 8, linewidth = 1.5, facecolors='white')
            axs[1, 1].set_title(f'Quantile-quantile plot {sndwch_text}', fontsize = fontsize_title_small_graphs)
            axs[1, 1].plot(x, x, 'r--', linewidth = 1.5, label = '1:1 relationship')
            axs[1, 1].set_xlabel(x_lab, fontsize = fontsize_xylabel)
            if Lik_info['t_start'] > 0: 
                axs[1,1].set_ylabel(f'Quantiles of st. partial residuals', fontsize = fontsize_xylabel)
            else:
                axs[1,1].set_ylabel(f'Quantiles of st. raw residuals', fontsize = fontsize_xylabel)
            # Customize ticks
            axs[1, 1].tick_params(axis = 'both', direction = 'out', length = 5, width=1, colors = 'black')
            axs[1, 1].set_xlim([np.min(x), np.max(x)])
            axs[1, 1].set_ylim([np.min(y), np.max(y)])
            legend_handles is None
            legend_handles = [Line2D([0], [0], color = 'red', lw = 3, label = '1:1 relationship')]
            symbol_ext = ['red']
            legend = axs[1, 1].legend(handles = legend_handles, fontsize = fontsize_legend, loc='upper right')
            # Match text color to line color
            for z, text in enumerate(legend.get_texts()):
                text.set_color(symbol_ext[z])
            # Add letter
            label_plot = get_label(5)  # This converts 0 -> 'A', 1 -> 'B', etc.
            axs[1, 1].text(0.01,0.98, f'({label_plot})', transform = axs[1, 1].transAxes, fontsize = fontsize_text, horizontalalignment = 'left', va = 'top', ha = 'left', color = 'black')
            axs[1, 1].minorticks_on()
            axs[1, 1].tick_params(which='major', length = length_major_ticks, color = 'black')
            axs[1, 1].tick_params(which='minor', length = length_minor_ticks, color = 'black')
            axs[1, 1].tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
            axs[1, 1].tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)
            pdf.savefig()
            plt.show() # (block=False)
  
            # Compute scoring rules after Vrugt, 2024
            fx_mod = fx_mod.T
            FX_post = FX_post.T
            # Compute scoring rules of forecast density
            LS, _, num_zeroLS = log_score(fx_mod, Y_meas)                       # Logarithmic score
            CRPS, _, __ = CRP_score(fx_mod, Y_meas)                             # Continuous Ranked Probability Score
            SS, _, num_zero = spherical_score(fx_mod, Y_meas)                   # Spherical score (zeta = 2)
            DSS, _, _, m_F, s_F = dawid_sebas_score(fx_mod, Y_meas)             # Dawid-Sebastiani score
            IS, _ = interval_score(tot_unc, Y_meas, p_alfa)                     # Interval score
            # Compute summary metrics of forecast density
            p_val, _ = p_values(fx_mod, Y_meas)                                 # Compute p-values
            RLBL, eCDF, uCDF = rlbl(p_val)                                      # Reliability from p-values
            CV = (1 / Meas_info['n']) * np.sum(s_F / m_F)                       # Coefficient of variation
                                                                                # CHECK: dimensions are fine but CV differs from MATLAB
            width_par = np.mean(par_unc[:, 1] - par_unc[:, 0])                  # Mean spread/width of confidence interval
            width_mod = np.mean(tot_unc[:, 1] - tot_unc[:, 0])                  # Mean spread/width of prediction interval
            prec_par = np.mean(np.std(FX_post, axis = 1)) / np.mean(Y_meas)     # Mean precision of parameters
            prec_mod = np.mean(np.std(fx_mod, axis = 1)) / np.mean(Y_meas)      # Mean precision of total
            pbias_par = 100 * np.sum(np.mean(FX_post, axis = 1) - Y_meas) / np.sum(Y_meas)      # Mean percentage bias
            pbias_mod = 100 * np.sum(np.mean(fx_mod, axis = 1) - Y_meas) / np.sum(Y_meas)
            meanLogP = np.mean(np.sum(P_post[:, DREAMPar['d']:DREAMPar['d']+2], axis = 1))      # Mean log-density of posterior realizations
            stdLogP = np.std(np.sum(P_post[:, DREAMPar['d']:DREAMPar['d']+2], axis = 1))        # Std deviation log-density of posterior realizations
            # Compute summary metrics of MAP solution
            maxLogL = np.max(P_post[:, DREAMPar['d']+1])                                        # Max of log-likelihood
            maxLogP = np.max(np.sum(P_post[:, DREAMPar['d']:DREAMPar['d']+2], axis = 1))        # Max of log posterior density

            ## --------------------------------------
            # A. Quantile - quantile plot of empirical and theoretical CDF of p-values
            ## --------------------------------------
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.plot(eCDF, uCDF, color = light_gray, label = str_leg)
            ax.plot([0, 1], [0, 1], 'r--', linewidth = 1, label = "1:1 line")
            ax.set_xlabel('Theoretical quantiles', fontsize = fontsize_xylabel)
            ax.set_ylabel('Empirical quantiles', fontsize = fontsize_xylabel)
            ax.set_title(f'Quantile-quantile plot empirical & theoretical CDF of p-values {sndwch_text}', fontsize = fontsize_title)
            ax.legend()
            ax.minorticks_on()
            ax.tick_params(which='major', length = length_major_ticks, color = 'black')
            ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
            ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
            ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)
            pdf.savefig()
            plt.show() # (block=False)

            ## --------------------------------------
            # B. Rank Histogram
            ## --------------------------------------
            rank_hist = np.ceil(p_val * fx_mod.shape[1]) + 1   # Compute rank histogram
            # nr, xr = np.histogram(rank_hist, bins = np.arange(1, fx_mod.shape[1]+2), density = True)
            nr, edges = np.histogram(rank_hist, bins = maxbins) #, density = False)
            # Bin centers and normalize with upper bin value to get Xr in [0,1]
            Xr = 0.5 * (edges[:-1] + edges[1:])
            Xr = Xr / np.max(Xr)
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.bar(Xr, nr / np.trapezoid(nr, Xr))
            ax.set_xlabel('Scaled rank of observation', fontsize = fontsize_xylabel)
            ax.set_ylabel('Probability density', fontsize = fontsize_xylabel)
            ax.tick_params(axis='both', labelsize = fontsize_axis_numbers, direction = 'out', length = 6, width = 1)
            ax.set_title(f'Rank histogram of the observations given PDF of total uncertainty {sndwch_text}', fontsize = fontsize_title)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))
            ax.minorticks_on()
            ax.tick_params(which='major', length = length_major_ticks, color = 'black')
            ax.tick_params(which='minor', length = length_minor_ticks, color = 'black')
            ax.tick_params(axis = 'x', labelsize = fontsize_axis_numbers_small_graphs)
            ax.tick_params(axis = 'y', labelsize = fontsize_axis_numbers_small_graphs)
            pdf.savefig()
            plt.show() # (block=False)

    # Open the final PDFs
    os.startfile(file_name)

    if FX_MAP is not None:
        # Print results
        if sndwch == 0:
            file_metric = f'{method_table}_forecast_metrics.txt'
            file_figs = f'{method_table}_figures.pdf'
            file_data = f'{method_table}_data'
        else:
            file_metric = f'{method_table}_forecast_metrics_sandwich.txt'
            file_figs = f'{method_table}_figures_sandwich.pdf'
            file_data = f'{method_table}_data_sandwich'

        if 'filename' in Lik_info:
            if isinstance(func_handle, str):
                # func_handle is a string [so that function is pickable for parallellizaton with IO = 'yes']
                name_mod = str(func_handle)
            else:
                # func_handle is a dynamic function [function not pickable but not a problem for parallellization if IO = 'no']
                name_handle = func_handle.__name__
                name_mod = str(name_handle)

            str_prt_nuis = ''
            if Lik_info['str_nuis'] is not None:
                for z in range(len(Lik_info['str_nuis'])):
                    str_prt = Lik_info['str_nuis'][z]
                    str_prt_nuis += f'_{str_prt[1:-1]}'

            str_prt_nuis = str_prt_nuis.replace('\\', '')
            name_lik_func = Lik_info['name_lik_func']
            file_name = f'{name_mod}_{name_lik_func}{str_prt_nuis}'
            # Define filenames for metrics, figures, and data files
            if sndwch == 0:
                file_metric = f'{file_name}_forecast_metrics_normalized.txt'
                file_figs = f'{file_name}_figures_normalized.pdf'
                file_data = f'{file_name}_data_normalized.mat'
            else:
                file_metric = f'{file_name}_forecast_metrics_normalized_sandwich.txt'
                file_figs = f'{file_name}_figures_normalized_sandwich.pdf'
                file_data = f'{file_name}_data_normalized_sandwich.mat'

        # Print Table to file
        with open(file_metric, 'w', encoding='utf-8') as fid_latex:
            fid_latex.write("="*103 + "\n")
            fid_latex.write(f"Likelihood {DREAMPar['lik']}\n")           
            if sndwch == 0:
                fid_latex.write("Sandwich correction is INACTIVE\n")
            elif sndwch == 1:
                fid_latex.write("Sandwich correction is ACTIVE\n")            
            fid_latex.write(f"{Lik_info['name_lik']}\n")            
            if Lik_info['str_nuis']:
                str_pr = str(Lik_info['str_nuis'][0])
                fid_latex.write(f" Nuisance variables: {str_pr[1:-1]}")
                for z in range(1, len(Lik_info['str_nuis'])):
                    str_pr = str(Lik_info['str_nuis'][z])
                    fid_latex.write(f", {str_pr[1:-1]}")
                fid_latex.write("\n")
            
            fid_latex.write(f" Total number of observations/simulated output: {Meas_info['n']}\n")
            fid_latex.write("="*103 + "\n")
            fid_latex.write(f"Score rules                           Abbrev.              Value   Unit        Reference                     \n")
            fid_latex.write("="*103 + "\n")
            fid_latex.write(f"Logarithmic Score                     LS              {LS:10.3f}   [log(1/y)]  Good 1952                     \n")
            fid_latex.write(f"Continuous Ranked Probability Score   CRPS            {CRPS:10.3f}   [y]         Matheson & Winkler 1972     \n")
            fid_latex.write(f"Spherical Score                       SS              {SS:10.3f}   [-]         Good 1971, Friedman 1983      \n")
            fid_latex.write(f"Dawid-Sebastiani Score                DSS             {DSS:10.3f}   [log(y^2)]  Dawid & Sebastiani 1999      \n")
            fid_latex.write(f"Interval Score                        IS              {IS:10.3f}   [y]         Gneiting & Raftery 2007       \n")
            fid_latex.write("="*103 + "\n")
            fid_latex.write(f"Summary metrics density forecast      Abbrev.              Value   Unit        Reference                     \n")
            fid_latex.write("="*103 + "\n")
            fid_latex.write(f"Reliability                           RLBL            {RLBL:10.3f}   [-]         Renard et al. 2011          \n")
            fid_latex.write(f"Coefficient of variation              CV              {CV:10.3f}   [-]         Evin et al. 2013              \n")
            fid_latex.write(f"Containment of 95% tot interval       ctnd_tot        {coverage_percent:10.3f}   [%]                         \n")
            fid_latex.write(f"Mean tot. width (spread at 95%)       wdth_mod        {width_mod:10.3f}   [y] ♫                              \n")
            fid_latex.write(f"Mean par. width (spread at 95%)       wdth_par        {width_par:10.3f}   [y] ♫                              \n")
            fid_latex.write(f"Mean par. precision                   prec_par        {prec_par:10.3f}   [-]         McInerney et al. 2017   \n")
            fid_latex.write(f"Mean tot. precision                   prec_mod        {prec_mod:10.3f}   [-]         McInerney et al. 2017   \n")
            fid_latex.write(f"Mean percentage bias par              pbias_par       {pbias_par:10.3f}   [%]         McInerney et al. 2017  \n")
            fid_latex.write(f"Mean percentage bias tot              pbias_mod       {pbias_mod:10.3f}   [%]         McInerney et al. 2017  \n")
            fid_latex.write(f"Mean log-posterior density            meanLogP        {meanLogP:10.3f}   [-]                                 \n")
            fid_latex.write(f"St. dev. of log-posterior density     stdLogP         {stdLogP:10.3f}   [-]                                  \n")
            fid_latex.write("="*103 + "\n")
            fid_latex.write(f"Summary metrics of MAP solution       Abbrev.              Value   Unit        Reference                     \n")
            fid_latex.write("="*103 + "\n")
            fid_latex.write(f"Maximum log-posterior density         maxLogP         {maxLogP:10.3f}   [-]                                  \n")
            fid_latex.write(f"Maximum log-likelihood                maxLogL         {maxLogL:10.3f}   [-]                                  \n")
            fid_latex.write(f"Root mean square error MAP output     RMSE_MAP        {RMSE_MAP:10.3f}   [y] ♫                               \n")
            if 'PBIAS_MAP' not in locals():
                PBIAS_MAP = np.nan
            fid_latex.write(f"Percentage bias MAP output            PBIAS_MAP       {PBIAS_MAP:10.3f}   [%]                                \n")
            fid_latex.write("="*103 + "\n")
            fid_latex.write(f"♫ Of the {Meas_info['n']} observations {num_zeroLS} have density less than real_min (= 2.22e-308) per forecast pdf \n")
            fid_latex.write(f"♫ Unit of observations in training data record \n")
        
        # Open the file in a text editor (for MATLAB, in Python we could use the default system editor)
        os.system(f'notepad {file_metric}')             # This will work on Windows
        
        # Set the maximum column width to be wider
        pd.set_option('display.width', 120)       
        # Print Table
        metrics = [
            'Logarithmic Score', 'Continuous Ranked Probability Score', 'Spherical Score',
            'Dawid-Sebastiani Score', 'Interval Score', "", 'Reliability', 'Coefficient of variation',
            'Percentage contained in 95% total uncertainty', 'Width of 95% total uncertainty', 
            'Width of 95% parameter uncertainty', 'Precision of 95% parameter uncertainty',
            'Precision of 95% total uncertainty', 'Percentage bias of 95% parameter uncertainty',
            'Percentage bias of 95% total uncertainty', 'Mean log-posterior density', 
            'Std. dev. of log-posterior density', "", 'Maximum log-posterior density', 
            'Maximum log-likelihood', 'Root Mean Square Error of MAP solution', 'Percentage bias of MAP solution'
        ]
        values = [
            LS, CRPS, SS, DSS, IS, "", RLBL, CV, coverage_percent, width_mod, width_par, 
            prec_par, prec_mod, pbias_par, pbias_mod, meanLogP, stdLogP, "", maxLogP, 
            maxLogL, RMSE_MAP, PBIAS_MAP
        ]
        # Convert metrics to dataframe for output
        df = pd.DataFrame({
            'Metric': metrics,
            'Value': values,
            'Unit': ['[log(1/y)]', '[y]', '[-]', '[log(y^2)]', '[y]', "", '[-]', '[-]', '[%]', '[y]', '[y]', '[-]', '[-]', '[%]', '[%]', '[-]', '[-]', "", '[-]', '[-]', '[y]', '[%]'],
            'Reference': ['Good 1952', 'Matheson & Winkler 1972', 'Good 1971, Friedman 1983', 'Dawid & Sebastiani 1999', 'Gneiting & Raftery 2007', "", 'Renard et al. 2011', 
                          'Evin et al. 2013', "", "", "", 'McInerney et al. 2017', 'McInerney et al. 2017', 'McInerney et al. 2017', 'McInerney et al. 2017', "", "", "", "", "", '[y] ♫', ""]
        })

        # Set pandas options to left-align the columns
        pd.set_option('display.colheader_justify', 'left')
        
        # Left-align the first column 'Metric'
        # df.style.set_properties(subset = ['Metric', 'Value', 'Unit', 'Reference'], **{'text-align': 'left'}).hide(axis = "index")
        # df.style.set_properties(**{'text-align': 'left'})
        # Left-align and pad the 'Metric' column for console output
        # df['Metric'] = df['Metric'].apply(lambda x: f"{x:<46}")
        # Round to 5 digits after the dot
        # df['Value'] = df['Value'].round(5)
        # Left-align and pad the 'Value' column for console output
        # df['Value'] = df['Value'].apply(lambda x: f"{x:<22}")
        # Left-align and pad the 'Unit' column for console output
        # df['Unit'] = df['Unit'].apply(lambda x: f"{x:<14}")
        # Left align the 'Unit' column
        # df.style.set_properties(subset = ['Unit'], **{'text-align': 'left'})

        # Left-align and pad the 'Reference' column for console output
        # df['Reference'] = df['Reference'].apply(lambda x: f"{x:<27}")
        # left align the column headers with a width of 30
#        df.columns = df.columns.str.ljust(30)
        # df.columns = df.columns.str.ljust(15)
        # print to screen without the first column [indexes: 0,1,2,3, etc.]
        print(df.to_string(index=False))
        
        # Export to PDF
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.axis('tight')
        ax.axis('off')
        ax.table(cellText=df.values, colLabels=df.columns, loc = 'center', cellLoc = 'center')
        plt.savefig(file_figs, format = 'pdf', bbox_inches = 'tight')
        
        # Saving the data to file (simulating .mat)
        np.save(file_data, {
            'chain': chain,
            'output': output,
            'MAP_un': MAP_un,
            'FX_MAP': FX_MAP,
            'DREAMPar': DREAMPar,
            'Par_info': Par_info,
            'Meas_info': Meas_info,
            'Lik_info': Lik_info,
            'options': options,
            'MAP_info': MAP_info,
            'par_unc': par_unc,
            'tot_unc': tot_unc
        })
        
    return FX, Z


def Calc_MTproposal(X, CR, DREAMPar, Table_gamma, Par_info, jump_method, step, gamma_jump=None, Z=None):
    ## ##################################################################### ##
    ##                                                                       ##
    ## Calculate candidate points using discrete proposal distribution       ##
    ##                                                                       ##
    ## SYNOPSIS:                                                             ##
    ##  [Xp,log_alfa_sn_X,v,Z,CR_mt,gamma] = Calc_MTproposal(X,CR, ...       ##
    ##   DREAMPar,Table_gamma,Par_info,jump_method,step,gamma,Z)             ##
    ## where                                                                 ##
    ##  X           [input] Nxd matrix of current chain states               ##
    ##  CR          [input] Nx1 vector of crossover values                   ##
    ##  DREAMPar    [input] Structure with algorithmic variables             ##
    ##  Table_gamma [input] Table with optimal DE-MC jumprates               ## 
    ##  Par_info    [input] Parameter structure: Ranges, initial/prior & bnd ##
    ##  jump_method [input] Scalar, 1: Parallel direct., 2: snooker update   ##
    ##  step        [input] Scalar, 1: Candidate points, 2: reference points ##
    ##  gamma       [input] Nx1 vector of jumprates if step 2                ##
    ##  Z           [input] archive of external samples if step 2            ##
    ##  X           [outpt] (mtxN)xd matrix of candidate points (step 1)     ##
    ##                      (mtx(N-1))xd matrix of reference points (step 2) ##
    ##  log_alfa_sn [outpt] (mtxN)x1 vector log snkr trans. dens. (step 1)   ##
    ##                      (mtx(N-1))x1 vector log snkr trans. dens. (st#2) ##
    ##  v           [outpt] (mtxN)x1 vector with 0/1 proposl in bound (st#1) ##
    ##                      (mtx(N-1))x1 vector 0/1 proposal in bound (st#2) ##
    ##  Z           [outpt] archive of external samples if step 1            ##
    ##  CR_mt       [outpt] mtxN matrix crossover multi-try props. (step 1)  ##
    ##  gamma       [outpt] Nx1 vector jumprates multi-try props. (step 1)   ##
    ##                                                                       ##
    ## Calculate candidate points using discrete proposal distribution       ##
    ## MATLAB function: Calc_MTproposal.mlx provides a detailed explanation  ##
    ## of the theory and implementation in DREAM_Suite                    ##
    ##                                                                       ##
    ##  © Written by Jasper A. Vrugt, Feb 2007                               ##
    ##  Los Alamos National Laboratory                                       ##
    ##                                                                       ##
    ## ##################################################################### ##

    if Z is None:
        Z = np.empty((DREAMPar['select'], DREAMPar['d'], DREAMPar['N']))
    if gamma_jump is None:
        gamma_jump = np.array([])

    Nmt = DREAMPar['mt'] * DREAMPar['N']
    if step == 1:  # Proposal Step
        MT = DREAMPar['mt']
        if jump_method == 1:
            gamma_de = np.random.choice([0, 1], size=DREAMPar['N'], p=[1 - DREAMPar['p_unit_gamma'], DREAMPar['p_unit_gamma']])
        elif jump_method == 2:
            gamma_sn = 1.2 + np.random.rand(DREAMPar['N'])

        row_size = (DREAMPar['d'] + 2) * np.dtype(np.float64).itemsize
        R = np.empty((DREAMPar['N'], DREAMPar['select']), dtype=int)
        # Open external archives with historical chain states
        with open('Z.bin', 'rb') as fid_Z:
            # added for python: but does this slow things down - reading the entire file?
            for i in range(DREAMPar['N']):
                # Without replacement draw rows from Z to create proposals: draw from 0 to m-1
                R[i, :DREAMPar['select']] = np.random.choice(DREAMPar['m'], DREAMPar['select'], replace=False)
                # Now read from the file
                for idx, row_idx in enumerate(R[i, :DREAMPar['select']]):
                    # Calculate the byte offset for the specific row
                    byte_offset = row_idx * row_size
                    # Move to the position of the selected row
                    fid_Z.seek(byte_offset, 0)  
                    # Read the row data from the file
                    Z[idx, :, i] = np.fromfile(fid_Z, dtype=np.float64, count=DREAMPar['d'])
                
    elif step == 2:  # Reference Step (Auxiliary Proposals)
        MT = DREAMPar['mt'] - 1
        if jump_method == 1:
            gamma_de = gamma_jump
        elif jump_method == 2:
            gamma_sn = gamma_jump

    # Copy crossover values MT times (columnwise - then transpose)
#    CR_mt = np.tile(CR, (MT, 1)).T
    CR_mt = np.tile(CR, (MT, 1))

    # Initialize log snooker space
    log_alfa_sn_X = np.full(Nmt, np.nan)

    # Initialize return matrix of candidate points
    Xp = np.full((Nmt, DREAMPar['d']), np.nan)

    # Multi-try sampling --> per chain
    for i in range(DREAMPar['N']):
        id_mt = np.arange(i * DREAMPar['mt'], (i + 1) * DREAMPar['mt'])
        # Set jump vectors of multiple trials equal to zero
        dX = np.zeros((MT, DREAMPar['d']))
        # Calculate the ergodicity perturbation
        eps = DREAMPar['zeta'] * np.random.randn(MT, DREAMPar['d'])
        # Randomness each individual chain
        rnd_jump = DREAMPar['lambda'] * (2 * np.random.rand(MT, DREAMPar['d']) - 1)
        # Unpack Z from 3d array to 2d array - easier to digest
        Z2 = Z[:, :, i]
        if jump_method == 1:  # Parallel Direction Proposal Distribution
            # number of chain pairs each individual MT proposal
            delta = np.random.randint(1, DREAMPar['delta'] + 1, MT)
            # Generate uniform random numbers for each MT proposal to determine which dimension to update
            rnd_cr = np.random.rand(MT, DREAMPar['d'])
            # Create MT proposals for chain i
            for mt in range(MT):
                # Define r1 and r2
                r1 = DREAMPar['R'][:delta[mt], mt]
                r2 = DREAMPar['R'][delta[mt]:2 * delta[mt], mt]
                if gamma_de[i] == 1:    # Update all dimensions and set jump rate to one
                    A = np.arange(DREAMPar['d'])
                    gamma_RWM = 1
                    CR_mt[mt, i] = DREAMPar['nCR']
                    A = np.array(A).flatten()  # Ensure 'A' is 1D
                    r1 = np.array(r1).flatten()  # Ensure 'r1' is 1D
                    r2 = np.array(r2).flatten()  # Ensure 'r2' is 1D
                else:                   # Derive subset A with dimensions to sample
                    A = np.where(rnd_cr[mt, :] < (CR_mt[mt, i] / DREAMPar['nCR']))[0]
                    if len(A) == 0:
                        A = np.random.permutation(DREAMPar['d'])[:1]
                    d_prime = len(A)
                    gamma_RWM = Table_gamma[d_prime - 1, delta[mt] - 1]         # gamma_RWM = Table_gamma[d_prime, delta[mt]]
                    A = np.array(A).flatten()  # Ensure 'A' is 1D                    
                    
                dX[mt, A] = (1 + rnd_jump[mt, A]) * gamma_RWM * np.sum(Z2[np.ix_(r1, A)] - Z2[np.ix_(r2, A)], axis = 0) + eps[mt, A]
                log_alfa_sn_X[id_mt[:MT]] = np.zeros(MT)

        elif jump_method == 2:  # Snooker Proposal Distribution
            # Loop over the individual chains
            for mt in range(MT):
                # Sample a, b, and c from external archive Z
                a = DREAMPar['R'][0, mt]
                b = DREAMPar['R'][1, mt]
                c = DREAMPar['R'][2, mt]
                # Define projection vector X(i,:) - Zc
                # F = X[i, :DREAMPar['d']] - Z[c, :DREAMPar['d'], i]
                F = X[i, :DREAMPar['d']] - Z2[c, :DREAMPar['d']]
                FF = max(np.dot(F, F), np.finfo(float).eps)
                # Orthogonally project zR1 and zR2 onto F
                zP = F * np.sum((Z2[a, :DREAMPar['d']] - Z2[b, :DREAMPar['d']]) * F) / FF
                dX[mt, :] = (1 + rnd_jump[mt, :]) * gamma_sn[i] * zP + eps[mt, :]
                CR_mt[mt, i] = DREAMPar['nCR']
                log_alfa_sn_X[id_mt[mt]] = (DREAMPar['d'] - 1) * np.linalg.norm(X[i, :DREAMPar['d']] - Z2[c, :DREAMPar['d']], 2) + np.finfo(float).eps

        Xp[id_mt[:MT], :] = np.tile(X[i, :DREAMPar['d']], (MT, 1)) + dX

    # Boundary handling ('bound'/'reflect'/'fold'/'reject')
    if 'boundhandling' in Par_info:
        Xp, v = Boundary_handling(Xp, Par_info)
        # added by JAV
        non_nan_rows = ~np.isnan(Xp).any(axis = 1)
        # remove boolean for nan rows [= only for reference set]
        v = v[non_nan_rows]            
    else:
        if step == 1:        
            v = np.zeros(Nmt, dtype=bool)                   # 0 for in-bound, 1 for out-of-bound
        elif step == 2:
            v = np.zeros(Nmt - DREAMPar['N'], dtype=bool)

    if jump_method == 1:
        return Xp, log_alfa_sn_X, v, Z, CR_mt, gamma_de
    elif jump_method == 2:
        return Xp, log_alfa_sn_X, v, Z, CR_mt, gamma_sn


def Calc_proposal(method, X, EX, std_mX, CR, DREAMPar, Table_gamma, Par_info, Meas_info):
    ## ####################################################################% ##
    ## This function generates candidate points using a discrete proposal    ##
    ##  distribution                                                         ##
    ##                                                                       ##
    ## SYNOPSIS: [Xp,log_alfa_sn,CR,v] = Calc_proposal(method,X,EX,...       ##
    ##              std_mX,CR,DREAMPar,Table_gamma,Par_info,Meas_info)       ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                ##
    ## Los Alamos National Laboratory                                        ##
    ##                                                                       ##
    ## ##################################################################### ##

    # Initialize jump vectors and candidate points
    Xp = np.zeros((DREAMPar['N'], DREAMPar['d']))
    dX = np.zeros((DREAMPar['N'], DREAMPar['d']))
    
    log_alfa_sn = np.zeros(DREAMPar['N'])                                                           # Define log of alfa snooker
    jump_method = np.random.choice([1, 2, 3], p = [DREAMPar['pparallel'], DREAMPar['psnooker'], \
        DREAMPar['pkalman']])                                                                       # Determine the method for the jump
    
    # Prepare arguments for creating candidate points
    if jump_method in [1, 2]:
        if method in ['dream', 'dream_d']:
            draw = np.argsort(np.random.rand(DREAMPar['N']-1, DREAMPar['N']), axis = 0)             # Randomly permute numbers [1,...,N-1] N times
            Z = X                                                                                   # External archive is simply equal to current chain states
        elif method in ['dream_zs', 'dream_dzs', 'dream_kzs']:
            draw = np.tile(np.arange(0, 6), (DREAMPar['N'], 1)).T                                   # Permute [1,...,N-1] N times so that Z samples differ
            R = np.random.choice(DREAMPar['m'], DREAMPar['select'], replace=False).astype(int)
            row_size = (DREAMPar['d'] + 2) * np.dtype(np.float64).itemsize
            Z = np.empty((DREAMPar['select'], DREAMPar['d']))
            with open('Z.bin', 'rb') as fid_Z:                                                      # Open external archives with historical chain states
                for idx, row_idx in enumerate(R):
                    byte_offset = row_idx * row_size                                                # Calculate the byte offset for the specific row
                    fid_Z.seek(byte_offset, 0)                                                      # Move to the position of the selected row
                    Z[idx, :] = np.fromfile(fid_Z, dtype=np.float64, count=DREAMPar['d'])           # Read the row data from the file

    elif jump_method == 3:                                                                          # dream_kzs: We use past states and their corresponding simulations
        R = np.random.choice(DREAMPar['m'], DREAMPar['M'], replace=False).astype(int)
        Z = np.empty((DREAMPar['M'], DREAMPar['d']))
        FXZ = np.empty((DREAMPar['M'], Meas_info['n']))
        row_size = (DREAMPar['d'] + 2) * np.dtype(np.float64).itemsize                              # Open external archive with historical chain states
        with open('Z.bin', 'rb') as fid_Z:
            for idx, row_idx in enumerate(R):
                byte_offset = row_idx * row_size                                                    # Calculate the byte offset for the specific row
                fid_Z.seek(byte_offset, 0)                                                          # Move to the position of the selected row
                Z[idx, :] = np.fromfile(fid_Z, dtype=np.float64, count=DREAMPar['d'])               # Read the row data from the file

        row_size = Meas_info['n'] * np.dtype(np.float64).itemsize                                   # Open external archives with simulated outputs of Z
        with open('FXZ.bin', 'rb') as fid_FXZ:
            for idx, row_idx in enumerate(R):
                byte_offset = row_idx * row_size                                                    # Calculate the byte offset for the specific row
                fid_FXZ.seek(byte_offset, 0)                                                        # Move to the position of the selected row
                FXZ[idx, :] = np.fromfile(fid_FXZ, dtype=np.float64, count=Meas_info['n'])          # Read the row data from the file

    # Determine how many chain pairs to use for each individual chain
    delta = np.random.randint(1, DREAMPar['delta'] + 1, size = DREAMPar['N'])                       # Draw delta from [1,...,Delta]
    rnd_cr = np.random.rand(DREAMPar['N'], DREAMPar['d'])                                           # Nxd matrix of U[0,1]
                                                                                                    # = np.random.uniform(0,1,size=(DREAMPar['N'], DREAMPar['d']))
    eps = DREAMPar['zeta'] * np.random.randn(DREAMPar['N'], DREAMPar['d'])                          # Nxd matrix of ζ x N(0,1)
                                                                                                    # = ζ x np.random.normal(0,1,size=(DREAMPar['N'], DREAMPar['d']))
    rnd_jump = DREAMPar['lambda'] * np.random.uniform(-1,1,size=(DREAMPar['N'], DREAMPar['d']))     # Nxd matrix of λ x U[-1,1]    
                                                                                                    # = λ x (2 x np.random.rand(DREAMPar['N'], DREAMPar['d']) - 1)
    rnd_lambda = 1 + rnd_jump                                                                       # Nxd matrix of λ x U[0,2]                                                                                              

    if jump_method == 1:        ## PARALLEL DIRECTION PROPOSAL DISTRIBUTION
        gamma_jump = np.random.choice([0, 1], size = DREAMPar['N'], p = [1 - DREAMPar['p_unit_gamma'], DREAMPar['p_unit_gamma']])
        for i in range(DREAMPar['N']):
            a = DREAMPar['R'][i, draw[:delta[i], i]]
            b = DREAMPar['R'][i, draw[delta[i]:2*delta[i], i]]
            if gamma_jump[i] == 0:  # Subspace sampling            
                A = np.where(rnd_cr[i, :DREAMPar['d']] < (CR[i] / DREAMPar['nCR']))[0]
                d_prime = len(A)
                if d_prime == 0:
                    A = np.argmin(rnd_cr[i, :])                                                     # id = np.random.permutation(DREAMPar['d']), A = [id[0]]
                    d_prime = 1
                gamma_d = Table_gamma[d_prime - 1, delta[i] - 1]                                    # Get jump rate from tabulated values: gamma_d = Table_gamma[d_prime, delta[i] ]
                A = np.array(A).flatten()                                                           # Ensure 'A' is 1D                    
            else:               # All dimensions update
                A = np.arange(DREAMPar['d'])
                gamma_d = 1
                CR[i] = DREAMPar['nCR']
                a = a[0]; b = b[0]                                                                  # JAV: this line was added in 2025, single pair if gamma_jump is 1, OK, no big deal
                a = np.array(a).flatten()                                                           # Ensure 'a' is 1D
                b = np.array(b).flatten()                                                           # Ensure 'b' is 1D
                A = np.array(A).flatten()                                                           # Ensure 'A' is 1D                    
            
            dX[i, A] = rnd_lambda[i, A] * gamma_d * np.sum(Z[np.ix_(a, A)] \
                                                           - Z[np.ix_(b, A)], axis = 0) + eps[i, A] # dX[i, A] = (1 + rnd_jump[i, A]) * gamma_d * np.sum(Z[np.ix_(a, A)] - Z[np.ix_(b, A)], axis = 0) + eps[i, A]
        Xp[:, :DREAMPar['d']] = X[:, :DREAMPar['d']] + dX                                           # Compute candidate point            

    elif jump_method == 2:      ## SNOOKER PROPOSAL DISTRIBUTION
        a = DREAMPar['R'][:, 0]
        b = DREAMPar['R'][:, 1]
        c_sn = DREAMPar['R'][:, 2]
        for i in range(DREAMPar['N']):
            F = X[i, :DREAMPar['d']] - Z[c_sn[i], :DREAMPar['d']]                                   # a = DREAMPar['R'][i, 0], b = DREAMPar['R'][i, 1], c = DREAMPar['R'][i, 2]
            FF = np.maximum(F @ F.T, np.finfo(float).eps)                                           # (1xd) * (dx1) = 1x1
            zP = F * (np.sum((Z[a[i], :DREAMPar['d']] - Z[b[i], :DREAMPar['d']]) * F) / FF)         # zp = np.dot(F, np.sum((Z[a, :DREAMPar['d']] - Z[b, :DREAMPar['d']]) * F, axis = 1) / FF)
            gamma_s = 1.2 + np.random.rand()
            dX[i, :DREAMPar['d']] = rnd_lambda[i, :DREAMPar['d']] * gamma_s * zP \
                + eps[i, :DREAMPar['d']]                                                            # dX[i, :DREAMPar['d']] = (1 + rnd_jump[i, :DREAMPar['d']]) * gamma_s * zP + eps[i, :DREAMPar['d']]
            # ################## Serial implementation of snooker correction ###################
            # Xp[i, :] = X[i, :DREAMPar['d']] + dX[i, :]
            # XpZ = np.linalg.norm(Xp[i, :DREAMPar['d']] - Z[c_sn[i], :DREAMPar['d']], 2) + np.finfo(float).eps
            # XZ = np.linalg.norm(X[i, :DREAMPar['d']] - Z[c_sn[i], :DREAMPar['d']], 2) + np.finfo(float).eps
            # log_alfa_sn[i] = (DREAMPar['d'] - 1) * np.log(XpZ / XZ)
            # ##################################################################################
        Xp[:, :DREAMPar['d']] = X[:, :DREAMPar['d']] + dX                                           # Compute candidate points
        log_alfa_sn = (DREAMPar['d'] - 1)/2 * ( np.log(np.sum((Xp[:DREAMPar['N'], :DREAMPar['d']] \
            - Z[c_sn, :DREAMPar['d']]) ** 2, axis = 1)) - np.log(np.sum((X[:DREAMPar['N'], \
                :DREAMPar['d']] - Z[c_sn, :DREAMPar['d']]) ** 2, axis = 1) + np.finfo(float).eps) ) # Compute snooker correction: vectorized

    # Python implementation of Kalman jump [= did not verify all different implementations with "R", "std_mX"]
    elif jump_method == 3:  # KALMAN PROPOSAL DISTRIBUTION
        m_Z = np.tile(np.mean(Z, axis = 0), (DREAMPar['M'], 1))                                     # Compute mean of each parameter
        m_FX = np.tile(np.mean(FXZ, axis = 0), (DREAMPar['M'], 1))                                  # Compute mean of output at each time/coordinate
        C_ZFXZ = (Z - m_Z).T @ (FXZ - m_FX) / (DREAMPar['M'] - 1)                                   # Compute cross-covariance between Z and FXZ, (Mxd)' * (Mxn) = dxn = C_ZFXZ
        C_FXZFXZ = ((FXZ - m_FX).T @ (FXZ - m_FX)) / (DREAMPar['M'] - 1)                            # Compute autocovariance of model output of Z, (Mxn)' * (Mxn) = nxn = C_FXZFXZ

        # In following: KXZ = dxn matrix = type of Kalman gain
        # Note: The loop below can be written differently to avoid repeated calculation of constant KXZ in some cases
        # Theory: C_ZFXZ/C_FXZFXZ in MATLAB is equivalent to (inv(C_FXZFXZ)*C_ZFXZ)' 
        for i in range(DREAMPar['N']):
            if 'R' not in Meas_info:
                if np.all(np.isnan(std_mX)):
                    R = 0                                                                           # nxn measurement error matrix
                    r = np.zeros(Meas_info['n'])                                                    # Residual vector is zeroth vector
                else:
                    R = std_mX[:Meas_info['n'], i] ** 2 * np.eye(Meas_info['n'])                    # nxn measurement error matrix
                    # KXZ = C_ZFXZ @ np.linalg.inv(C_FXZFXZ + np.cov(std_mX[:Meas_info['n'], i], \
                    #     rowvar = False) * np.eye(Meas_info['n']))                                 # dxn * (nxn) = dxn 
                    r = np.random.randn(Meas_info['n']) * std_mX[:Meas_info['n'], i]                # Residual vector is N_{n}(0,Σ = diag(std_mX^2))
            elif Meas_info['R']: 
                R = Meas_info['R']                                                                  # nxn measurement error matrix            
                if np.all(np.isnan(std_mX)):
                    r = np.random.multivariate_normal(np.zeros(Meas_info['n']), Meas_info['R'])     # Residual vector is N_{n}(0,Σ = diag(R))
                else:
                    r = np.random.randn(Meas_info['n']) * std_mX[:Meas_info['n'], i]                # Residual vector is N_{n}(0,Σ = diag(std_mX^2))
            elif Meas_info['R'] == []:
                if np.all(np.isnan(std_mX)):
                    R = 0                                                                           # nxn measurement error matrix
                    r = np.zeros((Meas_info['n'], 1))                                               # Residual vector is zeroth vector
                else:
                    R = std_mX[:Meas_info['n'], i] ** 2 * np.eye(Meas_info['n'])                    # nxn measurement error matrix
                    # KXZ = C_ZFXZ @ np.linalg.inv(C_FXZFXZ + np.cov(std_mX[:Meas_info['n'], i], \
                    #     rowvar = False) * np.eye(Meas_info['n']))                                 # dxn * (nxn) = dxn
                    r = np.random.randn(Meas_info.n, 1) * std_mX[:Meas_info['n'], i]                # Residual vector is N_{n}(0,Σ = diag(std_mX^2))

            # KXZ = C_ZFXZ @ np.linalg.inv(C_FXZFXZ + R)                                              # Kalman gain: (dxn) * (nxn) = dxn
            L = la.cholesky(C_FXZFXZ + R, lower = True)                                             # nxn matrix
            Y = la.solve(L, C_ZFXZ.T)                                                               # nxd matrix
            KXZ = la.solve(L.T, Y).T                                                                # dxn matrix = Kalman gain
            dx_KXZ = KXZ @ (r - EX[:Meas_info['n'], i])                                             # (dxn) * (nx1) = dx1 Kalman jump vector of parameters
            dX[i, :] = (1 + rnd_jump[i, :DREAMPar['d']]) * dx_KXZ.T + eps[i, :DREAMPar['d']]        # Compute actual 1xd shift vector of parameters: (1xd) = (1xd) .* (1xd) + (1xd)
            Xp[i, :] = X[i, :DREAMPar['d']] + dX[i, :]                                              # Compute candidate point
            
    # Boundary handling
    if 'boundhandling' in Par_info:
        Xp, v = Boundary_handling(Xp, Par_info)
    else:
        v = np.zeros(DREAMPar['N']) > 1
    
    if method in ['dream_d', 'dream_dzs']:
        Xp = Discrete_space(Xp, Par_info)
    
    return Xp, log_alfa_sn, CR, v


def Eval_prior(X, SX, Par_info, Meas_info, options):
    ## ##################################################################### ##
    ## This function evaluates the prior distribution and returns the        ##
    ##  log-prior of the N parameter vectors of X                            ##
    ##                                                                       ##
    ## SYNOPSIS: logPR_X = Eval_prior(X,SX,DREAMPar,Par_info,Meas_info,...   ##
    ##                         options)                                      ##
    ##  where                                                                ##
    ##   X         [input] N x d matrix of parameter vectors                 ##
    ##   SX        [input] m x N matrix of m summary metrics paramtr vectors ##
    ##   Par_info  [input] Parameter structure                               ##
    ##   Meas_info [input] Measurement data structure                        ##
    ##   options   [input] Structure with computational settings/options     ##
    ##   logPR_X   [outpt] Nx1 vector with logarithmic prior density of X    ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb. 2007                               ##
    ## Los Alamos National Laboratory                                        ##
    ##                                                                       ##
    ## ##################################################################### ##

    N, d = X.shape          # Number of parameter vectors (N) and their dimensions (d)
    logPR_X = np.zeros(N)   # Initialize log-prior

    # Check if prior information exists
    if 'prior' in Par_info:
        if Par_info['u'] == 'yes':                          # Univariate case
            PR_X = np.full((N, d), np.nan)                  # Initialize prior densities
            for ii in range(len(Par_info['prior'])):
                #par_name = Par_info['prior'][ii]
                #fnc = eval(f'{par_name}')                
                fnc = Par_info['prior'][ii]
                PR_X[:, ii] = np.maximum(fnc.pdf(X[:, ii]),np.finfo(np.float64).tiny)   

            if Par_info['pr'] == 'pdf':
                logPR_X = np.sum(np.log(PR_X), axis = 1)    # Prior = pdf --> take log density
            elif Par_info['pr'] == 'logpdf':
                logPR_X = np.sum(PR_X, axis = 1)            # Prior = logpdf --> no log needed

        elif Par_info['u'] == 'no':                         # Multivariate case
            PR_X = np.full((N, 1), np.nan)                  # Initialize prior densities
            fnc = Par_info['prior'][0]
            for jj in range(N):                             # Evaluate joint prior density at once
                PR_X[jj, 0] = fnc.pdf(X[jj, :])
                
            if Par_info['pr'] == 'pdf':
                logPR_X = np.log(PR_X)                      # Prior = pdf --> take log density
            elif Par_info['pr'] == 'logpdf':
                logPR_X = PR_X                              # Prior = logpdf --> no log needed

    # Check if we need to do diagnostic Bayes
    if options.get('DB') == 'yes':
        for ii in range(N):
            logPR_X[ii] = np.min(options['epsilon'] - np.abs(Meas_info['S'] - SX[:, ii]))

    return logPR_X


def Calc_likelihood(Xp, FXp, DREAMPar, Par_info, Meas_info, Lik_info, options, nargout, MAP_info = None):
    ## ################################################################################## ##
    ## This function computes the log-likelihood of the candidate points Xp               ##
    ##                                                                                    ##
    ## SYNOPSIS: logPR_Xp = Calc_likelihood(Xp,FXp,DREAMPar,Par_info,Meas_info,...        ##
    ##               Lik_info,options,MAP_info)                                           ##
    ##  where                                                                             ##
    ##   Xp        [input] N x d matrix of candidate points                               ##
    ##   FXp       [input] n x N matrix with simulated values candidate points            ##
    ##   DREAMPar  [input] Algorithmic structure with variables MCMC method               ##
    ##   Par_info  [input] Parameter structure                                            ##
    ##   Meas_info [input] Measurement data structure                                     ##
    ##   Lik_info  [input] Structure with information likelihood function                 ##
    ##   options   [input] Structure with computational settings/options                  ##
    ##   MAP_info  [input] OPT: Information MAP solution for sandwich estimator           ##
    ##   loglik_Xp [outpt] Nx1 vector with logarithmic value of likelihood of Xp          ##
    ##   E         [outpt] nxN matrix with residuals of Xp                                ##
    ##   std_m     [outpt] nxN matrix with measurement errors of Xp for Kalman proposal   ##
    ##   eps_n     [outpt] nxN matrix of standardized partial residuals of Xp             ##
    ##   f_eps_n   [outpt] nxN matrix of density standardized partial residuals of Xp     ##
    ##   ell       [outpt] nxN matrix of log-likelihoods of Xp                            ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                             ##
    ## Los Alamos National Laboratory                                                     ##
    ##                                                                                    ##
    ## ################################################################################## ##

    if MAP_info is None:
        MAP_info = {}

    n = Meas_info['n']                  # Number of measurements
    n_S = Meas_info['n_S']              # Number of summary metrics
    N = Xp.shape[0]                     # Number of candidate vectors (can differ from DREAMPar.N)
    loglik_Xp = np.nan * np.ones(N)     # Initialize log-likelihood for candidate points
    E = []                              # Initialize residuals
    if n > 0:                           # FXp is a simulation
        E = Meas_info['Y'].reshape(-1, 1) - FXp
        
    ell = np.full((max(max(n_S, n),1),N),np.nan)
    # For Python
    FXp = np.array(FXp)
    if 'rho' in options or DREAMPar['lik'] in [21, 22]:
        rho_name = options['rho']
        fnc_rho = eval(f'{rho_name}')

    # User-specified measurement error?
    if Meas_info['Sigma'] != None:
        std_e = Meas_info['Sigma']
        std_m = np.tile(std_e, (1, N))
    else:
        std_m = np.nan * np.ones((Meas_info['n'], N))
        
    for ii in range(N):  # Loop over each candidate point

        if DREAMPar['lik'] == 1:
            ell[:, ii] = np.log(FXp[:, ii])
            loglik_Xp[ii] = np.sum(ell[:, ii])
        
        elif DREAMPar['lik'] == 2:
            ell[:, ii] = FXp[:, ii]
            loglik_Xp[ii] = np.sum(ell[:, ii])
            
        elif DREAMPar['lik'] == 11:  # Measurement error integrated out
            loglik_Xp[ii] = -(n / 2) * np.log(np.sum(np.abs(E[:n, ii]) ** 2))
            ell[0, ii] = loglik_Xp[ii]

        elif DREAMPar['lik'] == 12:  # Normal likelihood with homos/heteroscedastic measurement error
            e_n = E[:n, ii] / std_e  # Standardize residuals
            ell[:n, ii] = -0.5 * np.log(2 * np.pi) - np.log(std_e) - 0.5 * e_n ** 2
            loglik_Xp[ii] = np.nansum(ell[:n, ii])

        elif DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:   # NL/GL/LAPL/SSTL/GL+/SGTL functions
            par = Lik_info['fpar'].copy()
            par[Lik_info['id_vpar']] = Xp[ii, :DREAMPar['d']]
            #for i, idx in enumerate(Lik_info['id_vpar']):
            #    par[idx] = Xp[ii,i]
            nuisvar = par[Lik_info['id_nuisvar']].copy()    # Isolate nuisance variables
            exec_scope = {                                  # Must add this executive dictionary
                'nuisvar': nuisvar,                         # Input variable
                'Lik_info': Lik_info,                       # The Lik_info dictionary with the function string
                'FXp': FXp,                                 # Output of candidate point
                'ii': ii,                                   # Index of candidate point
                'n': n,                                     # Number of data points of observation [data] vector
                'Meas_info': Meas_info,                     # Dictionary of measurement data/information
                'loglik': None,                             # Return argument 1: nx1 vector of log-likelihoods
                'std_e': None,                              # Return argument 2: nx1 vector of measurement error stds.
                'eps_n': None,                              # Return argument 3: nx1 vector of studentized partial residuals
                'f_eps_n': None,                            # Return argument 4: nx1 vector of density of studentized partial residuals
                '__': None}                                 # Return argument 5: matrix of replicates 

            exec_scope[Lik_info['name_lik_func']] = globals()[Lik_info['name_lik_func']]    # Add name of likelihood function
            exec(Lik_info['stringL'], exec_scope)                                           # Evaluate likelihood function
            loglik = exec_scope['loglik'].squeeze()
            std_e = exec_scope['std_e'].squeeze()
            eps_n = exec_scope['eps_n']
            f_eps_n = exec_scope['f_eps_n']
            ell[:, ii] = loglik 
            loglik_Xp[ii] = np.sum(loglik)
            std_m[:, ii] = std_e                            # For Kalman jump: std_e will equal Meas_info['Sigma'] if this has been defined by user            

        elif DREAMPar['lik'] == 15:  # Whittle quasi-likelihood function
            loglik_Xp[ii] = Whittle_loglik(FXp[:n, ii], Meas_info)
            ell[0, ii] = loglik_Xp[ii]
            
        elif DREAMPar['lik'] == 21:  # Approximate Bayesian Computation
            # phi1 = options['rho'](FXp[:n_S, ii], Meas_info['S']) + np.random.normal(0, options['epsilon'])
            phi1 = fnc_rho(FXp[:n_S, ii], Meas_info['S']) + np.random.normal(0, options['epsilon'])
            ell[:, ii] = -0.5 * np.log(2 * np.pi) - np.log(options['epsilon']) - 0.5 * (phi1 / options['epsilon']) ** 2
            loglik_Xp[ii] = np.sum(ell[:, ii])

        elif DREAMPar['lik'] == 22:  # Approximate Bayesian Computation
            # loglik_Xp[ii] = np.min(options['epsilon'] - options['rho'](FXp[:n_S, ii], Meas_info['S']))
            loglik_Xp[ii] = np.min(options['epsilon'] - fnc_rho(FXp[:n_S, ii], Meas_info['S']))
            ell[0, ii] = loglik_Xp[ii]

        elif DREAMPar['lik'] == 23:  # Approximate Bayesian Computation (limits of acceptability)
            ell[:, ii] = np.abs(FXp[:n_S, ii] - Meas_info['S']) <= options['epsilon']
            loglik_Xp[ii] = np.sum(ell[:, ii])

        elif DREAMPar['lik'] == 31:  # GLUE: Option a in Beven and Freer, 2001
            loglik_Xp[ii] = - DREAMPar['GLUE'] * np.log(np.var(E[:, ii]))

        elif DREAMPar['lik'] == 32:  # GLUE: Option b in Beven and Freer, 2001
            a = np.var(E[:, ii])
            b = np.var(Meas_info['Y'])
            if a < b:
                loglik_Xp[ii] = DREAMPar['GLUE'] * np.log(1 - a / b)
            else:
                loglik_Xp[ii] = - (a - b) * 1000 * DREAMPar['GLUE']

        elif DREAMPar['lik'] == 33:  # GLUE: Option c in Beven and Freer, 2001
            loglik_Xp[ii] = -DREAMPar['GLUE'] * np.var(E[:, ii])

        elif DREAMPar['lik'] == 34:  # GLUE: Beven and Binley, 1992
            loglik_Xp[ii] = -np.log(np.sum(np.abs(E[:, ii])))

        elif DREAMPar['lik'] == 52:  # Matrix implementation of Gaussian likelihood
            if 'C' in Meas_info:
                loglik_Xp[ii] = -(n / 2) * np.log(2 * np.pi) - 0.5 * np.log(Meas_info['detC']) - 0.5 * E[:, ii].T @ Meas_info['invC'] @ E[:, ii]
            else:
                loglik_Xp[ii] = -(n / 2) * np.log(2 * np.pi) - np.sum(np.log(std_e)) - 0.5 * E[:, ii].T @ np.diag(1 / std_e ** 2) @ E[:, ii]

        elif DREAMPar['lik'] == 61:  # Laplace likelihood with variable learning rate
            lambda_ = Xp[ii, -2]
            b = std_e / np.sqrt(2)
            f_lambda = 1 / (2 * b) ** lambda_ * np.exp(-np.abs(E[:, ii] / b)) ** lambda_
            Z_lambda = 1 / (lambda_ * (2 * b) ** (lambda_ - 1))
            loglik_Xp[ii] = np.sum(np.log(f_lambda)) - np.sum(np.log(Z_lambda))

        elif DREAMPar['lik'] == 62:  # Normal likelihood with variable learning rate
            lambda_ = Xp[ii, -2]
            f_lambda = 1 / (std_e * np.sqrt(2 * np.pi)) ** lambda_ * np.exp(-0.5 * (E[:, ii] / std_e) ** 2) ** lambda_
            Z_lambda = 1 / (np.sqrt(lambda_) * (std_e * np.sqrt(2 * np.pi)) ** (lambda_ - 1))
            loglik_Xp[ii] = np.sum(np.log(f_lambda)) - np.sum(np.log(Z_lambda))

    # Sandwich adjustment of log-likelihood
    if 'map' in MAP_info:
        a = np.nan * np.ones(N)
        loglik_adj = np.nan * np.ones(N)
        
        if Par_info['norm'] == 1:
            Xp_n = X_normalize(Xp, Par_info)
        else:
            Xp_n = Xp
        
        for ii in range(N):
            dx = Xp_n[ii, :DREAMPar['d']] - MAP_info['map']
            a[ii] = (dx @ MAP_info['Godambe'] @ dx) / (dx @ MAP_info['Fisher'] @ dx)
            loglik_adj[ii] = a[ii] * (loglik_Xp[ii] - MAP_info['loglik'])

        if np.any(a < 0):
            print("DREAM-Suite WARNING: Negative value of a(ii) in Sandwich correction")

        loglik_Xp = loglik_adj

    # Return the appropriate output based on the number of requested arguments
    if nargout == 1:
        return loglik_Xp
    elif nargout == 3:
        return loglik_Xp, E, std_m
    elif nargout == 5:
        if DREAMPar['lik'] in [1, 2, 21, 22, 23]:
            e_n = []
            f_e_n = []
        elif DREAMPar['lik'] in [11, 15, 31, 32, 33, 34]:
            std_e = np.std(E[:n, 0])
            e_n = E[:n, 0] / std_e
            f_e_n = norm.pdf(e_n, 0, 1)
        elif DREAMPar['lik'] == 12:
            f_e_n = norm.pdf(e_n, 0, 1)
        elif DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:
            pass  # Do nothing; eps_n and f_eps_n should be known
        elif DREAMPar['lik'] == 52:
            if 'C' in Meas_info:
                W = np.linalg.cholesky(Meas_info['invC'])
                e_n = W @ E[:n, 0]
            else:
                e_n = E[:n, 0] / std_e
            f_e_n = norm.pdf(e_n, 0, 1)
        else:
            e_n = E[:n, 0] / std_e
            f_e_n = norm.pdf(e_n, 0, 1)

        if DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:
            return loglik_Xp, E, std_m, eps_n, f_eps_n
        else:
            return loglik_Xp, E, std_m, e_n, f_e_n
    elif nargout == 6:
        eps_n = []
        f_eps_n = []
        return loglik_Xp, E, std_m, eps_n, f_eps_n, ell

def Boundary_handling(X, Par_info):
    ## ##################################################################### ##
    ## This function checks parameters are in prior bounds, corrects them    ##
    ##                                                                       ##
    ## SYNOPSIS: [Xr, v] = Boundary_handling(X, Par_info)                    ##
    ##                                                                       ##
    ##   X          [input] Nxd matrix of candidate points                   ##
    ##   Par_info   [input] Dictionary of parameter information              ##
    ##   Xr         [outpt] Nxd matrix of revised candidate points           ##
    ##   v          [outpt] Nx1 vector (0: fine, 1: revised)                 ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                ##
    ## Los Alamos National Laboratory 			                             ##
    ##                                                                       ##
    ## ##################################################################### ##
    
    Xr = np.copy(X)                                         # Create a copy of X for the revised values
    N, d = X.shape
    v = np.zeros(N, dtype=bool)                             # Logical array indicating out-of-bound values
    
    mn = np.tile(Par_info['min'], (N, 1))                   # Lower bounds replicated N times
    mx = np.tile(Par_info['max'], (N, 1))                   # Upper bounds replicated N times
    
    # Positions where X is below lower bound or above upper bound
    id_l = np.where(X < mn)                                 # Smaller than lower bound
    id_u = np.where(X > mx)                                 # Larger than upper bound
    
    # Boundary handling options
    if Par_info['boundhandling'] == 'reflect':              # Reflection method
        Xr[id_l] = 2 * mn[id_l] - X[id_l]                   # Reflect below the lower bound
        Xr[id_u] = 2 * mx[id_u] - X[id_u]                   # Reflect above the upper bound
    elif Par_info['boundhandling'] == 'bound':              # Bound method
        Xr[id_l] = mn[id_l]                                 # Set to lower bound
        Xr[id_u] = mx[id_u]                                 # Set to upper bound
    elif Par_info['boundhandling'] == 'fold':               # Folding method
        Xr[id_l] = mx[id_l] - (mn[id_l] - X[id_l])          # Fold below the lower bound
        Xr[id_u] = mn[id_u] + (X[id_u] - mx[id_u])          # Fold above the upper bound
    elif Par_info['boundhandling'] == 'reject':             # Reject method
        o = np.zeros_like(X, dtype=bool)                    # Initialize out-of-bound array
        o[id_l] = 1                                         # Mark positions below the lower bound
        o[id_u] = 1                                         # Mark positions above the upper bound
        v = np.sum(o, axis = 1) > 0                         # Identify rows with any out-of-bound values
    
    # Reflection or folding: Check if all elements are within bounds
    # Both methods can go out of bounds if violation exceeds |mx - mn|
    if Par_info['boundhandling'] in ['reflect', 'fold']:
        id_l = np.where(Xr < mn)                            # Smaller than lower bound
        id_u = np.where(Xr > mx)                            # Larger than upper bound
        Xr[id_l] = np.random.uniform(mn[id_l], mx[id_l])    # Random draw in [mn, mx]
        Xr[id_u] = np.random.uniform(mn[id_u], mx[id_u])    # Random draw in [mn, mx]

    # BMA model training if applicable
    # if 'unit_simplex' in Par_info:
    #     wght_sum = np.sum(Xr[:N, :int(K)], axis = 1)
    #     Xr[:, :int(K)] = Xr[:, :int(K)] / wght_sum[:, np.newaxis]  # Normalize weights in the unit simplex

    return Xr, v


def select_Xp(DREAMPar, log_PRL_Xp, log_sn_Xp):
    ## ##################################################################### ##
    ## Computes importance weights of multi-try candidate points each chain  ##
    ##                                                                       ##
    ## SYNOPSIS: [log_w,id_pop] = select_Xp(DREAMPar, log_PRL_Xp, log_sn_Xp) ##
    ##                                                                       ##
    ##   DREAMPar   [input] Dictionary containing algorithmic parameters     ##
    ##   log_PRL_Xp [input] mtxN matrix log-prior + log-likelihoods          ##
    ##   log_sn_Xp  [input] mtxN matrix log snooker corrections, log(g(x))   ##
    ##   log_w      [outpt] mtx1 vector of log_prior + loglik - log g(x)     ##
    ##                      based on Eq. 10 in the referenced paper          ##
    ##   id_pop     [outpt] mtx1 vector of indices best proposals mt samples ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                ##
    ## Los Alamos National Laboratory 			                             ##
    ##                                                                       ##
    ## ##################################################################### ##

    # Compute log_w based on Eq. 10 in the referenced paper
    log_w = np.sum(log_PRL_Xp, axis = 1) - log_sn_Xp            # Sum log_prior + loglik - log g(x)
    log_w = np.maximum(log_w, -np.finfo(np.float64).max)        # Avoid -Inf or NaN
    log_w = np.minimum(log_w, np.finfo(np.float64).max)         # Avoid Inf

    # Reshape log_w to separate columns for each chain
    # log_w = log_w.reshape(DREAMPar['mt'], DREAMPar['N'])
    log_w = log_w.reshape(DREAMPar['mt'], DREAMPar['N'], order = 'F')
    # Rescale log_w to avoid numerical underflow (subtract max of log densities each column)
    log_ws = log_w - np.max(log_w, axis = 0, keepdims = True)
    # Compute rescaled weights
    ws = np.exp(log_ws) 
    # need to do check here about NaN values - otherwise code may crash

    # Normalize weights for each chain
    w = ws / np.sum(ws, axis = 0, keepdims=True)
    # Initialize column vector id_pop
    id_pop = np.full(DREAMPar['N'], np.nan)  
    # Now select the best among the DREAMPar.mt proposals in each chain
    for i in range(DREAMPar['N']):
        # Select the "best" proposal from the normalized weights using resampling
        id = np.random.choice(DREAMPar['mt'], size = 1, p = w[:DREAMPar['mt'], i])[0]
        # Determine the corresponding index in the flattened list of candidate points
        id_pop[i] = id + (i) * DREAMPar['mt']

    # Python addition
    id_pop = id_pop.astype(int)

    return log_w, id_pop


def mtMetropolis_rule(log_wp, log_wr, DREAMPar, id_pop):
    ## ##################################################################### ##
    ## Multi-try Metropolis rule for acceptance or rejection of proposals    ##
    ##                                                                       ##
    ## SYNOPSIS: [accept,id,ii] = mtMetropolis_rule(log_wp, log_wr, ...      ##
    ##                                              DREAMPar, id_pop)        ##
    ##                                                                       ##
    ##   log_wp     [input] mtxN matrix of log densities                     ##
    ##   log_wr     [input] mtxN matrix of log densities                     ##
    ##   DREAMPar   [input] Dictionary containing algorithmic parameters     ##
    ##   id_pop     [input] Nx1 vector of selected proposals of MT samples   ##
    ##   accept     [outpt] Nx1 vector of 1 (accept) or 0 (reject) proposals ##
    ##   id         [outpt] Vector of indices of accepted proposals          ##
    ##   ii         [outpt] Vector of boolean indices for accepted proposals ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, Feb 2007                                ##
    ## Los Alamos National Laboratory 			                             ##
    ##                                                                       ##
    ## ##################################################################### ##

    accept = np.zeros(DREAMPar['N'], dtype=int)                     # Initialize accept
    log_wr = np.maximum(log_wr, -np.finfo(np.float64).max)          # Handle -Inf or NaN
    log_wr = np.minimum(log_wr, np.finfo(np.float64).max)           # Handle Inf
    # Maximum of log(w) in each chain
    mx = np.maximum(np.max(log_wp), np.max(log_wr))
    wp = np.maximum(np.exp(log_wp - mx), np.finfo(np.float64).tiny) # Candidate weights
    wr = np.maximum(np.exp(log_wr - mx), np.finfo(np.float64).tiny) # Reference weights
    # Metropolis acceptance probability
    alfa = np.sum(wp, axis = 0) / np.sum(wr, axis = 0)
    # Draw uniform random numbers
    Z = np.random.rand(DREAMPar['N'])                               
    # Accept proposals where alfa > Z
    ii = alfa > Z
    accept[ii] = 1          # Accepted chains receive a 1
    id = id_pop[ii]         # Corresponding samples id_pop

    return accept, id, ii


def Metropolis_rule(DREAMPar, Meas_info, log_alfa_sn, logLPR_Xp, logLPR_X, options):
    ## ##################################################################### ##
    ## Metropolis rule for acceptance or rejection of candidate points       ##
    ##                                                                       ##
    ## SYNOPSIS: [accept,id_acc] = Metropolis_rule(DREAMPar, Meas_info, ...  ##
    ##                                            log_alfa_sn,logLPR_Xp, ... ##
    ##                                            logLPR_X, options)         ##
    ##                                                                       ##
    ##   DREAMPar   [input] Dictionary containing algorithmic parameters     ##
    ##   Meas_info  [input] Dictionary of measurement data [fitting]         ##
    ##   log_alfa_sn[input] Nmtx1 vector of logs snooker correction          ##
    ##   logLPR_Xp  [input] Nmtx2 matrix log-prior+loglik candidate points   ##
    ##   logLPR_X   [input] Nmtx2 matrix log-prior+loglik current states     ##
    ##   options    [input] Dictionary of algorithmic settings               ##
    ##   accept     [outpt] Nmtx1 vector of 1 (accept) or 0 (reject) points  ##
    ##   id_acc     [outpt] Vector of indices of accepted proposals          ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, March 2005                              ##
    ## Los Alamos National Laboratory 			                             ##
    ##                                                                       ##
    ## ##################################################################### ##

    logPR_Xp = logLPR_Xp[:, 0]  # Log-prior density of candidate points
    logL_Xp = logLPR_Xp[:, 1]   # Log-likelihood of candidate points
    logPR_X = logLPR_X[:, 0]    # Log-prior density of current chain states
    logL_X = logLPR_X[:, 1]     # Log-likelihood of current chain states

    # a_L = np.exp(logL_Xp - logL_X)  # Likelihood ratio
    # --> must rewrite acceptance rule: taking sum of 
    min_exp_value = -700; max_exp_value = 700               # Reasonable lower and upper bounds for float64
    Z = np.random.rand(DREAMPar['N'])                       # Draw standard uniform labels
    accept = np.full(DREAMPar['N'], np.nan)                 # Initialize acceptance vector

    # Check for ABC (Approximate Bayesian Computation) and Metropolis rule
    if options['ABC'] == 'no':                              # Regular MCMC with prior and likelihood
        if options['DB'] == 'no':                           # No Diagnostic Bayes
            #a_PR = np.exp(logPR_Xp - logPR_X)              # Prior ratio
            #a_S = np.exp(log_alfa_sn)                      # Snooker correction
            #alfa = a_S * a_L * a_PR                        # Product of ratios
            alfa = np.exp(np.clip((logPR_Xp - logPR_X) + (logL_Xp - logL_X) + log_alfa_sn, min_exp_value, max_exp_value))
            accept = alfa >= Z                              # Accept if alfa ≥ Z

        elif options['DB'] == 'yes':                        # Diagnostic Bayes
            a_L = np.exp(logL_Xp - logL_X)                  # Likelihood ratio
            for z in range(DREAMPar['N']):
                if logPR_Xp[z] >= logPR_X[z]:               # If proposal better than current chain state
                    if logPR_X[z] < 0:                      # If current state is outside epsilon
                        accept[z] = 1                       # Accept proposal
                    else:                                   # If current state is inside epsilon
                        accept[z] = a_L[z] >= Z[z]          # Accept with Metropolis probability
                else:                                       # If proposal worse than current chain state
                    if logPR_Xp[z] < 0:                     # If proposal is outside epsilon
                        accept[z] = 0                       # Reject proposal
                    else:                                   # If proposal inside epsilon
                        accept[z] = a_L[z] >= Z[z]          # Accept with Metropolis probability

    elif options['ABC'] == 'yes':                           # ABC or Limits of Acceptability
        if DREAMPar['lik'] == 21:                           # Turner approach
            #a_PR = np.exp(logPR_Xp - logPR_X)              # Prior ratio
            #a_S = np.exp(log_alfa_sn)                      # Snooker correction
            #alfa = a_S * a_L * a_PR                        # Product of ratios
            alfa = np.exp(np.clip((logPR_Xp - logPR_X) + (logL_Xp - logL_X) + log_alfa_sn, min_exp_value, max_exp_value))
            accept = alfa >= Z                              # Accept if alfa ≥ Z
        elif DREAMPar['lik'] == 22:     ## Sadegh and Vrugt (2014) approach
            accept = (logL_Xp >= logL_X) | (logL_Xp >= 0)   # Epsilon multiple values
        elif DREAMPar['lik'] == 23:     ## Limits of acceptability (Vrugt and Beven, 2018)
            accept = (logL_Xp >= logL_X) | (logL_Xp == Meas_info['n_S'])

    else:
        raise ValueError("DREAM-Suite ERROR:Metropolis_rule: Unknown option")

    id_acc = np.where(accept > 0)[0]    ## Indices accepted proposals (i.e., row numbers)

    return accept, id_acc


def Evaluate_target(X, DREAMPar, func_handle, Meas_info, options, base_dir, plugin, printed_warnings, verbose = 0):
    ## ##################################################################### ##
    ## Metropolis rule for acceptance or rejection of candidate points       ##
    ##                                                                       ##
    ## SYNOPSIS: [FX, S] = Evaluate_target(X, DREAMPar, func_handle, ...     ##
    ##                                     Meas_info, options, base_dir, ... ##
    ##                                     plugin, printed_warnings, ...     ## 
    ##                                     verbose)                          ##
    ##                                                                       ##
    ##   X          [input] Nxd matrix of candidate points                   ##
    ##   DREAMPar   [input] Dictionary containing algorithmic parameters     ##
    ##   func_handle[input] Function handle passed to worker                 ##
    ##   Meas_info  [input] Dictionary of measurement data [fitting]         ##
    ##   options    [input] Dictionary of algorithmic settings               ##
    ##   base_dir   [input] String with location of model files              ##
    ##   plugin     [input] Variable/dictionary passed to function_handle    ##
    ##   prwarnings [input] List of past printed warnings                    ##
    ##   verbose    [input] Print simulation progress to screen or not       ##
    ##   FX         [outpt] Nxn matrix of simulations candidate points       ##
    ##   S          [outpt] Nxs matrix of summary metrics candidate points   ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, March 2005                              ##
    ## Los Alamos National Laboratory 			                             ##
    ##                                                                       ##
    ## ##################################################################### ##

    if X.ndim == 1:
        X = X[:, None]  # So we can refer to X as X[ii, :] rather than X[ii]

    N = X.shape[0]
    n = Meas_info['n'] + Meas_info['n_S']       # Total number of Y and S

    if DREAMPar['CPU'] == 1:    ## Sequential evaluation - CORRECT in Python
        for ii in range(N):
            if plugin is None:
                results = func_handle(X[ii, :])
            else:
                results = func_handle(X[ii, :], plugin)

            # results can consists or more than one array depending on # return arguments func_handle
            if isinstance(results, (tuple, list)):
                if ii == 0:
                    Z = results[0]
                    if isinstance(Z, (int, float)):
                        ndim = 1
                    elif isinstance(Z, np.ndarray):
                        ndim = Z.shape[0]
                        
                    FX = np.empty((ndim,N))                 # Initialize FX        
                    FX[:,0] = results[0]
                else:
                    FX[:,ii] = results[0]
                # Now check whether we have to store model output
                if options['modout'] == 'yes':
                    if ii == 0:                             # Initialize Y on first iteration
                        Z = results[1]                      # Extract 2nd element
                        if isinstance(Z, (int, float)):
                            ndim = 1
                        elif isinstance(Z, np.ndarray):                     
                            ndim = Z.shape[0]
                        Y = np.empty((N, ndim))             # Initialize Y 
                        Y[0, :] = Z                         # First row of Y
                    else:
                        Y[ii, :] = results[1]               # Subsequent iterations
                else:
                    Y = None                                # Y is returned but modout is 'no'
            else:
                if ii == 0:
                    if isinstance(results, (int, float)):
                        ndim = 1
                    elif isinstance(results, np.ndarray):                     
                        ndim = results.shape[0]

                    FX = np.empty((ndim,N))                # Initialize FX        
                    FX[:,0] = results
                else:
                    FX[:,ii] = results
                if options['modout'] == 'yes':
                    warning_msg = f"Evaluate_target WARNING: Did not expect one output argument from the function as options['modout'] == 'yes'. Setting Y to None"
                    if warning_msg not in printed_warnings:
                        print(warning_msg)
                        printed_warnings.add(warning_msg)                    
                Y = None
 
            if verbose:
                # Print progress if verbose flag is set
                print(f'Posterior simulation, %% done: {100 * ((ii+1) / N):3.2f}', end = '\r')

    elif DREAMPar['CPU'] > 1:   ## Parallel evaluation - CORRECT in Python if IO = 'No' and 'Yes'
        task_ranges = distribute_tasks(N, DREAMPar['CPU'])      # Task distribution for each worker (divide work)
        # task_ranges = [(i * (N // DREAMPar['CPU']), (i + 1) * (N // DREAMPar['CPU'])) for i in range(DREAMPar['CPU'])]
        if isinstance(plugin, dict):                            # Convert to regular numpy array so that it can be shared with workers [= pickable]
            plugin = convert_memoryview_to_array(plugin)
            
        with mp.Pool(processes=DREAMPar['CPU']) as pool:
            results = pool.starmap(worker_task, [(worker_id, start_idx, end_idx, X, func_handle, plugin, base_dir)
                                                    for worker_id, (start_idx, end_idx) in enumerate(task_ranges)])

        if isinstance(results, (tuple, list)):      
            # Initialize FX and set Y to be None
            FX = np.full((n,N),np.nan)
            Y = None
            # Unpack the results from the workers
            ct = 0
            i = 0
            for worker_result in results:
                sample_worker = worker_result[0]
                # print(type(sample_worker))
                if isinstance(sample_worker, (int, float)):
                    nvar = 1
                elif isinstance(sample_worker, np.ndarray):
                    nvar = 1
                elif isinstance(sample_worker, (tuple)):
                    nvar = len(sample_worker)
                if isinstance(worker_result, list):         ## This is when N > DREAMPar['CPU']
                    if nvar == 2:    ## 2 outputs, FX and Y    
                        for fx, y in worker_result:
                            FX[:, ct] = fx.flatten() 
                            if options['modout'] == 'no':
                                warning_msg = f"Evaluate_target WARNING: Did not expect two output arguments from the function. Setting Y to None as options['modout'] == 'no'."
                                if warning_msg not in printed_warnings:
                                    print(warning_msg)
                                    printed_warnings.add(warning_msg)
                            else:
                                if ct == 0:
                                    Y = np.full((N,len(y)),np.nan)
                                Y[ct, :] = y.flatten()
                            # update counter                                
                            ct = ct + 1
    
                    elif nvar == 1:  ## 1 output, FX
                        for fx in worker_result:
                            FX[:, ct] = fx.flatten()
                            ct = ct + 1
                            if options['modout'] == 'yes':
                                warning_msg = f"Evaluate_target WARNING: Did not expect one output argument from the function as options['modout'] == 'yes'. Setting Y to None"
                                if warning_msg not in printed_warnings:
                                    print(warning_msg)
                                    printed_warnings.add(warning_msg)
    
                    else:           ## 0 or more than 2 outputs
                        warning_msg = f"Evaluate_target ERROR: No output arguments from the function. Setting FX and Y to None."
                        if warning_msg not in printed_warnings:
                            print(warning_msg)
                            printed_warnings.add(warning_msg)
 
                elif isinstance(worker_result, tuple):      
                    if nvar == 2:    ## 2 outputs, FX and Y
                        fx, y = worker_result
                        if options['modout'] == 'no':
                            warning_msg = f"Evaluate_target WARNING: Did not expect two output arguments from the function. Setting Y to None as options['modout'] == 'no'."
                            y = []
                            if warning_msg not in printed_warnings:
                                print(warning_msg)
                                printed_warnings.add(warning_msg)
 
                    elif nvar == 1:  ## 1 output, FX
                        fx = worker_result
                        y = []
                        if options['modout'] == 'yes':
                            warning_msg = f"Evaluate_target WARNING: Did not expect one output argument from the function as options['modout'] == 'yes'. Setting Y to None"
                            if warning_msg not in printed_warnings:
                                print(warning_msg)
                                printed_warnings.add(warning_msg)

                    else:           ## 0 or more than 2 outputs
                        warning_msg = f"Evaluate_target ERROR: No output arguments from the function. Setting FX and Y to None."
                        if warning_msg not in printed_warnings:
                            print(warning_msg)
                            printed_warnings.add(warning_msg)

                    if options['modout'] == 'yes':
                        FX[:, ct] = fx.flatten()
                        if ct == 0:
                            Y = np.full((N,len(y)),np.nan)
                        Y[ct, :] = y.flatten()
                    else:
                        FX[:, ct] = fx.flatten()
                    # update counter        
                    ct = ct + 1  

        # If there is only one return value, unpack it directly
        elif not isinstance(results, (tuple, list)):
            FX = results
            Y = None

    # Check if diagnostic Bayes is used
    if options['DB'] == 'yes':
        S = FX[Meas_info['n']: , :]     # Extract summary metrics
        FX = FX[:Meas_info['n'], :]     # Keep only the first part (model output)
    else:
        S = None
        
    if verbose:
        print("\nModel simulation ... done")
    
    return FX, S


def Calc_pCR(DREAMPar, sdX_Xp, TsdX_Xp, cCR, CR):
    ## ##################################################################### ##
    ## This function updates selection probabilities of nCR crossover values ##
    ##                                                                       ##
    ## SYNOPSIS: [TsdX_Xp, cCR, pCR] = Calc_pCR(DREAMPar, sdX_Xp, ...        ##
    ##                                          TsdX_Xp, cCR, CR)            ##
    ##                                                                       ##
    ##   DREAMPar   [input] Dictionary containing algorithmic parameters     ##
    ##   sdX_Xp     [input] Array of standardized Euclidean distances        ##
    ##   TsdX_Xp    [input] 1D array of total traveled distances each CR     ##
    ##   cCR        [input] 1D array of # times each crossover has been used ##
    ##   CR         [input] 2D array of crossover values                     ##
    ##   TsdX_Xp    [outpt] Update traveled standardized Euclidean distances ##
    ##   cCR        [outpt] Update of counts # times crossover values used   ##
    ##   pCR        [outpt] 1D array of crossover selection probabilities    ##
    ##                                                                       ##
    ## © Written by Jasper A. Vrugt, March 2005                              ##
    ## Los Alamos National Laboratory 			                             ##
    ##                                                                       ##
    ## ##################################################################### ##
    
    # Python: CR = 1,...,nCR, whereas index of TsdX_Xp and cCR starts at zero
    # Loop through each individual (DREAMPar.N) and each step (DREAMPar.steps)
    for i in range(DREAMPar['N']):
        for j in range(DREAMPar['steps']):
            # Update the traveled standardized Euclidean distance for the crossover
            TsdX_Xp[CR[i, j]-1] += sdX_Xp[i, j]
            # Count how many times the crossover has been used
            cCR[CR[i, j]-1] += 1
    
    # Calculate the new selection probability for each crossover
    pCR = TsdX_Xp / cCR     # Selection probabilities
    pCR /= np.sum(pCR)      # Normalize probabilities so that sum pCR equals 1
    
    return TsdX_Xp, cCR, pCR


def checkfile_T(file_name):
    # ####################################################################### #
    # Check the content of the restart budget file 'T.txt'                    #
    # ####################################################################### #

    # Check if file exists
    if not os.path.exists(file_name):
        raise FileNotFoundError(f"DREAM-Suite ERROR: File '{file_name}' does not exist.")

    # Read the content of the file
    with open(file_name, 'r') as f:
        content = f.read().strip()

    # Check if the content is empty
    if not content:
        raise ValueError(f"DREAM-Suite ERROR: File '{file_name}' is empty --> Please store an integer in file '{file_name}'.")

    # Try to convert the content to a numeric value
    try:
        T_new = int(content)
    except ValueError:
        raise ValueError(f"DREAM-Suite ERROR: File '{file_name}' does not store a numerical value --> Store only a single value (integer) in file '{file_name}'.")

    # Check if the value is a positive integer
    if T_new <= 0:
        raise ValueError(f"DREAM-Suite ERROR: File '{file_name}' stores negative integers (or zero) --> Store a single positive integer in file '{file_name}'.")

    return T_new


def prepare_output(DREAMPar, Par_info, Meas_info, n_chains, MAP_info, method):
    # ####################################################################### #
    # This function prepares variables for postprocessor DREAM_Suite       #
    # ####################################################################### #

    str_par = []
    if 'names' in Par_info:
        for name in Par_info['names']:
            str_par.append(f"${name}$")
        str_par = [s.replace("{\rm ", "\\text{").replace("\\;", "") for s in str_par]
        # Now create parameter strings for the Tables
        str_table = [s.replace('$', '').replace('\\text','').replace('\\','') for s in str_par]
    else:
        for i in range(DREAMPar['d']):
            str_par.append(f"$x_{{{i+1}}}$")
        # Prepare table string (for output formatting)
        str_table = [s.replace('$', '').replace('\\;', '').replace('\\rm', '') for s in str_par]

    # Create legend/label string for different crossover values
    str_CR = [f"$\\;n_{{\\rm CR}} = {i / DREAMPar['nCR']:4.2f}$" for i in range(1, DREAMPar['nCR'] + 1)]
    # chain string
    str_chain = ['chain 1']
    for i in range(1,n_chains):
        str_chain.append(f'chain {i + 1}')
    str_chain.append('MAP')

    # Create legend/label string for different summary statistics
    str_S = [f"$S_{{{i+1}}}$" for i in range(0, Meas_info['n_S'])]

    # Rename method for correct subscript printing
    method_table = method
    ii = method.find('_')
    if ii != -1:
        method_fig = f"{method[:ii].upper()}$_{{\\rm ({method[ii+1:].upper()})}}$"
    else:
        method_fig = f"{method.upper()}"

    # Look at MAP_info to see whether we use sandwich correction or not
    sndwch = 1 if 'map' in MAP_info else 0

    # Switch based on sndwch
    if sndwch == 0:
        sndwch_text = ': Sandwich Inactive'
    elif sndwch == 1:
        sndwch_text = ': Sandwich Active'

    return str_par, str_table, str_CR, str_chain, str_S, method_table, method_fig, sndwch, sndwch_text 


def tabulate_output(method, DREAMPar, ML, MAP, MEAN, MED, STD, CORR, str_table, sndwch):
    # ####################################################################### #
    # This function tabulates posterior moments of sampled chains             #
    # ####################################################################### #

    method = method.upper()

    # Open the output file for writing
    filename = f"{method}_output.txt"
    
    max_length = max(len(name) for name in str_table)   # maximum length for table writing  
    max_length = max(5,max_length)                      # must at least be 5 characters
    # Write the results to an output file
    with open(filename, 'w') as fid:
        fid.write('--------------------------- DREAM-Suite output file --------------------------- \n')
        fid.write('\n')
        if sndwch == 0:
            fid.write('Sandwich correction is INACTIVE \n')
        elif sndwch == 1:
            fid.write('Sandwich correction is ACTIVE \n')
        fid.write('\n')
        fid.write(f" Table 1. {method}: Maximum likelihood and posterior median, mean, and standard deviation of parameters \n")
        fid.write('          ===================================================================== \n')
        fid.write('           Parameter         ML        MAP       MEAN      MEDIAN       STD     \n')
        fid.write('          --------------------------------------------------------------------- \n')
        # Print parameter statistics
        fmt_1 = f'          %-{max_length}s \t %7.3f    %7.3f    %7.3f    %7.3f    %7.3f\n'
        for j in range(DREAMPar['d']):
            fid.write(fmt_1 % (str_table[j], ML[j], MAP[j], MED[j], MEAN[j], STD[j]))
        fid.write('          ===================================================================== \n')
        fid.write("          ML: Maximum Likelihood values                                         \n")
        fid.write("          MAP: Maximum A-Posteriori density values                              \n")
        fid.write("          MEAN: Posterior mean values                                           \n")
        fid.write("          STD: Posterior standard deviation                                     \n")
        fid.write('\n')
        fid.write('\n')
        
        table_width = (max_length+2) * ( DREAMPar['d'] + 1)     # Each column is 7 characters wide
        top_line = '=' * table_width                            # Line of '=' with the correct length
        # Print correlation matrix
        fid.write(f" Table 2. {method}: Pearson''s correlation coefficients of posterior parameter values \n")
        fid.write(f"          {top_line} \n")
        fid.write(f"{'':{max_length+15}}")
        for i in range(DREAMPar['d']):
            fid.write(f"{str_table[i]:{max_length+2}}")
        fid.write('\n')
        for i in range(DREAMPar['d']):
            fid.write(f"          {str_table[i]:{max_length+2}}")
            for j in range(DREAMPar['d']):
                if DREAMPar['d'] == 1:
                    fid.write(f"{CORR:{max_length+2}.3f}")
                else:
                    fid.write(f"{CORR[i, j]:{max_length+2}.3f}")
            fid.write('\n')
        fid.write(f"          {top_line} \n")
        fid.write('\n')
        fid.write('------------------------- End DREAM-Suite output file ------------------------- \n')

    # Optionally open the file (not on Linux/Unix)
    if os.name != 'posix':  # Not on Unix-like systems
        os.startfile(filename)


def tabulate_diagnostics(method, DREAMPar, options, chain, output, iloc, sndwch):
    # ####################################################################### #
    # This function tabulates within-chain convergence diagnostics            #
    # ####################################################################### #

    # Open the file for writing diagnostics
    filename = f"{method}_diagnostics.txt"

    # 1. Open warning file
    with open('warning_file.txt', 'a+') as fid_w:
        if sndwch == 0:
            fid_w.write('Sandwich correction is INACTIVE\n')
        else:
            fid_w.write('Sandwich correction is ACTIVE\n')

        # Tabulate within-chain convergence diagnostics [only if chain has at least 100 samples]
        if options['diagnostics'] == 'yes':
            if iloc > 200:
                with open(filename, 'w') as fid_d:
                    # Print header to file
                    evalstr_file = f"{method} diagnostics file"
                    ii_string = len(evalstr_file)
                    d_ii = 76 - ii_string - 2
                    d_half = d_ii / 2
                    d1 = int(np.floor(d_half))
                    d2 = int(d_ii - d1)
                    plot_str = '=' * (d1 - 1) + f" {evalstr_file} " + '=' * (d2 - 1)
                    fid_d.write(plot_str + '\n')
                    # Sandwich correction status
                    if sndwch == 0:
                        fid_d.write('Sandwich correction is INACTIVE\n')
                    elif sndwch == 1:
                        fid_d.write('Sandwich correction is ACTIVE\n')
                    
                    fid_d.write('\n')
                    fid_d.write(f" Table 3. {method}: Single chain convergence diagnostics \n")
                    # Calculate convergence diagnostics for individual chains
                    for j in range(DREAMPar['N']):
                        # Calculate diagnostics using coda (assumed function)
                        diagnostic_info = coda(chain[int(0.5*iloc):iloc, :DREAMPar['d'], j])             # Slice the chain data
                        diagnostic_info['chain_number'] = j + 1
                        prt_coda(diagnostic_info, None, fid_d)
                    # End of file footer
                    d_ii = 76 - ii_string - 4
                    d_half = d_ii / 2
                    d1 = int(np.floor(d_half))
                    d2 = int(d_ii - d1)
                    plot_str = '=' * (d1 - 2) + f" End {evalstr_file} " + '=' * (d2 - 2)
                    fid_d.write(plot_str + '\n')
                    flag = 1
            else:
                evalstr = f'DREAM-Suite WARNING: Cannot compute coda diagnostics for chains due to insufficient chain samples -> Use at least DREAMPar.T = {(200 * DREAMPar["thinning"])}\n'
                print(evalstr)
                fid_w.write(evalstr)
                flag = 0
        # 4. Check convergence
        if not np.all(output['R_stat'][-1, 1:DREAMPar['d'] + 1] < 1.2):
            evalstr = 'DREAM-Suite WARNING: Chains did not converge according to \\hat{R} scale reduction factor\n'
            print(evalstr)
            fid_w.write(evalstr)
        # 5. Write final lines of warning file
        fid_w.write('----------- End of DREAM-Suite WARNING file ---------\n')
  
    # Optionally open the file on screen (works on Windows and macOS, not Linux/Unix)
    if os.name != 'posix' and flag == 1:  # Not on Unix-like systems
        os.startfile(filename)


def define_lik(DREAMPar, Par_info, par_names, fpar_mod, parmin_mod, parmax_mod):
    # ####################################################################### #
    # This function defines the likelihood function (variables, names, etc.)  #
    # ####################################################################### #

    # Handling likelihood function selection with a switch-like structure
    if DREAMPar['lik'] == 13:  # Normal distribution
        # index:      0   1    2   3
        # parname:   s0  s1  phi1 phi2    
        fpar_nuis = [1e-4, 0, 0, 0]
        parmin_nuis = [0, 0, 0, 0]
        parmax_nuis = [1, 1, 1, 1]
        filename = 'Normal'
        id_nuis = [0, 1, 2, 3]
    elif DREAMPar['lik'] == 14:  # GL (OBSOLETE)
        # index:      0   1    2    3   4   5   6    7    8    9  10  
        # parname:  std0 std1 beta xi  mu1 phi1 phi2 phi3 phi4 K lambda
        fpar_nuis = [0.1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
        parmin_nuis = [0, 0, -1, 0.1, 0, 0, 0, 0, 0, 0, 0.1]
        parmax_nuis = [1, 1, 1, 10, 100, 1, 1, 1, 1, 1, 1]
        filename = 'GL'
        # id_nuis = [0, 1, 2, 3, 5]
    elif DREAMPar['lik'] == 16:  # Laplace distribution
        # index:      0   1   2
        # parname:   s0  s1  phi1
        fpar_nuis = [1e-4, 0, 0]
        parmin_nuis = [0, 0, 0]
        parmax_nuis = [1, 1, 1]
        filename = 'Laplace'
        id_nuis = [0, 1, 2]
    elif DREAMPar['lik'] == 17:  # SST (Skewed Student-t)
        # index:      0   1   2    3   4    5
        # parname:   s0  s1  nu   xi phi1  phi2
        fpar_nuis = [1e-4, 0, 1e10, 1, 0, 0]
        parmin_nuis = [0, 0, 2, 0.1, 0, 0]
        parmax_nuis = [1, 1, 100, 10, 1, 1]
        filename = 'SL'
        id_nuis = [2, 3, 4]
    elif DREAMPar['lik'] == 44:  # GL+
        # index:      0   1   2    3   4    5
        # parname:   s0  s1  beta xi phi1  phi2
        fpar_nuis = [1e-4, 0, 0, 1, 0, 0]
        parmin_nuis = [0, 0, -1, 0.1, 0, 0]
        parmax_nuis = [1, 1, 1, 10, 1, 1]
        filename = 'GL_plus'
        id_nuis = [2, 3, 4]
    elif DREAMPar['lik'] == 45:  # Universal likelihood function
        # index:      0   1    2    3   4   5    6   
        # parname:    s0  s1  labda p   q  phi1 phi2
        fpar_nuis = [1e-4, 0, 0, 2, 1e10, 0, 0]
        parmin_nuis = [0, 0, -1, 0.5, 2, 0, 0]
        parmax_nuis = [1, 1, 1, 100, 100, 1, 1]
        filename = 'UL'
        id_nuis = [2, 3, 4, 5]
    else:
        idx_nuis = fname = fpar_nuis = parmin_nuis = parmax_nuis = None

    # Define nuisance variable names
    nuis_names = nuis_var_names(DREAMPar)
    names = par_names + nuis_names

    # Model parameters and nuisance variables
    n_modpar = len(fpar_mod)
    id_vpar = list(range(n_modpar)) + [x + n_modpar for x in id_nuis]
    fpar = fpar_mod + fpar_nuis
    parmin = parmin_mod + parmin_nuis
    parmax = parmax_mod + parmax_nuis

    # Min/max values of parameter selection
    Par_info['min'] = [parmin[i] for i in id_vpar]
    Par_info['max'] = [parmax[i] for i in id_vpar]
    Par_info['names'] = [names[i] for i in id_vpar]

    return DREAMPar, Par_info, id_vpar, fpar, filename 


def setup_lik(Func_name, DREAMPar, Par_info, Meas_info, LV, unc_method = 1):
    # ####################################################################### #
    # This function defines Lik_info dictionary with likelihood formulation   #
    # ####################################################################### #

    # Initialization of output variables
    Lik_info = {}

    # Name of likelihood function
    if DREAMPar['lik'] == 1:  # User supplies likelihood
        name_lik = 'Own likelihood'
        name_lik_func = Func_name
    elif DREAMPar['lik'] == 2:  # User supplies log-likelihood
        name_lik = 'Own log-likelihood'
        name_lik_func = Func_name
    elif DREAMPar['lik'] == 11:  # Box and Tiao (1953) variant of log-likelihood
        name_lik = 'Gaussian likelihood without sigma'
        name_lik_func = 'Normal'
    elif DREAMPar['lik'] == 12:  # Normal likelihood with homos/heteroscedastic meas. error
        name_lik = 'Gaussian likelihood'
        name_lik_func = 'Normal'
    elif DREAMPar['lik'] == 13:  # Normal likelihood with AR(2) & (non)constant residual var.
        name_lik = 'Normal likelihood with AR(2)-process'
        name_lik_func = 'NL'
        lik_names = ['s_{0}', 's_{1}', '\\phi_{1}', '\\phi_{2}']
        ar_id = [2, 3]
    elif DREAMPar['lik'] == 14:  # Generalized likelihood with AR(2) & (non)constant res. var.
        name_lik = 'Generalized Likelihood'
        name_lik_func = 'GL'
        lik_names = ['s_{0}', 's_{1}', '\\beta', '\\xi', '\\mu_{1}', '\\phi_{1}', '\\phi_{2}', '\\phi_{3}', '\\phi_{4}', 'K', '\\lambda']
        ar_id = [5, 6, 7, 8]
    elif DREAMPar['lik'] == 15:  # Spectral likelihood of Whittle, 1952, 1953
        name_lik = 'Whittle (spectral) likelihood'
        name_lik_func = 'Whittle'
    elif DREAMPar['lik'] == 16:  # Laplace likelihood with AR(1) & (non)constant residual var.
        name_lik = 'Laplacian likelihood with AR(1) process'
        name_lik_func = 'LAPL'
        lik_names = ['s_{0}', 's_{1}', '\\phi_{1}']
        ar_id = [2]
    elif DREAMPar['lik'] == 17:  # Student t likelihood with AR(2) & (non)constant residual var.
        name_lik = 'Skewed Student T likelihood'
        name_lik_func = 'SL'
        lik_names = ['s_{0}', 's_{1}', '\\nu', '\\xi', '\\phi_{1}', '\\phi_{2}']
        ar_id = [4, 5]
    elif DREAMPar['lik'] == 21:  # Approximate Bayesian Computation: Turner formulation
        name_lik = 'Approximate Bayesian Computation'
        name_lik_func = 'ABC'
    elif DREAMPar['lik'] == 22:  # Approximate Bayesian Computation: Sadegh and Vrugt, 2015
        name_lik = 'Approximate Bayesian Computation - alternative kernel'
        name_lik_func = 'ABC-kernel'
    elif DREAMPar['lik'] == 23:  # Approximate Bayesian Computation: Vrugt and Beven, 2018
        name_lik = 'Approximate Bayesian Computation - limits of acceptability'
        name_lik_func = 'ABC-LOA'
    elif DREAMPar['lik'] == 31:  # Informal likelihood of Beven and Freer, 2001
        name_lik = 'Informal likelihood: Beven and Freer, 2001'
        name_lik_func = 'Informal_BFa'
    elif DREAMPar['lik'] == 32:  # Informal likelihood of Beven and Freer, 2001
        name_lik = 'Informal likelihood: Option b in Beven and Freer, 2001'
        name_lik_func = 'Informal_BFb'
    elif DREAMPar['lik'] == 33:  # Informal likelihood of Beven and Freer, 2001
        name_lik = 'Informal likelihood: Option c in Beven and Freer, 2001'
        name_lik_func = 'Informal_BFc'
    elif DREAMPar['lik'] == 34:  # Informal likelihood of Beven and Binley, 1992
        name_lik = 'Informal likelihood: Last option, Page 284 of Beven and Binley, 1992'
        name_lik_func = 'Informal_BFd'
    elif DREAMPar['lik'] == 44:  # Generalized likelihood ++, AR(2) & (non)constant residual var.
        name_lik = 'Generalized Likelihood ++'
        name_lik_func = 'GL_plus'
        lik_names = ['s_{0}', 's_{1}', '\\beta', '\\xi', '\\phi_{1}', '\\phi_{2}']
        ar_id = [4, 5]
    elif DREAMPar['lik'] == 45:  # Universal likelihood ++, AR(2) & (non)constant residual var.
        name_lik = 'Universal likelihood'
        name_lik_func = 'UL'
        lik_names = ['s_{0}', 's_{1}', '\\lambda', 'p', 'q', '\\phi_{1}', '\\phi_{2}']
        ar_id = [5, 6]
    elif DREAMPar['lik'] == 52:  # Generalized least squares form of likelihood function
        name_lik = 'Matrix implementation of Gaussian likelihood using GLS form'
        name_lik_func = 'Normal_GLS'
    # Power Likelihoods
    elif DREAMPar['lik'] == 61:  # Laplace power likelihood with unit integral
        name_lik = 'Power likelihood: Laplace distribution'
        name_lik_func = 'Power_laplace'
    elif DREAMPar['lik'] == 62:  # Normal power likelihood with unit integral
        name_lik = 'Power likelihood: Normal distribution'
        name_lik_func = 'Power_normal'
    elif DREAMPar['lik'] == 63:  # Laplace power likelihood with fixed lambda
        name_lik = 'Power likelihood: Laplace distribution - fixed lambda'
        name_lik_func = 'Power_laplace'
    elif DREAMPar['lik'] == 64:  # Normal power likelihood with fixed lambda
        name_lik = 'Power likelihood: Normal distribution - fixed lambda'
        name_lik_func = 'Power_normal'

    # Now determine properties of likelihood function
    if DREAMPar['lik'] in [13, 14, 16, 17, 44, 45]:

        nest = len(LV['id_vpar'])               # number of estimable parameters
        ntot = len(LV['fpar'])                  # total number of parameters 
        LV['id_vpar'] = np.array(LV['id_vpar']) # turn field into numpy array
        LV['fpar'] = np.array(LV['fpar'])       # turn field into numpy array
        nf = ntot - nest                        # number of fixed nuisance variables
        nmod = nest - len(lik_names) + nf       # number of model parameters
        nmod = np.array(nmod).astype(int)
        # Check whether s0 and s1 are selected
        if len(LV['id_vpar']) >= nmod + 1:
            #if LV['id_vpar'][nmod] == nmod + 1:
            if LV['id_vpar'][nmod] == nmod:
                s0 = 1
                if len(LV['id_vpar']) >= nmod + 2:
                    # if LV['id_vpar'][nmod + 1] == nmod + 2:
                    if LV['id_vpar'][nmod + 1] == nmod + 1: 
                        s1 = 1
                    else: 
                        s1 = 0
                else:
                    s1 = 0
            # elif LV['id_vpar'][nmod] == nmod + 2:
            elif LV['id_vpar'][nmod] == nmod + 1:            
                s0 = 0
                s1 = 1
            else:
                s0 = 0
                s1 = 0
        else:
            s0 = 0
            s1 = 0
        # Check whether user specified Sigma or not
        # if 'Sigma' in Meas_info and np.size(Meas_info['Sigma']) > 0:
        if Meas_info['Sigma'] != None and np.size(Meas_info['Sigma']) > 0:
            method = 0  # no s0/s1 needed: Sigma is defined by user
        else:
            # Check treatment (method) of measurement error variance
            if Meas_info['sigma2'] == 'constant':
                if s0 == 0:
                    if s1 == 0:
                        method = 1
                    elif s1 == 1:
                        method = 2
                elif s0 == 1:
                    if s1 == 0:
                        method = 3
                    elif s1 == 1:
                        method = 4
            elif Meas_info['sigma2'] == 'nonconstant':
                if s0 == 0:
                    if s1 == 0:
                        method = 5
                    elif s1 == 1:
                        method = 6
                elif s0 == 1:
                    if s1 == 0:
                        method = 7
                    elif s1 == 1:
                        method = 8

        # Initialize the remaining variables
        id_rem = np.array([]).astype(int)
        if method == 0:
            if s0 == 0:
                if s1 == 1:
                    id_rem = [(nmod)] # [nmod + 1]
            elif s0 == 1:
                id_rem = [(nmod)] # [nmod + 1]
                if s1 == 1:
                    id_rem.append((nmod + 1))
        elif method == 2:
            id_rem = [(nmod)] # [nmod + 1]
        elif method == 4:
            id_rem = [(nmod + 1)] # [nmod + 2]
        # Remove id_rem entry from variable parameters & DREAMPar
        LV['id_vpar'] = np.delete(LV['id_vpar'], id_rem)
        DREAMPar['d'] -= len(id_rem)
        # Remove from Par_info
        Par_info['min'] = np.delete(Par_info['min'], id_rem)
        Par_info['max'] = np.delete(Par_info['max'], id_rem)
        if 'names' in Par_info:
            Par_info['names'] = np.delete(Par_info['names'], id_rem)
        if 'steps' in Par_info:
            Par_info['steps'] = np.delete(Par_info['steps'], id_rem)
            Par_info['step_size'] = np.delete(Par_info['step_size'], id_rem)

        # Process and generate the likelihood strings and nuisance variables
        fpar = LV['fpar']
        id_vpar = LV['id_vpar']
        par = np.array(fpar.copy())
        par[id_vpar] = np.nan
        nest = len(id_vpar)                             # Number of estimable parameters
        nf = ntot - nest                                # Number of fixed nuisance variables
        id_nuisvar = range(nmod, ntot)                  # Index of nuisance variables
        nuisvar = par[id_nuisvar]                       # Isolate nuisance variables
        id_ar = nmod + ar_id                            # Index of AR coefficients
        id_par = np.zeros(ntot)                         # Vector of variable parameters
        id_par[id_vpar] = 1                             # Set variable parameters to 1
#        t1 = 1 + np.sum(id_par[id_ar])                  # 1st measurement AR model
        t1 = np.sum(id_par[id_ar])                      # Python: 1st measurement AR model
        lik_names = np.array(lik_names)
        str_nuis = lik_names[np.isnan(nuisvar)]         # Get the likelihood parameters
        str_nuis = ['$' + s + '$' for s in str_nuis]    # Add $$ signs in front and end of str_lik for LaTeX print

        # Convert likelihood strings
        lik_strL = f"loglik, std_e, eps_n, f_eps_n, _ = {name_lik_func}('est', nuisvar, {method}, FXp[0:n, ii], Meas_info['Y'], Meas_info['Sigma'])"

        if unc_method == 1:
            lik_strF = f"_, _, _, _, Y_r = {name_lik_func}('sim', nuisvar, {method}, Ufx.T, Meas_info['Y'], Meas_info['Sigma'], Nr)"  
        else:
            lik_strF = f"_, _, _, _, Y_r = {name_lik_func}('sim', nuisvar, {method}, Ufx[ii, 0:Meas_info['n']].T, Meas_info['Y'], Meas_info['Sigma'], Nr[ii])"  

        if 'filename' in LV:
            filename = 1
        else:
            filename = None

        Lik_info = {'name_lik': name_lik,
            'name_lik_func': name_lik_func,
            'fpar': fpar,
            'id_vpar': id_vpar,
            'index': nmod,
            'stringL': lik_strL,
            'stringF': lik_strF,
            'str_nuis': str_nuis,
            't_start': int(t1),
            'id_ar': id_ar,
            'id_nuisvar': id_nuisvar,
            'method': method,
            'filename': filename}

    else: # other likelihood functions

        Lik_info = {
            'name_lik': name_lik,
            'name_lik_func': name_lik_func,
            'fpar': None,
            'id_vpar': None,
            'index': (0),
            'stringL': '',
            'stringF': '',
            'str_nuis': None,
            't_start': (0),
            'id_ar': None,
            'id_nuisvar': None,
            'method': (0),
            'filename': (0)
        }

    return Lik_info, DREAMPar, Par_info


def nuis_var_names(DREAMPar):
    # ####################################################################### #
    # This function returns nuisance variables chosen likelihood formulation  #
    # ####################################################################### #

    nuis_names = []  # Initialize empty list
    
    # Switch-like structure using if-elif statements
    if DREAMPar['lik'] == 13:
        # Normal AR(2)-likelihood with (non)constant residual variance
        nuis_names = ['s_{0}', 's_{1}', '\\phi_{1}', '\\phi_{2}']
    elif DREAMPar['lik'] == 14:
        # Generalized Likelihood function [= OBSOLETE]
        nuis_names = ['s_{0}', 's_{1}', '\\beta', '\\xi', '\\mu_{1}', 
                      '\\phi_{1}', '\\phi_{2}', '\\phi_{3}', '\\phi_{4}', 'K', '\\lambda']
    elif DREAMPar['lik'] == 16:
        # Laplace AR(1)-likelihood with (non)constant residual variance
        nuis_names = ['s_{0}', 's_{1}', '\\phi_{1}']
    elif DREAMPar['lik'] == 17:
        # Student AR(2)-likelihood with (non)constant residual variance
        nuis_names = ['s_{0}', 's_{1}', '\\nu', '\\xi', '\\phi_{1}', '\\phi_{2}']
    elif DREAMPar['lik'] == 44:
        # Generalized likelihood ++, AR(2) & (non)constant residual var.
        nuis_names = ['s_{0}', 's_{1}', '\\beta', '\\xi', '\\phi_{1}', '\\phi_{2}']
    elif DREAMPar['lik'] == 45:
        # Universal likelihood, AR(2) & (non)constant residual variance
        nuis_names = ['s_{0}', 's_{1}', '\\lambda', 'p', 'q', '\\phi_{1}', '\\phi_{2}']
    
    return nuis_names


def Discrete_space(X, Par_info):
    # ####################################################################### #
    # This function transforms a continuous X-space to a discrete space       #
    # ####################################################################### #
    
    method = 2      # Use the latest MATLAB release method (since the two methods differ only slightly)
    N = X.shape[0]  # Number of candidate vectors

    # Step 1: Transform continuous X to integer between 0 and number of steps
    # Step 2: Back transform to discrete space
    
    if method == 1:
        # Proper method for all MATLAB releases (older version)
        X_min = np.tile(Par_info['min'], (N, 1))  # Replicate min values N times
        X_int = np.round(np.tile(Par_info['steps'], (N, 1)) * ((X - X_min) / np.tile(Par_info['max'] - Par_info['min'], (N, 1))))
        X_dis = X_min + X_int * np.tile(Par_info['step_size'], (N, 1))
    
    elif method == 2:
        # New MATLAB method (for later releases)
        X_int = np.round(Par_info['steps'] * ((X - Par_info['min']) / (Par_info['max'] - Par_info['min'])))
        X_dis = Par_info['min'] + X_int * Par_info['step_size']
    
    return X_dis


def Remove_outlier(method, DREAMPar, X, t, loglik, options):
    # ####################################################################### #
    # This function identifies and removes outlier chains                     #
    # ####################################################################### #

    # Initial flag
    flag = 0
    
    # Check if diagnostic Bayes is used
    if options['DB'] == 'yes':
        logPR = X[:DREAMPar['N'], DREAMPar['d'] + 1]
        if np.all(logPR > 0):  # --> outlier based on likelihood not prior
            y = np.mean(loglik, axis = 0)
            flag = 1
        else:
            y = logPR  # --> outlier first based on prior only
    else:
        y = np.mean(loglik, axis = 0)

    # Choose outlier detection method
    if DREAMPar['outlier'] == 'iqr':
        chain_out = iqr(y)
    elif DREAMPar['outlier'] == 'grubbs':
        chain_out = grubbs(y)
    elif DREAMPar['outlier'] == 'peirce':
        chain_out = peirce(y)
    elif DREAMPar['outlier'] == 'chauvenet':
        chain_out = chauvenet(y)
    
    # Number of outlier chains
    N_out = len(chain_out)
    if N_out > 0:
        outlier = np.column_stack([np.full(N_out, t), chain_out])
        
        # Select good chains to replace outliers
        chain_in = list(range(DREAMPar['N']))
        chain_in = [i for i in chain_in if i not in chain_out]
        chain_select = random.sample(chain_in, N_out)
        for j in range(N_out):
            # Replace loglikelihood of outlier chain with a selected chain
            if options['DB'] == 'no' or flag == 1:
                loglik[:, chain_out[j]] = loglik[:, chain_select[j]]
            
            # Replace the state of outlier chain with selected chain
            X[chain_out[j], :DREAMPar['d'] + 2] = X[chain_select[j], :DREAMPar['d'] + 2]
            
            # Write warning to file
            with open('warning_file.txt', 'a+') as fid:
                fid.write(f"{method} WARNING: Irreversible jump chain {chain_out[j]} at generation {t}\n")
    else:
        outlier = None
    
    return X, loglik, outlier


# Secondary functions for outlier detection
def iqr(data):
    
    Q1, Q3 = np.percentile(data, [75, 25])
    IQR = Q1 - Q3
    return np.where(data < (Q3 - 2 * IQR))[0]


def grubbs(data, alpha=0.05):

    # Number of samples (chains)
    N = len(data)   
    # Calculate Grubbs statistic (for minimum only - one-sided interval)
    G = (np.mean(data) - np.min(data)) / np.std(data)
    # Compute critical t value for one-sided interval (1 - alpha same result!)
    # t_crit = t.ppf(1 - alpha / N, N - 2) ** 2
    t_crit = t.ppf(alpha / N, N - 2) ** 2
    # Now calculate Grubbs critical value
    T_c = (N - 1) / np.sqrt(N) * np.sqrt(t_crit / (N - 2 + t_crit))
    # Check whether to reject null-hypothesis (whether the min is an outlier)
    if G > T_c:
        # Minimum of data is an outlier
        id_outlier = np.argmin(data)
    else:
        id_outlier = None
    
    return id_outlier


def peirce(data):
    # Peirce's table (r values for different sample sizes)
    peirce_r = np.array([
        [-1, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        [3, 1.196, -1, -1, -1, -1, -1, -1, -1, -1],
        [4, 1.383, 1.078, -1, -1, -1, -1, -1, -1, -1],
        [5, 1.509, 1.200, -1, -1, -1, -1, -1, -1, -1],
        [6, 1.610, 1.299, 1.099, -1, -1, -1, -1, -1, -1],
        [7, 1.693, 1.382, 1.187, 1.022, -1, -1, -1, -1, -1],
        [8, 1.763, 1.453, 1.261, 1.109, -1, -1, -1, -1, -1],
        [9, 1.824, 1.515, 1.324, 1.178, 1.045, -1, -1, -1, -1],
        [10, 1.878, 1.570, 1.380, 1.237, 1.114, -1, -1, -1, -1],
        [11, 1.925, 1.619, 1.430, 1.289, 1.172, 1.059, -1, -1, -1],
        [12, 1.969, 1.663, 1.475, 1.336, 1.221, 1.118, 1.009, -1, -1],
        [13, 2.007, 1.704, 1.516, 1.379, 1.266, 1.167, 1.070, -1, -1],
        [14, 2.043, 1.741, 1.554, 1.417, 1.307, 1.210, 1.120, 1.026, -1],
        [15, 2.076, 1.775, 1.589, 1.453, 1.344, 1.249, 1.164, 1.078, -1],
        [16, 2.106, 1.807, 1.622, 1.486, 1.378, 1.285, 1.202, 1.122, 1.039],
        [17, 2.134, 1.836, 1.652, 1.517, 1.409, 1.318, 1.237, 1.161, 1.084],
        [18, 2.161, 1.864, 1.680, 1.546, 1.438, 1.348, 1.268, 1.195, 1.123],
        [19, 2.185, 1.890, 1.707, 1.573, 1.466, 1.377, 1.298, 1.226, 1.158],
        [20, 2.209, 1.914, 1.732, 1.599, 1.492, 1.404, 1.326, 1.255, 1.190],
        [21, 2.230, 1.938, 1.756, 1.623, 1.517, 1.429, 1.352, 1.282, 1.218],
        [22, 2.251, 1.960, 1.779, 1.646, 1.540, 1.452, 1.376, 1.308, 1.245],
        [23, 2.271, 1.981, 1.800, 1.668, 1.563, 1.475, 1.399, 1.332, 1.270],
        [24, 2.290, 2.000, 1.821, 1.689, 1.584, 1.497, 1.421, 1.354, 1.293],
        [25, 2.307, 2.019, 1.840, 1.709, 1.604, 1.517, 1.442, 1.375, 1.315],
        [26, 2.324, 2.037, 1.859, 1.728, 1.624, 1.537, 1.462, 1.396, 1.336],
        [27, 2.341, 2.055, 1.877, 1.746, 1.642, 1.556, 1.481, 1.415, 1.356],
        [28, 2.356, 2.071, 1.894, 1.764, 1.660, 1.574, 1.500, 1.434, 1.375],
        [29, 2.371, 2.088, 1.911, 1.781, 1.677, 1.591, 1.517, 1.452, 1.393],
        [30, 2.385, 2.103, 1.927, 1.797, 1.694, 1.608, 1.534, 1.469, 1.411],
        [31, 2.399, 2.118, 1.942, 1.812, 1.710, 1.624, 1.550, 1.486, 1.428],
        [32, 2.412, 2.132, 1.957, 1.828, 1.725, 1.640, 1.567, 1.502, 1.444],
        [33, 2.425, 2.146, 1.971, 1.842, 1.740, 1.655, 1.582, 1.517, 1.459],
        [34, 2.438, 2.159, 1.985, 1.856, 1.754, 1.669, 1.597, 1.532, 1.475],
        [35, 2.450, 2.172, 1.998, 1.870, 1.768, 1.683, 1.611, 1.547, 1.489],
        [36, 2.461, 2.184, 2.011, 1.883, 1.782, 1.697, 1.624, 1.561, 1.504],
        [37, 2.472, 2.196, 2.024, 1.896, 1.795, 1.711, 1.638, 1.574, 1.517],
        [38, 2.483, 2.208, 2.036, 1.909, 1.807, 1.723, 1.651, 1.587, 1.531],
        [39, 2.494, 2.219, 2.047, 1.921, 1.820, 1.736, 1.664, 1.600, 1.544],
        [40, 2.504, 2.230, 2.059, 1.932, 1.832, 1.748, 1.676, 1.613, 1.556],
        [41, 2.514, 2.241, 2.070, 1.944, 1.843, 1.760, 1.688, 1.625, 1.568],
        [42, 2.524, 2.251, 2.081, 1.955, 1.855, 1.771, 1.699, 1.636, 1.580],
        [43, 2.533, 2.261, 2.092, 1.966, 1.866, 1.783, 1.711, 1.648, 1.592],
        [44, 2.542, 2.271, 2.102, 1.976, 1.876, 1.794, 1.722, 1.659, 1.603],
        [45, 2.551, 2.281, 2.112, 1.987, 1.887, 1.804, 1.733, 1.670, 1.614],
        [46, 2.560, 2.290, 2.122, 1.997, 1.897, 1.815, 1.743, 1.681, 1.625],
        [47, 2.568, 2.299, 2.131, 2.006, 1.907, 1.825, 1.754, 1.691, 1.636],
        [48, 2.577, 2.308, 2.140, 2.016, 1.917, 1.835, 1.764, 1.701, 1.646],
        [49, 2.585, 2.317, 2.149, 2.026, 1.927, 1.845, 1.775, 1.711, 1.657],
        [50, 2.593, 2.326, 2.157, 2.035, 1.937, 1.855, 1.785, 1.721, 1.667],
	    [51, 2.600, 2.334, 2.167, 2.044, 1.945, 1.863, 1.792, 1.730, 1.675],
    	[52, 2.608, 2.342, 2.175, 2.052, 1.954, 1.872, 1.802, 1.740, 1.685],
    	[53, 2.615, 2.350, 2.184, 2.061, 1.963, 1.881, 1.811, 1.749, 1.694],
    	[54, 2.622, 2.358, 2.192, 2.069, 1.972, 1.890, 1.820, 1.758, 1.703],
    	[55, 2.629, 2.365, 2.200, 2.077, 1.980, 1.898, 1.828, 1.767, 1.711],
    	[56, 2.636, 2.373, 2.207, 2.085, 1.988, 1.907, 1.837, 1.775, 1.720],
    	[57, 2.643, 2.380, 2.215, 2.093, 1.996, 1.915, 1.845, 1.784, 1.729],
    	[58, 2.650, 2.387, 2.223, 2.101, 2.004, 1.923, 1.853, 1.792, 1.737],
    	[59, 2.656, 2.394, 2.230, 2.109, 2.012, 1.931, 1.861, 1.800, 1.745],
    	[60, 2.663, 2.401, 2.237, 2.116, 2.019, 1.939, 1.869, 1.808, 1.753] ])

    # Number of samples (chains)
    N = len(data)
    # Find the row index to use in the table for this sample
    if 2 < N < 61:
        n_ind = np.where(peirce_r[:, 0] == N)[0][0]
    else:
        if N >= 61:
            print("WARNING: DREAMPar.N > 60; using Peirce r-values for N = 60")
            # We continue with N = 60 (last row of peirce_r)
            n_ind = peirce_r.shape[0] - 1

    # Find the current r value
    r_curr = peirce_r[n_ind, 1]
    # One-sided interval! (thus negative distance)
    max_neg_dev_allowed = -r_curr * np.std(data)
    # Calculate distance to mean of each data point
    dev_L = data - np.mean(data)
    # Now apply the test (one-sided)
    id_outlier = np.where(dev_L < max_neg_dev_allowed)[0]

    return id_outlier


def chauvenet(data):

    # Number of samples (chains)
    N = len(data)
    # Calculate deviation from mean
    dev_L_ratio = (data - np.mean(data)) / np.std(data)
    # Define table with critical deviations
    n_sample = np.array([3, 4, 5, 6, 7, 10, 15, 25, 50, 100, 300, 500, 1000])
    max_dev_ratio = np.array([1.38, 1.54, 1.65, 1.73, 1.80, 1.96, 2.13, 2.33, 2.57, 2.81, 3.14, 3.29, 3.48])
    # Interpolate (linearly) the max deviation allowable (one-sided & negative)
    max_neg_dev_allowed = -np.interp(N, n_sample, max_dev_ratio)
    # Apply test (one-sided)
    id_outlier = np.where(dev_L_ratio < max_neg_dev_allowed)[0]
    
    return id_outlier


def Gelman(chain, t, method):
    ## ################################################################################## ##
    ## This function computes univariate/multivariate scale reduction factors             ##
    ##                                                                                    ##
    ##  SYNOPSIS: hatR, hatRd = Gelman(chain, t, method)                                  ##
    ##  where                                                                             ##
    ##   chain      [input] nxdxN array of chains                                         ##
    ##   t          [input] Generation [= sample number of MCMC chain]                    ##
    ##   method     [input] MCMC sampling method (= string)                               ##
    ##   hatR       [outpt] n x d matrix of univariate scale-reduction factors parameters ##
    ##   hatRd      [outpt] n x 1 vector of multivariate scale-reduction factors          ##
    ##                                                                                    ##
    ##  Reference:                                                                        ##
    ##   Gelman, A. and D.R. Rubin, (1992) Inference from iterative simulation using      ##
    ##       multiple chains, Statistical Science, Volume 7, Issue 4, 457-472.            ##
    ##   Brooks, S.P. and A. Gelman, (1998) General methods for monitoring convergence    ##
    ##       of iterative simulations, Journal of Computational and Graphical Statistics, ##
    ##       Volume 7, 434-455.                                                           ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Dec. 2016                                            ##
    ## DREAM-Suite                                                                      ##
    ## ################################################################################## ##

    n, d, N = chain.shape
    
    # Early exit if there are fewer than 10 iterations
    if n < 10:
        return np.nan * np.ones(d), np.nan

    # STEP 0: Compute the chain means and store in a N x d matrix
    mu_chains = np.mean(chain, axis = 0).T  # N x d
    # STEP 1: Compute the N within-chain variances
    s2_chains = np.array([np.var(chain[:, :, i], axis = 0, ddof = 1) for i in range(N)])  # N x d
    # STEP 2: Compute the N within-chain covariances
    cov_chains = np.array([np.cov(chain[:, :, i], rowvar=False) for i in range(N)])  # N x d x d
    
    # Univariate hatR diagnostics
    # STEP 1: Compute variance B of N chain means
    B = n * np.var(mu_chains, axis = 0)  # d
    # STEP 2: Compute 1xd vector W with mean of within-chain variances
    W = np.mean(s2_chains, axis = 0)  # d
    # STEP 3: Estimate target variance = sum of within- and between-chain s2's
    sigma2 = ((n - 1) / n) * W + (1 / n) * B
    # STEP 4: Compute univariate hatR diagnostic for each parameter
    hatR = np.sqrt((N + 1) / N * (sigma2 / W) - (n - 1) / (N * n))

    # Multivariate hatRd diagnostic
    # STEP 1: Compute dxd matrix with mean W of within-chain covariances
    W_cov = np.mean(cov_chains, axis = 0) + np.finfo(float).eps * np.eye(d)
    # STEP 2: Compute covariance B of N chain means
    B_cov = np.cov(mu_chains.T) + np.finfo(float).eps * np.eye(d)
    # Now check the covariance matrix, C, is it singular or not?
    if np.linalg.det(W_cov) == 0:
        with open("warning_file.txt", "a+") as fid:
            warning_message = (f"DREAM-Suite WARNING: Singular covariance matrix detected: W_cov"
                               f"R-statistic of Brooks and Gelman at iteration {t}.\n")
            fid.write(warning_message)
        # Apply Tikhonov regularization
        W_inv = np.linalg.inv(W_cov + 1e-6 * np.eye(d))  
    else:
        W_inv = np.linalg.inv(W_cov)

    # STEP 3: Compute multivariate scale reduction factor, hatRd
    hatRd = np.sqrt((N + 1) / N * np.max(np.abs(np.linalg.eigvals(np.linalg.inv(W_cov) @ B_cov))) + (n - 1) / n)

    return hatR, hatRd


# Latin Hypercube Sampling function
def LH_sampling(mn, mx, N):
    # ####################################################################### #
    # This function performs Latin Hypercube sampling                         #
    # ####################################################################### #

    if len(mn.shape) == 2:
        d = mn.shape[1]                                             # Number of parameters
    else:
        d = len(mn)
    rng = np.array(mx) - np.array(mn)                               # 1 x d vector with parameter ranges
    y =  np.random.rand(N, d)                                       # N x d matrix with uniform random labels
    # really important change below so that X stays in bound! as list is from 0 - N-1 rather than 1 to N
    id_matrix = 1 + np.argsort(np.random.rand(N, d), axis = 0)      # Random sort (1:N without replacement)
    M = (id_matrix - y) / N                                         # Multiplier matrix (y introduces randomness)
    R = np.add(np.multiply(M, rng), mn)                             # N x d matrix of stratified LH samples

    return R


def get_label(counter):
    # ####################################################################### #
    # This function turns an integer into a letter/letters                    #
    # ####################################################################### #

    label = ""
    while counter >= 0:
        label = chr(65 + (counter % 26)) + label
        counter = counter // 26 - 1
    return label


def convert_memoryview_to_array(plugin):
    # ####################################################################### #
    # Convert memoryview objects in plugin to numpy arrays to make picklable  #
    # ####################################################################### #

    for key, value in plugin.items():
        if isinstance(value, memoryview):
            plugin[key] = np.array(value)  # Convert memoryview to ndarray
    
    return plugin


def distribute_tasks(N, CPU):
    # ####################################################################### #
    # Split the task ranges evenly across workers                             #
    # ####################################################################### #

    chunk_size = N // CPU
    task_ranges = []
    for i in range(CPU):
        start_idx = i * chunk_size
        # Ensure the last worker gets the remaining tasks
        end_idx = N if i == CPU - 1 else (i + 1) * chunk_size
        task_ranges.append((start_idx, end_idx))
    return task_ranges


def X_unnormalize(Xn, Par_info):
    # ####################################################################### #
    # This function backtransforms to original space normalzed parameter val. #
    # ####################################################################### #
    
    if Par_info['norm'] == 0:
        # No normalization
        Xun = Xn
    elif Par_info['norm'] == 1:
        # Safest implementation: applies transformation
        Xun = Xn * (Par_info['maxun'] - Par_info['minun']) + Par_info['minun']
    
    return Xun


def X_normalize(Xun, Par_info):
    # ####################################################################### #
    # This function transforms to normalized parameter space sampled values   #
    # ####################################################################### #
    
    if Par_info['norm'] == 0:
        # No normalization
        Xn = Xun
    elif Par_info['norm'] == 1:
        # Safest implementation: applies normalization
        Xn = (Xun - Par_info['minun']) / (Par_info['maxun'] - Par_info['minun'])
    
    return Xn


###############################
def check_prior(Par_info, DREAMPar, M=100):
    # ####################################################################### #
    # This function Determines whether the prior handle returns pdf or logpdf #
    # ####################################################################### #

    r_arg = 'pdf'  # Default is 'pdf'
    
    # Draw samples from prior and then evaluate prior handle of pdf (or logpdf)
    # This may not work for highly peaked priors as pdf > 1 too often
    if Par_info['u'] == 'yes':  # Univariate case
        PrX = np.full((M, DREAMPar['d']), np.nan)  # Initialize array to store prior samples
        for j in range(DREAMPar['d']):
            fnc = Par_info['prior'][j]
            for zz in range(M):
                # Sample from prior
                z = fnc.rvs()
                # Evaluate the prior for each dimension and sample
                PrX[zz,j] = fnc.pdf(z) 

        # Check if any of the samples are negative
        if np.sum(PrX < 0) > 0:
            r_arg = 'logpdf'

    else:  # Multivariate case
        PrX = np.full((M, 1), np.nan)  # Initialize array for multivariate case
        fnc = Par_info['prior'][0]
        for zz in range(M):
            # Sample from prior
            z = fnc.rvs()
            # Evaluate the prior for each dimension and sample
            PrX[zz,0] = fnc.pdf(z) 
        # Check if any of the samples are negative
        if np.sum(PrX < 0) > 0:
            r_arg = 'logpdf'

    return r_arg


def genparset(chain):
    # ####################################################################### #
    # This function generates a matrix P from sampled chain trajectories      #
    # ####################################################################### #

    T, d, N = chain.shape  # Get dimensions: #samples, #parameters, #chains

    if T == 0:
        return np.array([])  # If no samples, return an empty array
    else:
        id_ = np.arange(1, T + 1)  # ID for each chain sample (1, 2, ..., T)
        P = np.full((N * T, d + 1), np.nan)  # Initialize matrix P with NaNs

        for z in range(N):  # For each chain
            # Copy each chain's data to P and add sample IDs
            P[z * T:(z + 1) * T, 0:d] = chain[:, :, z]  # Parameters for chain z
            P[z * T:(z + 1) * T, d] = id_  # Sample IDs for chain z
        
        # Sort P based on the last column (the sample ID), and remove the ID column
        P_sorted = P[np.argsort(P[:, d]), :]  # Sort by the last column (sample ID)
        P_sorted = P_sorted[:, 0:d]  # Remove the ID column

        return P_sorted


def calcnbins(x, method='middle', minb=1, maxb=np.inf):
    # ####################################################################### #
    # Computes the "ideal" number of bins for a histogram                     #
    # ####################################################################### #
    
    # Input checking
    if not isinstance(x, (np.ndarray, list, np.generic)):
        raise ValueError('The x argument must be numeric or logical.')

    x = np.asarray(x)
    
    # Ensure the array is real, discard imaginary part
    if np.iscomplexobj(x):
        x = np.real(x)
        print('Warning: Imaginary parts of x will be ignored.')
    
    # Ensure x is a vector (1D array)
    if x.ndim != 1:
        x = x.flatten()
        print('Warning: x will be coerced to a vector.')
    
    # Remove NaN values
    x = x[~np.isnan(x)]
    if len(x) == 0:
        raise ValueError("x must contain at least one valid number.")
    
    # Choose method if not specified
    valid_methods = ['fd', 'scott', 'sturges', 'all', 'middle']
    if method not in valid_methods:
        raise ValueError(f"Unknown method: {method}")
    
    # Method selection
    if method == 'fd':
        nbins = calc_fd(x)
    elif method == 'scott':
        nbins = calc_scott(x)
    elif method == 'sturges':
        nbins = calc_sturges(x)
    elif method == 'middle':
        nbins = [calc_fd(x), calc_scott(x), calc_sturges(x)]
        nbins = np.median(nbins)
    elif method == 'all':
        nbins = {
            'fd': calc_fd(x),
            'scott': calc_scott(x),
            'sturges': calc_sturges(x)
        }
    
    # Confine number of bins to the acceptable range
    nbins = confine_to_range(nbins, minb, maxb)
    
    return nbins


def calc_fd(x):
    # ####################################################################### #
    # Freedman-Diaconis rule for number of bins for given array               #
    # ####################################################################### #

    h = np.subtract(*np.percentile(x, [75, 25]))  # Interquartile range (IQR)
    if h == 0:
        h = 2 * np.median(np.abs(x - np.median(x)))  # Median absolute deviation (MAD)
    
    if h > 0:
        nbins = np.ceil((np.max(x) - np.min(x)) / (2 * h * len(x) ** (-1/3)))
    else:
        nbins = 1
    return nbins


def calc_scott(x):
    # ####################################################################### #
    # Scott's method for number of bins for given array                       #
    # ####################################################################### #

    h = 3.5 * np.std(x) * len(x) ** (-1/3)
    if h > 0:
        nbins = np.ceil((np.max(x) - np.min(x)) / h)
    else:
        nbins = 1
    return nbins


def calc_sturges(x):
    # ####################################################################### #
    # Sturges' method for number of bins for given array                      #
    # ####################################################################### #

    nbins = np.ceil(np.log2(len(x)) + 1)
    return nbins


def confine_to_range(x, lower, upper):
    # ####################################################################### #
    # Function ensures that bin count is within specified range               #
    # ####################################################################### #

    x = np.maximum(x, lower)
    x = np.minimum(x, upper)
    return np.floor(x)


def safe_int(value):
    return int(value) if value else 0  # or return a default value if empty


def safe_float(value):
    return float(value) if value else 0.0  # or return a default value if empty


############################### Interface with user-defined function ############################### 

def get_function_handle(Func_name):
    # ####################################################################### #
    # This function returns a function handle (reference to the function) for #
    # the specified Func_name. This avoids repeated imports and lookups       #
    # ####################################################################### #

    # Check whether Func_name is in duplicate form ("a.a" rather than "a" alone)
    Func_name = create_duplicate_with_dot(Func_name)

    # Check if Func_name contains a dot (.) to separate module and function   
    if '.' in Func_name:
        module_name, function_name = Func_name.rsplit('.', 1)
    else:
        # If no dot is present, assume Func_name is the function name in the current module
        module_name = __name__      # Use the current module
        function_name = Func_name
    
    try:
        # Dynamically import the module if it's not already imported
        module = importlib.import_module(module_name)
    except ModuleNotFoundError:
        raise ValueError(f"Module '{module_name}' not found.")
    
    # Get the function from the module
    Func = getattr(module, function_name, None)
    
    if Func is None:
        raise ValueError(f"Function '{function_name}' not found in module '{module_name}'.")
    
    # Ensure the function is callable
    if not callable(Func):
        raise ValueError(f"The object '{function_name}' in '{module_name}' is not callable.")
    
    # Return the function handle (reference to the function)
    return Func


# Function to recursively convert memoryview objects to bytearrays
def convert_memoryview(obj):
    if isinstance(obj, memoryview):
        # Convert memoryview to bytearray
        return bytearray(obj)
    elif isinstance(obj, dict):
        # Recursively convert memoryview objects inside a dictionary
        return {key: convert_memoryview(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        # Recursively convert memoryview objects inside a list
        return [convert_memoryview(item) for item in obj]
    elif isinstance(obj, tuple):
        # Recursively convert memoryview objects inside a tuple
        return tuple(convert_memoryview(item) for item in obj)
    elif isinstance(obj, np.ndarray):
        # Convert memoryview from NumPy arrays to ndarray (if necessary)
        return np.array(obj)
    else:
        # Return other objects as they are (assumed to be picklable)
        return obj


def copy_model_files(source_dir, target_dir):
    # ####################################################################### #
    # Copy model files from source directory to the target worker directory   #
    # ####################################################################### #    

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    
    # List all files in the source directory and copy them to the target directory
    for filename in os.listdir(source_dir):
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)
        
        if os.path.isfile(source_file):
            shutil.copy2(source_file, target_file)  # Copy file with metadata
        elif os.path.isdir(source_file):
            shutil.copytree(source_file, target_file)  # Copy directory recursively


def worker_task(worker_id, start_idx, end_idx, X, func_handle, plugin, base_dir):
    
    if base_dir is not None:
        # Create a worker-specific directory
        worker_dir = os.path.join(base_dir, f"worker_{worker_id}")   
        # Change to the worker's specific directory
        os.chdir(worker_dir)

    # Execute model
    results = []
    for idx in range(start_idx, end_idx):
        if plugin is not None:
            result = func_handle(X[idx, :], plugin)
        else:
            result = func_handle(X[idx, :])

        results.append(result)

    # Return the results for this worker
    return results


def cleanup_worker_directories(base_dir, N):
    # ####################################################################### #
    # Clean up (delete) worker directories after all generations are complete #
    # ####################################################################### #

    for worker_id in range(N):
        worker_dir = os.path.join(base_dir, f'worker_{worker_id}')
        # Delete the worker directory
        if os.path.exists(worker_dir):
            shutil.rmtree(worker_dir)  


############################# End Interface with user-defined function ############################# 

##################################### Miscellaneous functions ######################################

def UL(iflag, nuisvar, method, y_sim, y_meas, sigma=None, N=1):
    ## ################################################################################## ##
    ## Universal likelihood (UL) function: correlated, heteroscedastic and non-Gaussian   ##
    ## errors                                                                             ##
    ##  1. Partial stand. residuals expected to follow skew generalized Student's t (SGT) ##
    ##     distribution                                                                   ##
    ##  2. Serial correlation described using an autoregressive model - up to order 2     ##
    ## UL function performs standardization before treatment serial correlation           ##
    ##                                                                                    ##
    ##  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r] = GL_plus(iflag,nuisvar,method,...     ##
    ##                 y_sim,Y_meas,sigma,N)                                              ##
    ##  where                                                                             ##
    ##   iflag      [input] Estimation ('est') or simulation ('sim')                      ##
    ##   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         ##
    ##    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       ##
    ##    s1:  nuisvar(2) Slope of linear heteroscedastic model                           ##
    ##    lbd: nuisvar(3) Skewness (-1,1): (-1,0) neg. skew, 0 symmetric, (0,1) pos. skew ##
    ##    p:   nuisvar(4) Kurtosis (> 0) The larger p is, the more uniform distribution   ##
    ##    q:   nuisvar(5) Kurtosis (> 0) The larger q is, the peakier the distribution    ##
    ##    fi1: nuisvar(5) First-order AR coefficient (0,1: check book)                    ##
    ##    fi2: nuisvar(7) Second-order AR coefficient (0,1: check book)                   ##
    ##   method     [input] Treatment of s0 and/or s1: [1-8]                              ##
    ##   y_sim      [input] n x 1 vector with simulated values                            ##
    ##   Y_meas     [input] n x 1 vector with observed values                             ##
    ##   sigma      [input] Optional: Measurement sigma defined by user                   ##
    ##   N          [input] Optional: Number of replicates - resampling                   ##
    ##   loglik     [outpt] n x 1 vector of log-likelihood values                         ##
    ##   std_e      [outpt] n x 1 vector with standard deviation of raw residuals         ##
    ##   eps_n      [outpt] n x 1 vector with standardized decorrelated residuals         ##
    ##   f_eps_n    [outpt] n x 1 vector with density stand. decorrelated residuals       ##
    ##   Y_r        [outpt] n x N matrix of replicate simulations (for Bayes_pdf)         ##
    ##                                                                                    ##
    ##  Reference:                                                                        ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             ##
    ##       https://www.sciencedirect.com/science/article/pii/S002216942201112X          ##
    ##                                                                                    ##
    ##  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  ##
    ##         2. Measurement error computed from simulated data                          ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Dec. 2016                                            ##
    ## DREAM-Suite                                                                      ##
    ## ################################################################################## ##

    # Unpack nuisance variables
    s0, s1, lbd, p, q, fi1, fi2 = nuisvar

    # Number of samples
    n = len(y_meas) 

    # Initialize output variables
    loglik = -np.inf * np.ones((n,1))
    std_e = eps_n = f_eps_n = np.zeros((n,1))
    Y_r = None

    # Compute kappa and mu for the SGT distribution
    if p * q < 2:
        return loglik, std_e, eps_n, f_eps_n, Y_r

    kappa_lpq = sp.beta(1/p, q/p) / np.sqrt((1 + 3 * lbd**2) * sp.beta(1/p, q/p) * sp.beta(3/p, (q-2)/p) - 4 * lbd**2 * sp.beta(2/p, (q-1)/p)**2)
    mu_lpq = 2 * kappa_lpq * lbd * sp.beta(2/p, (q-1)/p) / sp.beta(1/p, q/p)

    e = y_meas - y_sim  # Raw residuals

    # Standardization based on the method parameter
    if method == 0:
        std_e = sigma  # User-specified
    elif method == 1:
        std_e = s0 * np.ones(n)  # Constant s0
    elif method == 2:
        std_e = np.std(e) * np.ones(n)  # Constant based on raw residuals
    elif method == 3:
        std_e = s0 * np.ones(n)  # Constant s0, no s1
    elif method == 4:
        std_e = s0 * np.ones(n)  # Constant s0, no s1
    elif method == 5:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 6:
        std_e = s0 + s1 * y_sim  # Heteroscedastic model with s1 fixed
    elif method == 7:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 8:
        std_e = s0 + s1 * y_sim  # Heteroscedastic model with s1 estimated

    # Compute theoretical std. of partial residuals
    std_eps = np.sqrt((1 + fi1**2 - fi2**2 - 2 * fi1**2 / (1 - fi2)))

    # Standardized residuals
    fi_p = [1, -fi1, -fi2]                              # Coefficients of AR filter
    e_n = e / std_e                                     # Studentized raw residuals (unit variance)
    # eps needs to be a vector with 1 dimension!
    # then lfilter is equal to filter in MATLAB
    eps = lfilter(fi_p, 1, e_n)                         # Partial residuals
    eps_n = eps / std_eps                               # Standardized partial residuals

    # Note: MATLAB returns a complex number if fi1 + fi2 too large
    #       Python returns nan instead
    # Likelihood estimation: For model training
    if iflag == 'est':      
        if not np.isreal(std_eps) or np.isnan(std_eps):  # Check for non-real values    
            return loglik, std_e, eps_n, f_eps_n, Y_r

        eps_sgt = eps_n + mu_lpq                        # Shift for SGT distribution
        loglik = -np.log(std_eps) - np.log(std_e) + np.log(p) - np.log(2) - np.log(kappa_lpq) - sp.betaln(1/p, q/p) - \
                 ((q + 1) / p) * np.log(1 + np.abs(eps_sgt / (kappa_lpq * (1 + lbd * np.sign(eps_sgt)))) ** p)
        # Next sentence is not required - only for postprocessing
        f_eps_n = f_SGT(eps_n, lbd, p, q)               # SGT density of eps_n
    
    # Simulation: Replicates for uncertainty quantification
    elif iflag == 'sim':
        eps_n = SGTrnd(0, 1, lbd, p, q, n, N)           # Draw nxm matrix of standardized partial residuals from SGT distribution
        dl_eps = np.zeros((n, N))                       # nxm matrix of zero-mean correction: Scharnagl 2015: dl_eps = fi1 * mean_of_e_n
        eps = std_eps * eps_n + dl_eps                  # nxm matrix of non-standardized partial residuals
        e_n = lfilter([1], fi_p, eps, axis = 0)         # nxN matrix of studentized raw residuals, var(e_n) = ones(1, N)
        E_r = std_e[:, None] * e_n                      # nxm matrix of non-standardized raw residuals
        Y_r = y_sim[:, None] + E_r                      # nxm matrix of replicates of simulation

    return loglik, std_e, eps_n, f_eps_n, Y_r


def GL_plus(iflag, nuisvar, method, y_sim, y_meas, sigma=None, N=1):
    ## ################################################################################## ##
    ## Modified Generalized likelihood function: correlated, heteroscedastic and          ##
    ## non-Gaussian errors                                                                ##
    ##  1. Partial stand. residuals expected to follow skew exponential power (SEP)       ##
    ##     distribution                                                                   ##
    ##  2. Serial correlation described using an autoregressive model - up to order 2     ##
    ## GL_plus function performs standardization before treatment serial correlation      ##
    ##                                                                                    ##
    ##  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r] = GL_plus(iflag,nuisvar,method,...     ##
    ##                 y_sim,y_meas,sigma,N)                                              ##
    ##  where                                                                             ##
    ##   iflag      [input] Estimation ('est') or simulation ('sim')                      ##
    ##   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         ##
    ##    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       ##
    ##    s1:  nuisvar(2) Slope of linear heteroscedastic model                           ##
    ##    ba:  nuisvar(3) Kurtosis (-1: uniform, 0: normal; 1: Laplace)                   ##
    ##    xi:  nuisvar(4) Skewness (1: symmetric; <1: negative skew; >1: positive skew)   ##
    ##    fi1: nuisvar(5) First-order AR coefficient (0,1)                                ##
    ##    fi2: nuisvar(6) Second-order AR coefficient (0,1)                               ##
    ##   method     [input] Treatment of s0 and/or s1: [1-8]                              ##
    ##   y_sim      [input] n x 1 vector with simulated values                            ##
    ##   y_meas     [input] n x 1 vector with observed values                             ##
    ##   sigma      [input] Optional: Measurement sigma defined by user                   ##
    ##   N          [input] Optional: Number of replicates - resampling                   ##
    ##   loglik     [outpt] n x 1 vector of log-likelihood values                         ##
    ##   std_e      [outpt] n x 1 vector with standard deviation of raw residuals         ##
    ##   eps_n      [outpt] n x 1 vector with standardized decorrelated residuals         ##
    ##   f_eps_n    [outpt] n x 1 vector with density stand. decorrelated residuals       ##
    ##   Y_r        [outpt] n x N matrix of replicate simulations (for Bayes_pdf)         ##
    ##                                                                                    ##
    ##  Reference:                                                                        ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             ##
    ##       https://www.sciencedirect.com/science/article/pii/S002216942201112X          ##
    ##   Schoups, G., and J. A. Vrugt (2010), A formal likelihood function for parameter  ##
    ##       and predictive inference of hydrologic models with correlated,               ##
    ##       heteroscedastic, and non-Gaussian errors, Water Resources Research, 46,      ##
    ##       W10531, doi:10.1029/2009WR008933                                             ##
    ##                                                                                    ##
    ##  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  ##
    ##         2. Measurement error computed from simulated data                          ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Dec. 2014                                            ##
    ## DREAM-Suite                                                                      ##
    ## ################################################################################## ##
    
    # Unpack nuisance variables
    s0, s1, ba, xi, fi1, fi2 = nuisvar
    
    # Number of samples
    n = len(y_meas)              

    # Initialize output variables
    loglik = -np.inf * np.ones((n,1))
    std_e = eps_n = f_eps_n = np.zeros((n,1))
    Y_r = None

    fi_p = [1, -fi1, -fi2]      # Coefficients of AR filter
    e = y_meas - y_sim

    # Standardization based on the method parameter
    if method == 0:
        std_e = sigma  # User-specified
    elif method == 1:
        std_e = s0 * np.ones(n)  # Constant s0
    elif method == 2:
        std_e = np.std(e) * np.ones(n)  # Constant based on raw residuals
    elif method == 3:
        std_e = s0 * np.ones(n)  # Constant s0, no s1
    elif method == 4:
        std_e = s0 * np.ones(n)  # Constant s0, no s1
    elif method == 5:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 6:
        std_e = s0 + s1 * y_sim  # Heteroscedastic model with s1 fixed
    elif method == 7:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 8:
        std_e = s0 + s1 * y_sim  # Heteroscedastic model with s1 estimated
    
    # Compute theoretical std. of partial residuals
    std_eps = np.sqrt((1 + fi1**2 - fi2**2 - 2 * fi1**2 / (1 - fi2)) * 1)
    # Note: MATLAB returns a complex number if fi1 + fi2 too large
    #       Python returns nan instead
    # Likelihood estimation: For model training
    if iflag == 'est':      
        if not np.isreal(std_eps) or np.isnan(std_eps):     # Check for non-real values    
            return loglik, std_e, eps_n, f_eps_n, Y_r

        e_n = e / std_e                                     # Studentized raw residuals
        # eps needs to be a vector with 1 dimension!
        # then lfilter is equal to filter in MATLAB
        eps = lfilter(fi_p, 1, e_n)                         # Apply AR filter (AR(1)/AR(2))
        eps_n = eps / std_eps                               # Standardize partial residuals
        f_eps_n, logf_eps_n = f_SEP(eps_n, ba, xi)          # Compute density
        loglik = logf_eps_n - np.log(std_e) - np.log(std_eps)
    # Simulation: Replicates for uncertainty quantification
    elif iflag == 'sim':    
        eps_n = SEPrnd(ba, xi, n, N)                        # Draw nxm matrix of standardized partial residuals from SEP distribution
        dl_eps = np.zeros((n, N))                           # nxm matrix of zero-mean correction: Scharnagl 2015: dl_eps = fi1 * mean_of_e_n
        eps = std_eps * eps_n + dl_eps                      # nxm matrix of non-standardized partial residuals
        e_n = lfilter([1], fi_p, eps, axis = 0)             # nxN matrix of studentized raw residuals, var(e_n) = ones(1, N)
        E_r = std_e[:, None] * e_n                          # nxm matrix of non-standardized raw residuals
        Y_r = y_sim[:, None] + E_r                          # nxm matrix of replicates of simulation

    return loglik, std_e, eps_n, f_eps_n, Y_r


def SL(iflag, nuisvar, method, y_sim, y_meas, sigma, N=1):
    ## ################################################################################## ##
    ## Student Likelihood (SL) function: correlated, non-constant and non-Gaussian errors ##
    ##  1. Partial stand. residuals expected to follow skew Student t (SST) distribution  ##
    ##  2. Serial correlation treated using an AR(1) or AR(2) model                       ##
    ##                                                                                    ##
    ## This SST function performs standardization before treatment correlation            ##
    ##                                                                                    ##
    ##  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r] = SL(iflag,nuisvar,method,y_sim,...    ##
    ##                 y_meas,sigma,N)                                                    ##
    ##  where                                                                             ##
    ##   iflag      [input] Estimation ('est') or simulation ('sim')                      ##
    ##   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         ##
    ##    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       ##
    ##    s1:  nuisvar(2) Slope of linear heteroscedastic model                           ##
    ##    nu:  nuisvar(3) Degrees of freedom (q > 0)                                      ##
    ##    xi:  nuisvar(4) Skewness parameter (xi > 0)                                     ##
    ##    fi1: nuisvar(5) First-order AR coefficient (0,1: check book)                    ##
    ##    fi2: nuisvar(6) Second-order AR coefficient (0,1: check book)                   ##
    ##   method     [input] Treatment of s0 and/or s1: [1-8]                              ##
    ##   y_sim      [input] n x 1 vector with simulated values                            ##
    ##   y_meas     [input] n x 1 vector with observed values                             ##
    ##   sigma      [input] Optional: Measurement sigma defined by user                   ##
    ##   N          [input] Optional: Number of replicates - resampling                   ##
    ##   loglik     [outpt] n x 1 vector of log-likelihood values                         ##
    ##   std_e      [outpt] n x 1 vector with standard deviation of raw residuals         ##
    ##   eps_n      [outpt] n x 1 vector with standardized decorrelated residuals         ##
    ##   f_eps_n    [outpt] n x 1 vector with density stand. decorrelated residuals       ##
    ##   Y_r        [outpt] n x N matrix of replicate simulations (for Bayes_pdf)         ##
    ##                                                                                    ##
    ##  Reference:                                                                        ##
    ##   Scharnagl, B., S.C. Iden, W. Durner, H. Vereecken, and M. Herbst (2015), Inverse ##
    ##       modelling of in situ soil water dynamics: accounting for heteroscedastic,    ##
    ##       autocorrelated, and non-Gaussian distributed residuals, Hydrology and Earth  ##
    ##       System Sciences Discussions, vol. 12, pp. 2155-2199, 2015.                   ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             ##
    ##       https://www.sciencedirect.com/science/article/pii/S002216942201112X          ##
    ##                                                                                    ##
    ##  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  ##
    ##         2. Measurement error computed from simulated data                          ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Dec. 2016                                            ##
    ## DREAM-Suite                                                                      ##
    ## ################################################################################## ##

    # Unpack nuisance variables
    s0, s1, nu, xi, fi1, fi2 = nuisvar

    # Number of samples
    n = len(y_meas)    

    # Initialize output variables
    loglik = -np.inf * np.ones((n,1))
    std_e = eps_n = f_eps_n = np.zeros((n,1))
    Y_r = None
    
    # Compute first moment of SST distribution
    def M(j, nu):
        return sp.gamma((j + 1) / 2) * sp.gamma((nu - j) / 2) * (nu - 2) ** (j / 2) / (np.sqrt(np.pi) * sp.gamma(nu / 2))

    # Compute second moment of SST distribution
    def logM(j, nu):
        return sp.gammaln((j + 1) / 2) + sp.gammaln((nu - j) / 2) + (j / 2) * np.log(nu - 2) - 0.5 * np.log(np.pi) - sp.gammaln(nu / 2)
    
    M1 = M(1, nu) if method == 0 else np.exp(logM(1, nu))
    M2 = M(2, nu) if method == 0 else np.exp(logM(2, nu))

    # Now compute mu_xi and sig_xi
    mu_xi = M1 * (xi - 1 / xi)
    sig_xi = np.sqrt((M2 - M1 ** 2) * (xi ** 2 + xi ** -2) + 2 * M1 ** 2 - M2)

    fi_p = [1, -fi1, -fi2]      # Coefficients of AR filter
    e = y_meas - y_sim          # Raw residuals

    # Standardization based on the method parameter
    if method == 0:
        std_e = sigma  # User-specified
    elif method == 1:
        std_e = s0 * np.ones(n)  # Constant s0
    elif method == 2:
        std_e = np.std(e) * np.ones(n)  # Constant based on raw residuals
    elif method == 3:
        std_e = s0 * np.ones(n)  # Constant s0, no s1
    elif method == 4:
        std_e = s0 * np.ones(n)  # Constant s0, no s1
    elif method == 5:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 6:
        std_e = s0 + s1 * y_sim  # Heteroscedastic model with s1 fixed
    elif method == 7:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 8:
        std_e = s0 + s1 * y_sim  # Heteroscedastic model with s1 estimated
    
    # Compute theoretical std of partial residuals
    std_eps = np.sqrt(1 + fi1 ** 2 - fi2 ** 2 - 2 * fi1 ** 2 / (1 - fi2))
    # Note: MATLAB returns a complex number if fi1 + fi2 too large
    #       Python returns nan instead
    # Likelihood estimation: For model training
    if iflag == 'est':      
        if not np.isreal(std_eps) or np.isnan(std_eps):  # Check for non-real values    
            return loglik, std_e, eps_n, f_eps_n, Y_r
        
        e_n = e / std_e                                     # Studentized raw residuals; var(e_n) = 1
        # eps needs to be a vector with 1 dimension!
        # then lfilter is equal to filter in MATLAB
        eps = lfilter(fi_p, 1, e_n)                         # Partial residuals
        eps_n = eps / std_eps                               # Standardized partial residuals
        a_sst = (mu_xi + sig_xi * eps_n) / (xi ** np.sign(mu_xi + sig_xi * eps_n))

        # Compute log-likelihood
        loglik = -np.log(std_e) - np.log(std_eps) + np.log(2) + np.log(sig_xi) + sp.gammaln((nu + 1) / 2) \
            - np.log(xi + 1 / xi) - sp.gammaln(nu / 2) - 0.5 * np.log(np.pi) - 0.5 * np.log(nu - 2) \
                - (nu + 1) / 2 * np.log(1 + 1 / (nu - 2) * a_sst ** 2)
        # Next sentence is not required - only for postprocessing
        f_eps_n = f_SST(eps_n, nu, xi)                      # SST density of eps_n

    # Simulation: Replicates for uncertainty quantification
    elif iflag == 'sim':
        eps_n = np.random.standard_t(nu, size=(n, N))       # Draw nxm matrix of standardized partial residuals from SST(0,1,nu,xi)
        dl_eps = np.zeros((n, N))                           # nxm matrix of zero-mean correction: Scharnagl 2015: dl_eps = fi1 * mean_of_e_n
        eps = std_eps * eps_n + dl_eps                      # nxm matrix of non-standardized partial residuals
        e_n = lfilter([1], fi_p, eps, axis = 0)             # nxN matrix of studentized raw residuals, var(e_n) = ones(1, N)
        E_r = std_e[:, None] * e_n                          # nxm matrix of non-standardized raw residuals
        Y_r = y_sim[:, None] + E_r                          # nxm matrix of replicates of simulation

    return loglik, std_e, eps_n, f_eps_n, Y_r


def LAPL(iflag, nuisvar, method, y_sim, y_meas, sigma=None, N=1):
    ## ################################################################################## ##
    ## Laplacian Likelihood (LAPL) function: correlated and non-constant errors           ##
    ##  1. Partial stand. residuals expected to follow Laplace distribution distribution  ##
    ##  2. Serial correlation treated using an AR(1) or AR(2) model                       ##
    ##                                                                                    ##
    ## LAPL function performs residual standardization before treatment serial correlatn  ##
    ##                                                                                    ##
    ##  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r] = LAPL(iflag,nuisvar,method,y_sim,...  ##
    ##                 y_meas,sigma,N)                                                    ##
    ##  where                                                                             ##
    ##   iflag      [input] Estimation ('est') or simulation ('sim')                      ##
    ##   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         ##
    ##    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       ##
    ##    s1:  nuisvar(2) Slope of linear heteroscedastic model                           ##
    ##    fi1: nuisvar(3) First-order AR coefficient (0,1)                                ##
    ##   method     [input] Treatment of s0 and/or s1: [1-8]                              ##
    ##   y_sim      [input] n x 1 vector with simulated values                            ##
    ##   y_meas     [input] n x 1 vector with observed values                             ##
    ##   sigma      [input] Optional: Measurement sigma defined by user                   ##
    ##   N          [input] Optional: Number of replicates - resampling                   ##
    ##   loglik     [outpt] n x 1 vector of log-likelihood values                         ##
    ##   std_e      [outpt] n x 1 vector with standard deviation of raw residuals         ##
    ##   eps_n      [outpt] n x 1 vector with standardized decorrelated residuals         ##
    ##   f_eps_n    [outpt] n x 1 vector with density stand. decorrelated residuals       ##
    ##   Y_r        [outpt] n x N matrix of replicate simulations (for Bayes_pdf)         ##
    ##                                                                                    ##
    ##  Reference:                                                                        ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             ##
    ##       https://www.sciencedirect.com/science/article/pii/S002216942201112X          ##
    ##                                                                                    ##
    ##  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  ##
    ##         2. Measurement error computed from simulated data                          ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Dec. 2014                                            ##
    ## DREAM-Suite                                                                      ##
    ## ################################################################################## ##

    lik_calc = 'cond'       # Conditional normal log-likelihood
                            # Other option: exact

    # Unpack nuisance variables
    s0, s1, fi1 = nuisvar

    fi2 = 0                 # Not used in the original function   
    n = len(y_meas)         # Number of observations

    # Initialize output variables
    loglik = -np.inf * np.ones((n,1))
    std_e = eps_n = f_eps_n = np.zeros((n,1))
    Y_r = None

    fi_p = [1, -fi1, -fi2]      # Coefficients of AR filter
    e = y_meas - y_sim          # Residuals (y_meas - y_sim)

    # Standardization based on the method parameter
    if method == 0:
        std_e = sigma  # User-specified
    elif method == 1:
        std_e = s0 * np.ones(n)  # Constant s0
    elif method == 2:
        std_e = np.std(e) * np.ones(n)  # Constant based on raw residuals
    elif method == 3:
        std_e = s0 * np.ones(n)  # Constant s0, no s1
    elif method == 4:
        std_e = s0 * np.ones(n)  # Constant s0, no s1
    elif method == 5:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 6:
        std_e = s0 + s1 * y_sim  # Heteroscedastic model with s1 fixed
    elif method == 7:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 8:
        std_e = s0 + s1 * y_sim  # Heteroscedastic model with s1 estimated
    
    # Compute theoretical std. of partial residuals
    std_eps = np.sqrt((1 + fi1**2 - fi2**2 - 2 * fi1**2 / (1 - fi2)) * 1)
    # Note: MATLAB returns a complex number if fi1 + fi2 too large
    #       Python returns nan instead
    # Likelihood estimation: For model training
    if iflag == 'est':      
        if not np.isreal(std_eps) or np.isnan(std_eps):  # Check for non-real values    
            return loglik, std_e, eps_n, f_eps_n, Y_r
        
        e_n = e / std_e                 # Studentized raw residuals; var(e_n) = 1
        # eps needs to be a vector with 1 dimension!
        # then lfilter is equal to filter in MATLAB
        eps = lfilter(fi_p, 1, e_n)     # Partial residuals
        eps_n = eps / std_eps           # Standardize partial residuals
        
        def f_ell(x, sig):          ## Log-likelihood function
            return -0.5 * np.log(2) - np.log(sig) - np.sqrt(2) * np.abs(x) / sig
        
        # loglik = np.nan * np.ones(n)
        if lik_calc == 'exact':     ## Exact log-likelihood function AR(1)
            loglik[0] = f_ell(e_n[0], np.sqrt(std_eps**2 / (1 - fi1**2)))
            loglik[1:] = np.sum(f_ell(eps[1:], std_eps)) - np.sum(np.log(std_e))
        elif lik_calc == 'cond':    ## Conditional log-likelihood
            loglik = f_ell(eps[:n], std_eps) - np.log(std_e)
        # MATLAB: only if nargout  > 3 --> replaced by N > 1 - is this OK?
        # Next two sentences are not required - only for postprocessing
        lappdf = lambda x, m, s: 1 / (np.sqrt(2) * s) * np.exp(-np.sqrt(2) * np.abs(x - m) / s)
        f_eps_n = lappdf(eps_n, 0, 1)   # Laplace density of eps_n
    
    # Simulation: Replicates for uncertainty quantification
    elif iflag == 'sim':                
        eps_n = LAPrnd(0, 1, n, N)                  # Draw nxm matrix of standardized partial residuals from Laplace distribution
        dl_eps = np.zeros((n, N))                   # nxm matrix of zero-mean correction: Scharnagl 2015: dl_eps = fi1 * mean_of_e_n
        eps = std_eps * eps_n + dl_eps              # nxm matrix of non-standardized partial residuals
        e_n = lfilter([1], fi_p, eps, axis = 0)     # nxN matrix of studentized raw residuals, var(e_n) = ones(1, N)
        E_r = std_e[:, None] * e_n                  # nxm matrix of non-standardized raw residuals
        Y_r = y_sim[:, None] + E_r                  # nxm matrix of replicates of simulation
        
    return loglik, std_e, eps_n, f_eps_n, Y_r


def NL(iflag, nuisvar, method, y_sim, y_meas, sigma=None, N=1):
    ## ################################################################################## ##
    ## Normal Likelihood (NL) function: correlated and non-constant errors                ##
    ##  1. Partial stand. residuals expected to follow normal distribution                ##
    ##  2. Serial correlation treated using an AR(1) or AR(2) model                       ##
    ##                                                                                    ##
    ## NL function performs standardization before treatment correlation                  ##
    ##                                                                                    ##
    ##  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r] = NL(iflag,nuisvar,method,y_sim,...    ##
    ##                 y_meas,sigma,N)                                                    ##
    ##  where                                                                             ##
    ##   iflag      [input] Estimation ('est') or simulation ('sim')                      ##
    ##   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         ##
    ##    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       ##
    ##    s1:  nuisvar(2) Slope of linear heteroscedastic model                           ##
    ##    fi1: nuisvar(3) First-order AR coefficient (0,1)                                ##
    ##    fi2: nuisvar(4) Second-order AR coefficient (0,1)                               ##
    ##   method     [input] Treatment of s0 and/or s1: [1-8]                              ##
    ##   y_sim      [input] n x 1 vector with simulated values                            ##
    ##   y_meas     [input] n x 1 vector with observed values                             ##
    ##   sigma      [input] Optional: Measurement sigma defined by user                   ##
    ##   N          [input] Optional: Number of replicates - resampling                   ##
    ##   loglik     [outpt] n x 1 vector of log-likelihood values                         ##
    ##   std_e      [outpt] n x 1 vector with standard deviation of raw residuals         ##
    ##   eps_n      [outpt] n x 1 vector with standardized decorrelated residuals         ##
    ##   f_eps_n    [outpt] n x 1 vector with density stand. decorrelated residuals       ##
    ##   Y_r        [outpt] n x N matrix of replicate simulations (for Bayes_pdf)         ##
    ##                                                                                    ##
    ##  Reference:                                                                        ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             ##
    ##       https://www.sciencedirect.com/science/article/pii/S002216942201112X          ##
    ##                                                                                    ##
    ##  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  ##
    ##         2. Measurement error computed from simulated data                          ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Oct. 2013                                            ##
    ## DREAM-Suite                                                                      ##
    ## ################################################################################## ##

    lik_calc = 'cond'       # Conditional normal log-likelihood
                            # Other option: exact

    # Unpack nuisance variables
    s0, s1, fi1, fi2 = nuisvar

    n = len(y_meas)         # Number of observations

    # Initialize output variables
    loglik = -np.inf * np.ones((n,1))
    std_e = eps_n = f_eps_n = np.zeros((n,1))
    Y_r = None

    fi_p = [1, -fi1, -fi2]  # AR filter coefficients
    e = y_meas - y_sim      # Raw residuals (n x 1 vector)
    
    # Standardization based on the method parameter
    if method == 0:
        std_e = sigma               # User-specified
    elif method == 1:
        std_e = s0 * np.ones(n)     # Constant s0
    elif method == 2:
        std_e = np.std(e) * np.ones(n)  # Constant based on raw residuals
    elif method == 3:
        std_e = s0 * np.ones(n)     # Constant s0, no s1
    elif method == 4:
        std_e = s0 * np.ones(n)     # Constant s0, no s1
    elif method == 5:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 6:
        std_e = s0 + s1 * y_sim     # Heteroscedastic model with s1 fixed
    elif method == 7:
        std_e, _, exitflag = s1_phantom(s0, e, y_sim)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 8:
        std_e = s0 + s1 * y_sim     # Heteroscedastic model with s1 estimated
    
    # Compute theoretical std of partial residuals
    std_eps = np.sqrt(1 + fi1**2 - fi2**2 - (2 * fi1**2 / (1 - fi2)))

    # Note: MATLAB returns a complex number if fi1 + fi2 too large
    #       Python returns nan instead
    # Likelihood estimation: For model training
    if iflag == 'est':      
        if not np.isreal(std_eps) or np.isnan(std_eps):  # Check for non-real values    
            return loglik, std_e, eps_n, f_eps_n, Y_r
        
        e_n = e / std_e                         # Studentized residuals (variance = 1)
        # eps needs to be a vector with 1 dimension!
        # then lfilter is equal to filter in MATLAB
        eps = lfilter(fi_p, 1, e_n)             # AR(2) filter applied to residuals
        eps_n = eps / std_eps                   # Standardized partial residuals
        
        # Log-likelihood calculation
        if lik_calc == 'exact':     ## Exact log-likelihood function
            loglik = (-n / 2) * np.log(2 * np.pi) - n * np.log(std_eps) \
                     + 0.5 * np.log((1 + fi2)**2 * ((1 - fi2)**2 - fi1**2)) \
                     - (1 + fi2) / (2 * std_eps**2) * (
                         (1 - fi2) * eps[0]**2 - 2 * fi1 * eps[0] * eps[1] + (1 - fi2) * eps[1]**2) \
                     - 0.5 * np.sum(eps_n[2:n]**2) - np.sum(np.log(std_e))
        elif lik_calc == 'cond':    ## Conditional log-likelihood function
            loglik = - (1 / 2) * np.log(2 * np.pi) - np.log(std_eps) \
                     - 0.5 * eps_n[0:n]**2 - np.log(std_e)
        
        # Next sentence is not required - only for postprocessing
        f_eps_n = norm.pdf(eps_n, 0, 1)         # normal density of eps_n
    
    # Simulation: Replicates for uncertainty quantification
    elif iflag == 'sim':            
        eps_n = np.random.normal(0, 1, (n, N))  # Draw nxm matrix of standardized partial residuals from N(0,1)
        dl_eps = np.zeros((n, N))               # nxm matrix of zero-mean correction: Scharnagl 2015: dl_eps = fi1 * mean_of_e_n
        eps = std_eps * eps_n + dl_eps          # Apply correction
        e_n = lfilter([1], fi_p, eps, axis = 0) # nxN matrix of studentized raw residuals, var(e_n) = ones(1, N)
        E_r = std_e[:, None] * e_n              # nxm matrix of non-standardized raw residuals
        Y_r = y_sim[:, None] + E_r              # nxm matrix of replicates of simulation
                
    return loglik, std_e, eps_n, f_eps_n, Y_r


## OBSOLETE GENERALIZED LIKELIHOOD FUNCTION [--> USE GL_PLUS]
def GL(iflag, nuisvar, method, y_sim, y_meas, sigma=None, N=1):
    ## ################################################################################## ##
    ## Generalized likelihood function: correlated, heteroscedastic and                   ##
    ##  non-Gaussian errors                                                               ##
    ##  1. Partial standardized residuals are expected to follow a skew exponential       ##
    ##     power (SEP) distribution                                                       ##
    ##  2. Serial correlation described using an autoregressive model - up to order 2     ##
    ##                                                                                    ##
    ## This GL function treats serial correlation before standardization (= not ideal)    ##
    ## This work has been published in Schoups and Vrugt (2010)                           ##
    ##                                                                                    ##
    ##  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r,y_bc] = GL(iflag,nuisvar,method,...     ##
    ##                 y_sim,y_meas,N)                                                    ##
    ##  where                                                                             ##
    ##   iflag      [input] Estimation ('est') or simulation ('sim')                      ##
    ##   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         ##
    ##    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       ##
    ##    s1:  nuisvar(2) Slope of linear heteroscedastic model                           ##
    ##    ba:  nuisvar(3) Kurtosis (-1: uniform, 0: normal; 1: Laplace)                   ##
    ##    xi:  nuisvar(4) Skewness (1: symmetric; <1: negative skew; >1: positive skew)   ##
    ##    mu1: nuisvar(5) Bias correction parameter                                       ##
    ##    fi1: nuisvar(6) First-order AR coefficient (0,1)                                ##
    ##    fi2: nuisvar(7) Second-order AR coefficient (0,1)                               ##
    ##    fi3: nuisvar(8) Third-order AR coefficient (0,1)                                ##
    ##    fi4: nuisvar(9) Fourth-order AR coefficient (0,1)                               ##
    ##    K:   nuisvar(10) Box-Cox transformation parameter (skewness)                    ##
    ##    lba: nuisvar(11) Box-Cox transformation parameter (heteroscedasticity)          ##
    ##   method     [input] Treatment of s0 and/or s1: [1-8]                              ##
    ##   y_sim      [input] n x 1 vector with simulated values                            ##
    ##   y_meas     [input] n x 1 vector with observed values                             ##
    ##   sigma      [input] Optional: Measurement sigma defined by user                   ##
    ##   N          [input] Optional: Number of replicates - resampling                   ##
    ##   loglik     [outpt] n x 1 vector of log-likelihood values                         ##
    ##   std_e      [outpt] n x 1 vector with standard deviation of raw residuals         ##
    ##   eps_n      [outpt] n x 1 vector with standardized decorrelated residuals         ##
    ##   f_eps_n    [outpt] n x 1 vector with density stand. decorrelated residuals       ##
    ##   Y_r        [outpt] n x N matrix of replicate simulations (for Bayes_pdf)         ##
    ##                                                                                    ##
    ##  Reference:                                                                        ##
    ##  Schoups, G., and J. A. Vrugt (2010), A formal likelihood function for parameter   ##
    ##       and predictive inference of hydrologic models with correlated,               ##
    ##       heteroscedastic, and non-Gaussian errors, Water Resources Research, 46,      ##
    ##       W10531, doi:10.1029/2009WR008933                                             ##
    ##                                                                                    ##
    ##  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  ##
    ##         2. Measurement error computed from simulated data                          ##
    ##         3. This function is OBSOLETE/ERRONEOUS, please use GL_plus                 ##
    ##                                                                                    ##
    ## Written by Gerrit Schoups (TU Delft) and modified by Jasper A. Vrugt               ##
    ## DREAM-Suite                                                                      ##
    ## ################################################################################## ##

    # Unpack nuisance variables
    s0, s1, ba, xi, mu1, fi1, fi2, fi3, fi4, K, lba = nuisvar

    n = len(y_meas)                                             # Number of observations
    
    # Initialize output variables
    loglik = -np.inf * np.ones((n,1))
    std_e = eps_n = f_eps_n = np.zeros((n,1))
    Y_r = None

    y_bc = y_sim * np.minimum(10, np.exp(mu1 * y_sim))          # Bias-corrected simulation
    y_bc2 = ((y_bc + K) ** lba - 1) / lba                       # Box-Cox bias-corrected simulation
    e_bc = ((y_meas + K) ** lba - 1) / lba - y_bc2              # Box-Cox raw residuals
    
    # Standardization based on the method parameter
    if method == 0:
        std_e = sigma                       # User-specified
    elif method == 1:
        std_e = s0 * np.ones(n)             # Constant s0
    elif method == 2:
        std_e = np.std(e_bc) * np.ones(n)   # Constant based on raw residuals
    elif method == 3:
        std_e = s0 * np.ones(n)             # Constant s0, no s1
    elif method == 4:
        std_e = s0 * np.ones(n)             # Constant s0, no s1
    elif method == 5:
        std_e, _, exitflag = s1_phantom(s0, e_bc, y_bc)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 6:
        std_e = s0 + s1 * y_bc              # Heteroscedastic model with s1 fixed
    elif method == 7:
        std_e, _, exitflag = s1_phantom(s0, e_bc, y_bc)
        if exitflag != 1:
            return loglik, std_e, eps_n, f_eps_n, Y_r
    elif method == 8:
        std_e = s0 + s1 * y_bc              # Heteroscedastic model with s1 estimated
    
    # Autoregressive coefficients for AR filter
    fi_p = np.array([1, -fi1, -fi2, -fi3, -fi4])

    # Likelihood estimation: For model training
    if iflag == 'est':
        # eps needs to be a vector with 1 dimension!
        # then lfilter is equal to filter in MATLAB
        eps = lfilter(fi_p, 1, e_bc)                            # Partial residuals
        eps_n = eps / std_e                                     # Standardized residuals
        f_eps_n, logf_eps_n = f_SEP(eps_n, ba, xi)              # SEP densities
        loglik = logf_eps_n - np.log(std_e) + (lba - 1) * np.log(y_meas + K)

    # Simulation: Replicates for uncertainty quantification
    elif iflag == 'sim':  
        eps_n = SEPrnd(ba, xi, n, N)                                # Draw nxm matrix of standardized partial residuals from SEP(0,1,beta,xi)
        eps = std_e[:, None] * eps_n                                # nxm matrix of non-standardized partial residuals
        E_r = lfilter([1], fi_p, eps, axis = 0)                     # nxm matrix of Box-Cox residuals
        Y_r = (lba * (y_bc2[:, None] + E_r) + 1) ** (1 / lba) - K   # nxm matrix of replicates of simulation

    return loglik, std_e, eps_n, f_eps_n, Y_r


def Whittle_loglik(y_sim, Meas_info):
    ## ################################################################################## ##
    ## Function computes Whittle's log-likelihood function using spectral densities of    ##
    ## measurements and model output                                                      ##
    ##                                                                                    ##
    ## SYNOPSIS: loglik = Whittle_loglik(fx,Meas_info)                                    ##
    ##  where                                                                             ##
    ##   fx        [input]  REQUIRED: n x 1 vector of simulated model output              ##
    ##   Meas_info [input]  REQUIRED: Measurement structure with measured values          ##
    ##   loglik    [outpt]  Whittle quasi-log-likelihood                                  ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##  Whittle, P. (1957), Curve and Periodogram Smoothing, Journal of the Royal         ##
    ##      Statistical Society, Ser. B, 19, 38-63                                        ##
    ##  Whittle, P. (1962), Gaussian Estimation in Stationary Time Series, Bulletin of    ##
    ##      the International Statistical Institute, 39, 105-129                          ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2007                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    # Residual vector (difference between measurements and simulated output)
    e = Meas_info['Y'] - y_sim      # n x 1 residual vector
    # AR(1) model: estimate AR coefficient phi and partial variance s2_pe
    phi, s2_pe = armcov(e, 1)       # armcov: e_(t) = phi*e_(t-1) + pe_(t)  where pe ~ N(0,s2_pe) is partial residual 
    phi = - phi[1]                  # to be consistent with MATLAB
    # n_half: Half the length of the data, rounded down
    n_half = (Meas_info['n'] - 1) // 2
    # Compute the periodogram (spectral density)
    per_meas = np.abs(np.fft.fft(Meas_info['Y'])) ** (2 / (2 * np.pi * Meas_info['n']))
    # Taking relevant part of periodogram for the measured data
    per_meas = per_meas[1 : n_half + 1]
    # Generate frequency indices (sine and cosine components)
    id = np.arange(1, n_half + 1) * (2 * np.pi) / Meas_info['n']
    sin_id = np.sin(id)
    cos_id = np.cos(id)

    # Periodogram of simulated values (spectral density)
    per_sim = np.abs(np.fft.fft(y_sim)) ** (2 / (2 * np.pi * Meas_info['n']))
    # Taking relevant part of periodogram for the simulated data
    per_sim = per_sim[1 : n_half + 1]
    # Compute autoregressive components
    I_ar = phi * sin_id
    R_ar = phi * cos_id
    f_ar = (1 - R_ar) ** 2 + I_ar ** 2
    # Compute f_spec (spectral density of the autoregressive error model)
    f_spec = s2_pe / (2 * np.pi) * 1 / f_ar
    # Total spectral density (model + AR(1) noise)
    per_tot = per_sim + f_spec   
    # Ratio of spectral densities of measured data and total model
    y_f = per_meas / per_tot    
    # Identify which elements of per_tot are positive (avoid log(0) issues)
    id_valid = per_tot > 0
    # Whittle's log-likelihood
    loglik = -np.sum(np.log(per_tot[id_valid])) - np.sum(y_f[id_valid])
    
    return loglik


def SGTrnd(mu, sigma, lbd, p, q, n, m):     #def SGTrnd(mu, sigma, lbd, p, q, *args):
    ## ################################################################################## ##
    ## Function draws samples, R, from the skewed generalized Student's t (SGT)           ##
    ## distribution, R ~ SGT(µ, σ, λ, p, q)                                               ##
    ##                                                                                    ##
    ## SYNOPSIS: R = SGTrnd(mu, sigma, lbd, p, q, n, m)                                   ##
    ##  where                                                                             ##
    ##   mu        [input]  REQUIRED: Mean of the SGT distribution (µ ∈ R)                ##
    ##   sigma     [input]  REQUIRED: Standard deviation of the SGT distribution (σ > 0)  ##
    ##   lbd       [input]  REQUIRED: Skewness of the SGT distribution (-1 < λ < 1)       ##
    ##   p         [input]  REQUIRED: Controls kurtosis, p > 0                            ##
    ##   q         [input]  REQUIRED: Controls kurtosis, q > 2                            ##
    ##   n         [input]  REQUIRED: Number of rows of R                                 ##
    ##   m         [input]  REQUIRED: NUmber of columns of R                              ##
    ##   R         [outpt]  Samples from SGT(µ, σ, λ, p, q)                               ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2017                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    n = int(n)
    m = int(m)
    # Check for valid parameter values
    if sigma <= 0 or lbd <= -1 or lbd >= 1 or p <= 0 or q <= 2 or np.isnan(mu) \
        or np.isnan(sigma) or np.isnan(lbd) or np.isnan(p) or np.isnan(q):
        return np.nan * np.ones(n, m)

    # Generate uniform random numbers (u) between 0 and 1
    u = np.random.rand(n, m)
    # Use the inverse CDF (quantile function) to generate samples
    R = SGTinv(u, mu, sigma, lbd, p, q)
    
    return R


def SGTinv(P, mu, sigma, lbd, p, q):
    ## ################################################################################## ##
    ## Quantile function of the skewed generalized Student's t (SGT) distribution         ##
    ## [= inverse CDF], X ~ SGT^{-1}(P, µ, σ, λ, p, q)                                    ##
    ##                                                                                    ##
    ## SYNOPSIS: R = SGTinv(P, mu, sigma, lbd, p, q)                                      ##
    ##  where                                                                             ##
    ##   P         [input]  REQUIRED: Array of probabilities, P ∈ [0,1]                   ##
    ##   mu        [input]  REQUIRED: Mean of the SGT distribution (µ ∈ R)                ##
    ##   sigma     [input]  REQUIRED: Standard deviation of the SGT distribution (σ > 0)  ##
    ##   lbd       [input]  REQUIRED: Skewness of the SGT distribution (-1 < λ < 1)       ##
    ##   p         [input]  REQUIRED: Controls kurtosis, p > 0                            ##
    ##   q         [input]  REQUIRED: Controls kurtosis, q > 2                            ##
    ##   X         [outpt]  Array with inverse SGT at the values in P                     ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2017                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    # Initialize output array X
    if np.issubdtype(P.dtype, np.floating):
        X = np.zeros_like(P, dtype=np.float32) if np.issubdtype(P.dtype, np.floating) else np.zeros_like(P)
    else:
        X = np.zeros_like(P)

    # Parameters
    q_max = 200

    # For q > 200, the GST becomes equivalent to GED
    # The inverse cdf of 0 is -Inf, and the inverse cdf of 1 is Inf
    # X[P == 1] = np.inf
    # X[P == 0] = -np.inf

    # Find values between zero and one
    # ii = np.where((P > 0) & (P < 1))[0]

    # Special case when lambda = 0 and x = 0.5, SGTinv(1/2) = 0
    # k = np.where((P == 0.5) & (lbd == 0))[0]
    # if k.size > 0:
    #    X[k] = 0
    #    ii = np.setdiff1d(ii, k)  # Remove the indices where k occurs from ii

    # Check whether SGT or GED
    if q <= q_max:  # SGT distribution
        theta = (sp.beta(1/p, q/p).mean() / np.sqrt((1 + 3 * lbd**2) * sp.beta(1/p, q/p).mean() * \
                                                    sp.beta(3/p, (q - 2)/p).mean() - 4 * lbd**2 * \
                                                        sp.beta(2/p, (q - 1)/p).mean()**2))
        mu = mu - 2 * theta * lbd * sigma * sp.beta(2/p, (q - 1)/p).mean() / sp.beta(1/p, q/p).mean()
        # Betaincinv function in Python equivalent to MATLAB's betaincinv but first input arg becomes last input arg
        # w = sp.betaincinv(1/p, q/p, (2 * P[ii] - (1 - lbd)) / (lbd + np.sign(P[ii] - 0.5 * (1 - lbd))))
        w = sp.betaincinv(1/p, q/p, (2 * P - (1 - lbd)) / (lbd + np.sign(P - 0.5 * (1 - lbd))))
        # X[ii] = theta * sigma * (np.sign(P[ii] - 0.5 * (1 - lbd)) + lbd) / (w**(-1) - 1)**(1/p) + mu
        X = theta * sigma * (np.sign(P - 0.5 * (1 - lbd)) + lbd) / (w**(-1) - 1)**(1/p) + mu

    else:  # GED distribution
        theta = sp.gamma(1/p) / np.sqrt((1 + 3*lbd**2) * sp.gamma(1/p) * sp.gamma(3/p) - 4*lbd**2 * sp.gamma(2/p)**2)
        mu = mu - 2 * theta * lbd * sigma * sp.gamma(2/p) / sp.gamma(1/p)
        # gammaincinv function in Python equivalent to MATLAB's gammaincinv but first input arg becomes second input arg
        # w = sp.gammaincinv(1/p, (2 * P[ii] - (1 - lbd)) / (lbd + np.sign(P[ii] - 0.5 * (1 - lbd))))
        w = sp.gammaincinv(1/p, (2 * P - (1 - lbd)) / (lbd + np.sign(P - 0.5 * (1 - lbd))))
        # X[ii] = theta * sigma * (lbd + np.sign(P[ii] - 0.5 * (1 - lbd))) * w**(1/p) + mu
        X = theta * sigma * (lbd + np.sign(P - 0.5 * (1 - lbd))) * w**(1/p) + mu

    # After the fact 
    X[P == 1] = np.inf
    X[P == 0] = -np.inf
    k = np.where((P == 0.5) & (lbd == 0))[0]
    if k.size > 0:
        X[k] = 0

    return X


def SSTrnd(nu, xi, n, m):   #def SSTrnd(nu, xi, *shape):
    ## ################################################################################## ##
    ## Function draws samples, R, from the standardized skewed Student's t (SST)          ##
    ## distribution, R ~ SST(0, 1, ν, ξ)                                                  ##
    ##                                                                                    ##
    ## SYNOPSIS: R = SSTrnd(nu, xi, n, m)                                                 ##
    ##  where                                                                             ##
    ##   nu        [input]  REQUIRED: Degrees of freedom (ν > 2)                          ##
    ##   xi        [input]  REQUIRED: Kurtosis parameter (ξ > 0)                          ##
    ##   n         [input]  REQUIRED: Number of rows of R                                 ##
    ##   m         [input]  REQUIRED: NUmber of columns of R                              ##
    ##   R         [outpt]  nxm matrix with samples from SST(0, 1, ν, ξ, n, m)            ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Scharnagl, B., S.C. Iden, W. Durner, H. Vereecken, and M. Herbst (2015), Inverse ##
    ##       modelling of in situ soil water dynamics: accounting for heteroscedastic,    ##
    ##       autocorrelated, and non-Gaussian distributed residuals, Hydrology and Earth  ##
    ##       System Sciences Discussions, vol. 12, pp. 2155-2199, 2015.                   ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2016                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    if nu <= 2 or xi <= 0 or np.isnan(nu) or np.isnan(xi):
        return np.nan * np.ones(n, m)

    # Compute first and second moments of the SST distribution
    logM = lambda j, nu: (np.log(sp.gamma((j + 1) / 2)) +
                          np.log(sp.gamma((nu - j) / 2)) +
                          (j / 2) * np.log(nu - 2) - 0.5 * np.log(np.pi) - np.log(sp.gamma(nu / 2)))
    M1 = np.exp(logM(1, nu))
    M2 = np.exp(logM(2, nu))
    # Compute mu_xi and sig_xi
    mu_xi = M1 * (xi - 1 / xi)
    sig_xi = np.sqrt((M2 - M1 ** 2) * (xi ** 2 + xi ** -2) + 2 * M1 ** 2 - M2)
    # Generate random samples from SST distribution
    ST_rnd = t.rvs(nu, size=(n,m)) * np.sqrt((nu - 2) / nu)  # Student's t with scale factor
    # Create random signs with probabilities based on xi
    w = xi / (xi + 1 / xi)
    signrndw = np.sign(np.random.rand(n, m) - w)
    # Compute SST samples
    SST_rnd = -signrndw * np.abs(ST_rnd) / (xi ** signrndw)
    # Normalize to get standardized SST samples
    R = (SST_rnd - mu_xi) / sig_xi

    return R


def SSTinv(P, nu, xi, calc_method = 1):
    ## ################################################################################## ##
    ## Quantile function of the standardized Student's t (SST) distribution               ##
    ## [= inverse CDF], X ~ SST^{-1}(P, ν, ξ)                                             ##
    ##                                                                                    ##
    ## SYNOPSIS: X = SSTinv(P, nu, xi, method)                                            ##
    ##  where                                                                             ##
    ##   P         [input]  REQUIRED: Array of probabilities, P ∈ [0,1]                   ##
    ##   nu        [input]  REQUIRED: Degrees of freedom (ν > 2)                          ##
    ##   xi        [input]  REQUIRED: Kurtosis parameter (ξ > 0)                          ##
    ##   method    [input]  OPTIONAL: Computation method                                  ##
    ##   X         [outpt]  Array with inverse SST at the values in P                     ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Scharnagl, B., S.C. Iden, W. Durner, H. Vereecken, and M. Herbst (2015), Inverse ##
    ##       modelling of in situ soil water dynamics: accounting for heteroscedastic,    ##
    ##       autocorrelated, and non-Gaussian distributed residuals, Hydrology and Earth  ##
    ##       System Sciences Discussions, vol. 12, pp. 2155-2199, 2015.                   ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2016                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    if nu <= 2 or xi <= 0:
        raise ValueError("SSTinv: nu must be greater than 2 and xi must be greater than 0.")
    
    P = np.array(P).flatten()   # Ensure P is a 1D array
    nP = len(P)                 # Number of elements in P
    X = np.full(nP, np.nan)     # Initialize output array with NaNs

    # Inverse CDF calculation (inverse of the CDF at given quantiles P)
    if calc_method == 0:        ## Fast but approximate brute force method
        x_vals = np.linspace(-50, 50, 100000)               # Create a range of x values
        pdf_vals = f_SST(x_vals, nu, xi)                    # Compute PDF at x
        cdf_vals = np.cumsum(pdf_vals) / np.sum(pdf_vals)   # Compute CDF as cumulative sum
        # Remove duplicate points for interpolation
        ii = np.where(np.diff(cdf_vals) > 0)[0] + 1         # Index where CDF changes
        x_vals = x_vals[ii]
        cdf_vals = cdf_vals[ii]
        # Interpolate to find the inverse CDF values for P
        interp_func = interp1d(cdf_vals, x_vals, bounds_error = False, fill_value = "extrapolate")
        X = interp_func(P)                                  # Inverse CDF using interpolation
    
    elif calc_method == 1:      ## Accurate but numerical method
        for i in range(nP):
            # Define the function to find the root of, based on the CDF
            f = lambda x: integrate.quad(lambda x: f_SST(x, nu, xi), -np.inf, x)[0] - P[i]
            # Use fsolve to find the root (inverse CDF)
            X[i] = fsolve(f, 0)[0]

    return X


def SEPrnd(beta_, xi, n, m):            #def SEPrnd(beta_, xi, *args):
    ## ################################################################################## ##
    ## Function draws samples, R, from the standardized exponential power (SEP)           ##
    ## distribution, R ~ SEP(0, 1, ß, ξ)                                                  ##
    ##                                                                                    ##
    ## SYNOPSIS: R = SEPrnd(beta_, xi, n, m)                                              ##
    ##  where                                                                             ##
    ##   beta_     [input]  REQUIRED: Skewness parameter (-1 < ß <= 1)                    ##
    ##   xi        [input]  REQUIRED: Kurtosis parameter (ξ > 0)                          ##
    ##   n         [input]  REQUIRED: Number of rows of R                                 ##
    ##   m         [input]  REQUIRED: NUmber of columns of R                              ##
    ##   R         [outpt]  nxm matrix with samples from SEP(ß, ξ, n, m)                  ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##  Schoups, G., and J. A. Vrugt (2010), A formal likelihood function for parameter   ##
    ##       and predictive inference of hydrologic models with correlated,               ##
    ##       heteroscedastic, and non-Gaussian errors, Water Resources Research, 46,      ##
    ##       W10531, doi:10.1029/2009WR008933                                             ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2016                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    # Error checking for input parameters
    if beta_ <= -1 or beta_ > 1 or xi <= 0 or np.isnan(beta_) or np.isnan(xi):
        return np.nan * np.ones(n, m)

    # Compute SEP variables so that mu = 0 and sigma = 1
    A1 = sp.gamma(3*(1 + beta_) / 2)
    A2 = sp.gamma((1 + beta_) / 2)
    M1 = sp.gamma(1 + beta_) / np.sqrt(A1 * A2)
    M2 = 1
    mu_xi = M1 * (xi - 1 / xi)
    sig_xi = np.sqrt((M2 - M1**2) * (xi**2 + 1 / xi**2) + 2 * M1**2 - M2)
    # p parameter for the gamma distribution
    p = 2 / (1 + beta_)
    # Step 1: Generate random variates from gamma distribution
    # CHECK (= gamrnd) MATLAB/Python
    # grnd = sp.gamma.rvs(1/p, size=args)
    grnd = np.random.gamma(1/p, 1, size=(n,m))
    # Step 2: Generate random signs (+1 or -1) with equal probability
    signrnd = np.sign(np.random.rand(n, m) - 0.5)
    # Step 3: Compute random variates from EP(0,1,beta_)
    EP_rnd = signrnd * (np.abs(grnd) ** (1/p)) * np.sqrt(sp.gamma(1/p)) / np.sqrt(sp.gamma(3/p))
    # Step 4: Generate random signs (+1 or -1) with probability 1-w and w
    w = xi / (xi + 1 / xi)
    signrndw = np.sign(np.random.rand(n, m) - w)
    # Step 5: Compute random variates from SEP(mu_xi, sig_xi, xi, beta_)
    SEP_rnd = -signrndw * np.abs(EP_rnd) / (xi ** signrndw)
    # Standardize the results
    R = (SEP_rnd - mu_xi) / sig_xi
    
    return R


def SEPinv(P, beta_, xi, calc_method = 1):
    ## ################################################################################## ##
    ## Quantile function of the standardized skewed exponential power (SEP) distribution  ##
    ## [= inverse CDF], X ~ SEP^{-1}(P, ß, ξ)                                             ##
    ##                                                                                    ##
    ## SYNOPSIS: X = SEPinv(P, beta_, xi, method)                                         ##
    ##  where                                                                             ##
    ##   P         [input]  REQUIRED: Array of probabilities, P ∈ [0,1]                   ##
    ##   beta_     [input]  REQUIRED: Skewness parameter (-1 < ß <= 1)                    ##
    ##   xi        [input]  REQUIRED: Kurtosis parameter (ξ > 0)                          ##
    ##   method    [input]  OPTIONAL: Computation method                                  ##
    ##   X         [outpt]  Array with inverse SEP at the values in P                     ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Scharnagl, B., S.C. Iden, W. Durner, H. Vereecken, and M. Herbst (2015), Inverse ##
    ##       modelling of in situ soil water dynamics: accounting for heteroscedastic,    ##
    ##       autocorrelated, and non-Gaussian distributed residuals, Hydrology and Earth  ##
    ##       System Sciences Discussions, vol. 12, pp. 2155-2199, 2015.                   ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2016                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##
   
    if beta_ <= -1 or beta_ > 1 or xi <= 0:
        raise ValueError("SEPinv: beta_ must be between -1 and 1, and xi must be greater than 0.")
    
    P = np.array(P).flatten()                   # Ensure P is a 1D array
    nP = len(P)                                 # Number of elements in P
    X = np.full(nP, np.nan)                     # Initialize output array with NaNs
    
    if calc_method == 0:    ## Faster but approximate brute force method
        x_vals = np.linspace(-50, 50, 100000)   # Create a range of x values
        pdf_vals = f_SEP(x_vals, beta_, xi)     # Compute PDF at x
        cdf_vals = np.cumsum(pdf_vals) / np.sum(pdf_vals)  # Compute CDF as cumulative sum
        # Remove duplicate points for interpolation
        ii = np.where(np.diff(cdf_vals) > 0)[0] + 1  # Index where CDF changes
        x_vals = x_vals[ii]
        cdf_vals = cdf_vals[ii]
        # Interpolate to find the inverse CDF values for P
        interp_func = interp1d(cdf_vals, x_vals, bounds_error=False, fill_value="extrapolate")
        X = interp_func(P)                      # Inverse CDF using interpolation
    
    elif calc_method == 1:  ## Accurate, but slower numerical solution
        for i in range(nP):
            # Define the function to find the root of, based on the CDF
            f = lambda x: integrate.quad(lambda x: f_SEP(x, beta_, xi)[0], -np.inf, x)[0] - P[i]
            # Use fsolve to find the root (inverse CDF)
            X[i] = fsolve(f, 0)[0]
    
    return X


def LAPrnd(mu, sigma, n, m):        #def LAPrnd(mu, sigma, *shape):
    ## ################################################################################## ##
    ## Function draws samples, R, from the standardized Laplacian (LAP) distribution      ##
    ## distribution, R ~ LAP(µ, σ)                                                        ##
    ##                                                                                    ##
    ## SYNOPSIS: R = LAPrnd(mu, sigma, n, m)                                              ##
    ##  where                                                                             ##
    ##   mu        [input]  REQUIRED: Mean of the LAP distribution (µ ∈ R)                ##
    ##   sigma     [input]  REQUIRED: Standard deviation of the LAP distribution (σ > 0)  ##
    ##   n         [input]  REQUIRED: Number of rows of R                                 ##
    ##   m         [input]  REQUIRED: NUmber of columns of R                              ##
    ##   R         [outpt]  nxm matrix with samples from LAP(µ, σ)                        ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2016                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##
    
    if sigma <= 0 or np.isnan(sigma):
        raise ValueError("LAPrnd: sigma must be greater than 0")
    
    # Generate uniform random numbers in the range [-0.5, 0.5]
    u = np.random.rand(n, m) - 0.5
    # Compute the scale parameter for the Laplace distribution
    b = sigma / np.sqrt(2)
    # Generate Laplace-distributed samples using the formula
    R = mu - b * np.sign(u) * np.log(1 - 2 * np.abs(u))
    
    return R


def LAPinv(P, mu, sigma, calc_method = 2):
    ## ################################################################################## ##
    ## Quantile function of the Laplace (LAP) distribution                                ##
    ## [= inverse CDF], X ~ LAP^{-1}(P, µ, σ)                                             ##
    ##                                                                                    ##
    ## SYNOPSIS: X = LAPinv(P, µ, σ, method)                                              ##
    ##  where                                                                             ##
    ##   P         [input]  REQUIRED: Array of probabilities, P ∈ [0,1]                   ##
    ##   mu        [input]  REQUIRED: Mean of the SGT distribution (µ ∈ R)                ##
    ##   sigma     [input]  REQUIRED: Standard deviation of the SGT distribution (σ > 0)  ##
    ##   method    [input]  OPTIONAL: Computation method                                  ##
    ##   X         [outpt]  Array with inverse LAP at the values in P                     ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2016                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##
    
    if sigma <= 0:
        raise ValueError("LAPinv: sigma must be greater than 0")
    
    # Define the scale parameter 'b' based on sigma
    b = sigma / np.sqrt(2)
    # Define the Laplace CDF function
    def F_lap(x):
        return 0.5 + 0.5 * np.sign(x - mu) * (1 - np.exp(-np.abs(x - mu) / b))
    
    # Convert P to a column vector
    P = np.array(P).flatten()
    nP = len(P)
    X = np.full(nP, np.nan)
    
    if calc_method == 0:    ## Approximate but brute force solution (interpolation)
        x = np.arange(-10, 10, 0.001)       # Fine grid of x values
        cdf_x = F_lap(x)                    # Calculate the CDF for each x
        unique_x = np.unique(x)             # Remove duplicate values
        cdf_x = F_lap(unique_x)
        X = np.interp(P, cdf_x, unique_x)   # Inverse CDF through interpolation
    
    elif calc_method == 1:  ## Exact but slower numerical solution (root finding)
        for i in range(nP):
            f = lambda x: F_lap(x) - P[i]       # Root-finding problem
            X[i] = optimize.fsolve(f, 0)[0]     # Find the root (zero point)
    
    elif calc_method == 2:  ## Exact analytic solution (quantile function)
        def Fi_lap(p):
            return mu - b * np.sign(p - 0.5) * np.log(1 - 2 * np.abs(p - 0.5))
        
        X = Fi_lap(P)  # Return the quantiles for the given probabilities
    
    return X


def f_SGT(a, lbd, p, q):
    ## ################################################################################## ##
    ## Function evaluates the density of the standardized skewed generalized Student's t  ##
    ## distribution at the samles of array a, pdf = SGT(a, 0, 1, λ, p, q)                 ##
    ##                                                                                    ##
    ## SYNOPSIS: pdf = f_SGT(a, lbd, p, q)                                                ##
    ##  where                                                                             ##
    ##   a         [input]  REQUIRED: Array of values for which we compute SGT density    ##
    ##   lbd       [input]  REQUIRED: Skewness of the SGT distribution (-1 < λ < 1)       ##
    ##   p         [input]  REQUIRED: Controls kurtosis, p > 0                            ##
    ##   q         [input]  REQUIRED: Controls kurtosis, q > 2                            ##
    ##   pdf       [outpt]  Density of SGT distribution at samples of a                   ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2017                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##
    
    # Check for valid inputs
    if p * q < 2:
        raise ValueError("f_SGT ERROR: The product of p and q must exceed two.")
    
    # Compute the scaling factor kappa_lpq and mean shift mu_lpq
    kappa_lpq = sp.beta(1/p, q/p) / np.sqrt((1 + 3*lbd**2) * sp.beta(1/p, q/p) * sp.beta(3/p, (q-2)/p) - 4*lbd**2 * sp.beta(2/p, (q-1)/p)**2)
    mu_lpq = 2 * kappa_lpq * lbd * sp.beta(2/p, (q-1)/p) / sp.beta(1/p, q/p)
    # Shift the input `a` by the mean `mu_lpq`
    eps_sgt = a + mu_lpq
    # Compute the SGT PDF
    pdf = p / (2 * kappa_lpq * sp.beta(1/p, q/p)) * (1 + (np.abs(eps_sgt) / (kappa_lpq * (1 + lbd * np.sign(eps_sgt))))**p)**(-(q+1)/p)
    
    return pdf


def f_SEP(a, beta_, xi):
    ## ################################################################################## ##
    ## Function evaluates the density of the standardized skewed exponential power (SEP)  ##
    ## distribution at the samples of array a, pdf = SEP(a, 0, 1, ß, ξ)                   ##
    ##                                                                                    ##
    ## SYNOPSIS: pdf, logpdf = f_SEP(a, beta_, xi)                                        ##
    ##  where                                                                             ##
    ##   a         [input]  REQUIRED: Array of values for which we compute SGT density    ##
    ##   beta_     [input]  REQUIRED: Skewness parameter (-1 < ß <= 1)                    ##
    ##   xi        [input]  REQUIRED: Kurtosis parameter (ξ > 0)                          ##
    ##   pdf       [outpt]  Density of SEP distribution at samples of a                   ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2017                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##
    
    # Compute SEP parameters
    A1 = sp.gamma(3 * (1 + beta_) / 2)
    A2 = sp.gamma((1 + beta_) / 2)
    Cb = (A1 / A2) ** (1 / (1 + beta_))
    Wb = np.sqrt(A1) / ((1 + beta_) * (A2 ** 1.5))
    M1 = sp.gamma(1 + beta_) / np.sqrt(A1 * A2)
    M2 = 1
    mu_xi = M1 * (xi - 1 / xi)
    sig_xi = np.sqrt((M2 - M1 ** 2) * (xi ** 2 + 1 / xi ** 2) + 2 * M1 ** 2 - M2)
    # Skewed SEP transformation
    a_SEP = mu_xi + sig_xi * a
    a_skew = a_SEP / xi ** (np.sign(a_SEP))
    # Small constant norm_cst to handle uniform case
    norm_cst = max(np.finfo(np.float64).tiny,(2 * sig_xi / (xi + 1 / xi)) * Wb)        # --> converges to Wb
    if Cb < 1e-100:                                         # Uniform distribution case
        # pdf = np.zeros_like(a)                            # Initialize density to zero
        pdf = np.finfo(np.float64).tiny * np.ones_like(a)   # Initialize density         
        ii = np.abs(a_skew) <= (1 / (2 * norm_cst))         # Range for uniform PDF
        pdf[ii] = norm_cst                                  # PDF equals Wb within the specified range
        logpdf = np.log(pdf)                                # Log of PDF
    else:  # SEP distribution with skewness
        logpdf = np.log(norm_cst) - Cb * np.abs(a_skew) ** (2 / (1 + beta_))
        pdf = np.exp(logpdf)
    
    return pdf, logpdf
    # Two output variables --> consequences for root finding: [0] must be added (SEPinv)


def f_SST(a, nu, xi, method = 1):
    ## ################################################################################## ##
    ## Function evaluates the density of the standardized skewed Student's t (SST)        ##
    ## distribution at the samples of array a, pdf = SST(a, 0, 1, ν, ξ)                   ##
    ##                                                                                    ##
    ## SYNOPSIS: pdf = f_SST(a, ν, xi)                                                    ##
    ##  where                                                                             ##
    ##   a         [input]  REQUIRED: Array of values for which we compute SST density    ##
    ##   nu        [input]  REQUIRED: Degrees of freedom (ν > 2)                          ##
    ##   xi        [input]  REQUIRED: Kurtosis parameter (ξ > 0)                          ##
    ##   method    [input]  OPTIONAL: Computation method                                  ##
    ##   pdf       [outpt]  Density of SST distribution at samples of a                   ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Scharnagl, B., S.C. Iden, W. Durner, H. Vereecken, and M. Herbst (2015), Inverse ##
    ##       modelling of in situ soil water dynamics: accounting for heteroscedastic,    ##
    ##       autocorrelated, and non-Gaussian distributed residuals, Hydrology and Earth  ##
    ##       System Sciences Discussions, vol. 12, pp. 2155-2199, 2015.                   ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, Feb. 2017                                            ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##
    
    # Ensure `method` is valid (0 or 1)
    if method not in [0, 1]:
        raise ValueError("f_SST ERROR: Method must be 0 or 1.")
    
    # Compute first and second moments of SST distribution
    if method == 0:     ## For small nu values (note: not defined for large nu)
        def M(j, nu):
            return (sp.gamma((j+1)/2) * sp.gamma((nu-j)/2) * (nu-2)**(j/2)) \
                / (np.sqrt(np.pi) * sp.gamma(nu/2))
        M1 = M(1, nu)
        M2 = M(2, nu)
    elif method == 1:   ## Log formulation for all nu to avoid numerical overflow
        def logM(j, nu):
            return (sp.gammaln((j+1)/2) + sp.gammaln((nu-j)/2) + (j/2)*np.log(nu-2) \
                    - 0.5*np.log(np.pi) - sp.gammaln(nu/2))
        M1 = np.exp(logM(1, nu))
        M2 = np.exp(logM(2, nu))
    
    # Now compute mu_xi and sig_xi
    mu_xi = M1 * (xi - 1 / xi)              # Matches Scharnagl paper
    sig_xi = np.sqrt((M2 - M1**2) * (xi**2 + xi**-2) + 2 * M1**2 - M2)
    # Standardize a
    a_sst = (mu_xi + sig_xi * a) / (xi ** np.sign(mu_xi + sig_xi * a))
    # Standardized Student's t density with skew
    if method == 0:     ## Method 0, valid for small degrees of freedom
        pdf = (2 * sig_xi * sp.gamma((nu + 1) / 2)) / ((xi + 1 / xi) * sp.gamma(nu / 2) * \
                np.sqrt(np.pi * (nu - 2))) * (1 + 1 / (nu - 2) * a_sst**2) ** (-(nu + 1) / 2)
        logpdf = np.log(pdf)
    elif method == 1:   ## Method 1, valid for all nu > 2
        logpdf = np.log(2 * sig_xi) + sp.gammaln((nu + 1) / 2) - np.log(xi + 1 / xi) - sp.gammaln(nu / 2) \
            - 0.5 * np.log(np.pi * (nu - 2)) - (nu + 1) / 2 * np.log(1 + 1 / (nu - 2) * a_sst**2)
        pdf = np.exp(logpdf)
    
    return pdf


def s1_phantom(s0, e, y):
    ## ################################################################################## ##
    ## Function determines the slope of the measurement error function so that the        ##
    ## studentized residuals have a unit variance                                         ##
    ##                                                                                    ##
    ## SYNOPSIS: std_e, s1, ier = s1_phantom(s0, e, y)                                    ##
    ##  where                                                                             ##
    ##   s0        [input]  REQUIRED: 1xN vector of intercepts measurement error function ##
    ##   e         [input]  REQUIRED: nxN matrix of raw residuals                         ##
    ##   y         [input]  REQUIRED: nx1 matrix of measured y-values                     ##
    ##   std_e     [outpt]  nxN matrix of measurement error standard devations            ##
    ##   s1        [outpt]  1xN vector of slopes of measurement error function            ##
    ##   ier       [outpt]  Exit flag of Nth sample of s0/e/y                             ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, April 2014                                           ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    if e.ndim == 1:
        n = e.shape[0]
        e = e.reshape(-1,1)
        N = 1
    else:
        n, N = e.shape  # Number of elements (residuals), Number of residual vectors
    
    s1 = np.nan * np.ones(N)            # Initialize slopes
    ier = np.nan * np.ones(N)           # Initialize exitflags
    std_e = np.full((n, N), np.nan)     # Initialize standard deviations
    # Make N copies of y if it is a column vector
    if y.ndim == 1:
        if N > 1:
            y = np.tile(y, (N, 1)).T
        else:
            y = y.reshape(-1,1)

    if np.isscalar(s0):                 # isinstance(s0, (int, float)):
        s0 = np.array([s0])
      
    for ii in range(N):
        # Use fsolve for root-finding: Newton's method
        s1[ii], info, ier[ii], msg = fsolve(lambda x: eval_func(x, s0[ii], e[:, ii], y[:, ii]), 0.2, full_output = True)
        # This gives a warning, why I do not understand; all inuts and outputs have the correct dimensions
        # DeprecationWarning: Conversion of an array with ndim > 0 to a scalar is deprecated, and will error in future. 
        # Ensure you extract a single element from your array before performing this operation
        if s1[ii] < 0 or ier[ii] !=1:
            # If slope s1 is improper, use fminsearch (Nelder-Mead method)
            # s1[ii], _, exflag = fminsearch(lambda x: get_opt(x, s0[ii], e[:, ii], y[:, ii]), 0.2)
            res = minimize(get_opt, 0.2, args = (s0[ii], e[:, ii], y[:, ii]), method = 'Nelder-Mead')
            s1[ii] = res.x
            if s1[ii] < 0:
                print("s1_phantom Warning: s1 < 0")
        # Calculate the measurement error standard deviations
        std_e[:, ii] = s1[ii] * y[:, ii] + s0[ii]
    if N == 1:
        if std_e.ndim > 1:  # Check if std_e is multidimensional
            std_e = std_e.flatten()
        if s1.ndim > 1:     # Check if s1 is multidimensional
            s1 = s1.squeeze()

    return std_e, s1, ier


def eval_func(s1, s0, e, y):
    ## ################################################################################## ##
    ## Function computes the distance (d ∈ R) from unity of the studentized variance of   ##
    ## the raw residuals                                                                  ##
    ##                                                                                    ##
    ## SYNOPSIS: d = eval_func(s1, s0, e, y)                                              ##
    ##  where                                                                             ##
    ##   s1        [input]  REQUIRED: Slope of measurement error function                 ##
    ##   s0        [input]  REQUIRED: Intercept of measurement error function             ##
    ##   e         [input]  REQUIRED: nx1 matrix of raw residuals                         ##
    ##   y         [input]  REQUIRED: nx1 matrix of measured y-values                     ##
    ##   d         [outpt]  Distance from unity of studentized variance of raw residuals  ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, April 2014                                           ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    std_e = s1 * y + s0             # Measurement error function
    e_n = e / std_e                 # Normalized residuals
    d = 1 - np.var(e_n, ddof=1)     # Distance = 1 - variance of studentized residuals

    return d


def get_opt(s1, s0, e, y, w = 100):
    ## ################################################################################## ##
    ## Function computes the error (err > 0) of the studentized variance of the raw       ##
    ## residuals. The error equals the sum of the squared distance from one of the        ##
    ## studentized variance and a penalty term for s1 < 0 which increases with magnitude  ##
    ## of the measurement error slope, s1                                                 ##
    ##                                                                                    ##
    ## SYNOPSIS: err = get_opt(s1, s0, e, y)                                              ##
    ##  where                                                                             ##
    ##   s1        [input]  REQUIRED: Slope of measurement error function                 ##
    ##   s0        [input]  REQUIRED: Intercept of measurement error function             ##
    ##   e         [input]  REQUIRED: nx1 matrix of raw residuals                         ##
    ##   y         [input]  REQUIRED: nx1 matrix of measured y-values                     ##
    ##   w         [input]  OPTIONAL: w = 100 is weight for penalty term                  ##
    ##   err       [outpt]  Error of measurement error function [= constraint violation]  ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, April 2014                                           ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    std_e = s1 * y + s0                 # Measurement error function
    e_n = e / std_e                     # Normalized residuals
    d = 1 - np.var(e_n, ddof=1)         # Distance = 1 - variance of studentized residuals
    pen = w * (s1 < 0) * np.abs(s1)     # Penalty if s1 < 0
    err = d**2 + pen                    # Sum of squared error (with penalty)

    return err


def Bayes_pdf(Up, Ufx, idU, Nr, RMSE_map, DREAMPar, Meas_info, Lik_info, p_gam = 0.95):
    ## ################################################################################## ##
    ## Function produces the Bayes posterior [output space] and determines the upper and  ##
    ## lower prediction limits of this distribution                                       ##
    ##                                                                                    ##
    ## SYNOPSIS: fx_mod, fx_tot = Bayes_pdf(Up, Ufx, idU, Nr, RMSE_map, DREAMPar,         ##
    ##                                      Meas_info, Lik_info, p_gam = 0.95)            ##
    ##  where                                                                             ##
    ##   Up        [input]  REQUIRED: M x d matrix M unique posterior parameter vectors   ##
    ##   Ufx       [input]  REQUIRED: M x n matrix of corresponding simulations           ##
    ##   idU       [input]  REQUIRED: Vector with indices of unique parameter vectors     ##
    ##   Nr        [input]  REQUIRED: Vector with # replicates unique parameter vectors   ##
    ##   RMSE_map  [input]  REQUIRED: Root Mean Square Error of map solution              ##
    ##   DREAMPar  [input]  REQUIRED: Dictionary with algorithmic variables               ##
    ##   Meas_info [input]  REQUIRED: Dictionary with measurement information (fitting)   ##
    ##   Lik_info  [input]  REQUIRED: Dictionary with information for likelihood function ##
    ##   p_gam     [input]  OPTIONAL: Prediction level (default 0.95)                     ##
    ##   fx_mod    [outpt]  Nxm matrix of simulated output of N posterior samples         ##
    ##   fx_tot    [outpt]  Nx2 matrix lower (def: 2.5%) & upper (def: 97.5%) pred. limts ##
    ##                                                                                    ##
    ## For more information please check                                                  ##
    ##   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of ##
    ##       distribution-adaptive likelihood functions: Generalized and universal        ##
    ##       likelihood functions, scoring rules and multi-criteria ranking, Journal of   ##
    ##       Hydrology, 615, Part B, 2022, https://doi.org/10.1016/j.jhydrol.2022.128542  ##
    ##                                                                                    ##
    ## © Written by Jasper A. Vrugt, April 2014                                           ##
    ## Los Alamos National Laboratory 			        	                              ##
    ##                                                                                    ##
    ## ################################################################################## ##

    # Define percentiles for prediction limits
    gam1 = 50 - 50*p_gam # 100 * (1 - p_gam) / 2        # Lower percentile
    gam2 = 50 + 50*p_gam # 100 * (1 + p_gam) / 2        # Upper percentile

    np.random.seed(1 + round(100 * np.random.rand()))   # Random seed
    
    if Up.ndim == 1:
        npars = 1
    else:
        npars = Up.shape[0]                             # How many parameter vectors?
    Ns = np.sum(Nr)                                     # Total number of simulations expected
    fx_mod = np.nan * np.ones((Ns, Meas_info['n']))     # Initial model output uncertainty

    # General code for different likelihood functions
    name_lik = Lik_info['name_lik_func']                # Name of likelihood 
    if DREAMPar['lik'] in {13, 14, 16, 17, 44, 45}:     # Likelihood function cases
        par = Lik_info['fpar'].copy()                   # Initialize fixed values
        for ii in range(npars):
            if Up.ndim == 1:
                par[Lik_info['id_vpar']] = Up               # Variable parameters
                nuisvar = par[Lik_info['id_nuisvar']]       # Nuisance variables
                exec_scope = {                              # Executive dictionary
                    'nuisvar': nuisvar,                     # Input variable
                    'Lik_info': Lik_info,                   # The Lik_info dictionary with the function string
                    'Ufx': Ufx,                             # Output of candidate point
                    'Meas_info': Meas_info,                 # Dictionary of measurement data/information
                    'Nr': Nr,                               # Number of replicates [simulations] requested
                    '_': None,                              # Return argument 1
                    '__': None,                             # Return argument 2
                    '___': None,                            # Return argument 3
                    '____': None,                           # Return argument 4
                    'Y_r': None,                            # Return argument: Ns[ii]xn matrix of simulated outputs
                    }
            else:
                par[Lik_info['id_vpar']] = Up[ii, :DREAMPar['d']]
                nuisvar = par[Lik_info['id_nuisvar']]   # Nuisance variables
                exec_scope = {                          # Executive dictionary
                    'nuisvar': nuisvar,                     # Input variable
                    'Lik_info': Lik_info,                   # The Lik_info dictionary with the function string
                    'Ufx': Ufx,                             # Output of candidate point
                    'Meas_info': Meas_info,                 # Dictionary of measurement data/information
                    'Nr': Nr,                               # Number of replicates [simulations] requested
                    'ii': ii,                               # index of Nr vector    
                    '_': None,                              # Return argument 1
                    '__': None,                             # Return argument 2
                    '___': None,                            # Return argument 3
                    '____': None,                           # Return argument 4
                    'Y_r': None,                            # Return argument: Ns[ii]xn matrix of simulated outputs
                }

            exec_scope[name_lik] = globals()[name_lik]  # Add name of likelihood function
            exec(Lik_info['stringF'], exec_scope)       # Evaluate likelihood function
            ii_id = (idU == ii+1)                       # Indices of replicates
            Y_r = exec_scope['Y_r']
            fx_mod[ii_id, :] = Y_r.T                    # Store replicates 

        if npars == 1:
            print(f"Bayes_PDF: MAP simulation only was used to calculate prediction uncertainty "
                  f"with {Lik_info['name_lik_func']} function!!!")
        else:
            print(f"Bayes_PDF: {Ns} simulations of {npars} unique parameter vectors were used to "
                  f"compute prediction uncertainty with {Lik_info['name_lik_func']} function!!!")
    else:   # Other cases (non-specific likelihood functions)
            # Assuming Meas_info is a dictionary-like object in Python
        if Meas_info['Sigma'] != None:
            if np.size(Meas_info['Sigma']) > 0:
                std_e = Meas_info['Sigma']
            else:
                print("Bayes_PDF WARNING: Field 'Sigma' of structure Meas_info is empty: I use the RMSE_map")
                std_e = RMSE_map
        else:
            std_e = RMSE_map

        # Now make n copies of std_e if it is a scalar
        if np.isscalar(std_e):
            std_e = std_e * np.ones((Meas_info['n'],1))

        # Generate replicates
        for ii in range(npars):
            if np.isscalar(Nr) and npars == 1:
                e_n = np.random.normal(0, 1, (Meas_info['n'], Nr))      # Standard normal residuals
                Uft = Ufx.T
            else:    
                e_n = np.random.normal(0, 1, (Meas_info['n'], Nr[ii]))  # Standard normal residuals
                Uft = Ufx[ii, :Meas_info['n']].T

            #e_r = std_e[:,None] * e_n                                  # Raw residuals
            e_r = std_e * e_n                                           # Raw residuals
            if np.isscalar(Nr) and npars == 1:
                fx = Uft[:, np.newaxis] + e_r
            else:    
                fx = Uft[:, None] + e_r                                 # Replicate simulations
            
            ii_id = (idU == ii+1)                                       # Indices of replicates
            fx_mod[ii_id, :] = fx.T                                     # Store replicates

    fx_tot = np.percentile(fx_mod, [gam1, gam2], axis = 0).T            # Compute 100*[(1-p_gam)/2, (1+p_gam)/2] prediction limits

    if Meas_info['n'] == 1:
        fx_tot = fx_tot.T                                           # Unusual case of only a single measurement

    return fx_tot, fx_mod


def armcov(e, order):
    ## ################################################################################## ##
    ## Autoregressive parameter estimation via modified covariance method                 ##
    ##                                                                                    ##
    ## SYNOPSIS: phi, s2_pe = armcov(e, order)                                            ##
    ##  where                                                                             ##
    ##   e         [input]  REQUIRED: nx1 vector of residuals                             ##
    ##   order     [input]  REQUIRED: order of the AR model: 1 for AR(1)                  ##
    ##   phi       [outpt]  AR coefficients (= 1 + order)                                 ##
    ##   s2_pe     [outpt]  Partial variance of the remaining noise                       ##
    ## ################################################################################## ##
    
    # Length of the residual vector
    n = len(e)
    
    # Compute the autocovariance sequence of the residuals
    gamma = np.correlate(e, e, mode='full') / n
    gamma = gamma[n-1:]  # Keep only non-negative lags

    # Autocovariance matrix for AR model
    R = np.zeros((order, order))
    for i in range(order):
        for j in range(order):
            R[i, j] = gamma[abs(i - j)]  # Using the covariance lag structure

    # Autocovariance vector for the AR model
    r = gamma[1:order+1]

    # Solve for the AR coefficients using the Yule-Walker equations
    phi = np.linalg.solve(R, r)
    
    # Compute the partial variance s2_pe
    s2_pe = gamma[0] - np.dot(phi, r)  # Variance of the residuals
    
    # Append 1 for the first coefficient (AR(0) term)
    phi = np.insert(phi, 0, 1)

    return phi, s2_pe


def extract_names(name):
    ## ################################################################################## ##
    ## Extract variable names univariate and multivariate priors from a function handle   ##
    ##                                                                                    ##
    ## SYNOPSIS: pr_var, idcp = extract_names(name)                                       ##
    ##  where                                                                             ##
    ##   name      [input]  REQUIRED: nx1 vector of residuals                             ##
    ##   pr_var    [outpt]  List of variable names within the function handle             ##
    ##   idcp      [outpt]  Index of closing parenthesis                                  ##
    ## ################################################################################## ##

    # Find the starting index of the first parenthesis '(' and the closing parenthesis ')'
    id_s = name.find('(') + 1
    idcp = [m.start() for m in re.finditer(r'\)', name)]
    id_e = idcp[0] - 1 if idcp else len(name)  # Index of closing parenthesis
    # Extract the part inside the parentheses (variable names separated by commas)
    name_inside_parens = name[id_s:id_e]
    # Find the comma-separated variable names inside the parentheses
    variables = name_inside_parens.split(',')
    # Return variable names
    pr_var = [var.strip() for var in variables]
    
    return pr_var, idcp


def consistent_Bn(Jn, q):
    ## ################################################################################## ##
    ## Newey and West (1987) estimator of the variability matrix [= consistent]           ##
    ##                                                                                    ##
    ## SYNOPSIS: Beta_n = consistent_Bn(Jn, q)                                            ##
    ##  where                                                                             ##
    ##   Jn        [input]  REQUIRED: nxd Jacobian matrix of the parameters               ##
    ##   q         [input]  Maximum correlation/dependence lag entries of score function  ##
    ##   Beta_n    [outpt]  dxd variability matrix according to Newey and West, 1987      ##
    ## ################################################################################## ##
    
    # How many samples, n, and parameters, d
    n, d = Jn.shape
    # Initialize covariance terms at lags 1, 2, ..., q
    Beta_j = np.full((d, d, q), np.nan)
    # Compute the mean variability matrix (= Bn in sandwich form)
    Bn = 1 / n * Jn[:n, :d].T @ Jn[:n, :d]
    # Now compute covariance terms at lags 1, ..., q (Newey and West, 1987)
    for j in range(q):
        t = j + 1
        # Compute mean variability matrix (= Beta_j/n in the paper)
        Beta_j[:, :, j] = 1 / (n - j) * Jn[t:n, :d].T @ Jn[:n - j, :d]
    
    # Now compute the correction term for lags
    dB = np.zeros((d, d))
    for j in range(q):
        dB += (1 - (j / (q + 1))) * (Beta_j[:, :, j] + Beta_j[:, :, j].T)
    
    # Compute the average variability matrix (= Beta_n in the paper)
    Beta_n = Bn + dB
    
    return Beta_n


def create_duplicate_with_dot(word):
    # Split the word at the dot
    parts = word.split('.')
    
    # If the word has only one part, create the duplicate with a dot
    if len(parts) == 1:
        return f"{word}.{word}"
    # If the word already has two identical parts separated by a dot, return the same word
    elif len(parts) == 2 and parts[0] == parts[1]:
        return word
    else:
        # In case it's not the expected format, just return the word as is
        return word
    
################################### End Miscellaneous functions ####################################

#################################### Single chain diagnostics ######################################

def coda(draws, vnames = None, info = None, fid = None):
    ## ################################################################################## ##
    ## Single chain covergence diagnostics [did not finish benchmarking against MATLAB]   ##
    ##                                                                                    ##
    ## SYNOPSIS: result = coda(draws)                                                     ##
    ##  where                                                                             ##
    ##   draws     [input]  REQUIRED: Nxd matrix of N chain samples of d parameters       ##
    ##   vnames    [input]  OPTIONAL: 1xd array with names (strings) of the parameters    ##
    ##   info      [input]  OPTIONAL: additional input information                        ##
    ##   fid       [input]  OPTIONAL: file identified for writing                         ##
    ##   result    [outpt]  Dictionary with Raftery-Lewis, ACF and Geweke diagnosteics    ##
    ## ################################################################################## ##

    num_draws = len(draws)
    
    if num_draws < 50:
        raise ValueError('coda: at least 50 draws are required')
    
    # Set defaults for options
    q, r, s, p1, p2 = 0.025, 0.01, 0.95, 0.2, 0.5
    pflag = 1 if fid is None else 0  # Decide if we need to print to a file or not
    # Handle info parameter (dictionary of options)
    if info is not None:
        if not isinstance(info, dict):
            raise ValueError('coda: must supply options as a dictionary')
        
        # Extract values from info dictionary
        q = info.get('q', q)
        r = info.get('r', r)
        s = info.get('s', s)
        p1 = info.get('p1', p1)
        p2 = info.get('p2', p2)
    
    # Handle variable names
    nflag = False
    if vnames is not None:
        nflag = True
        if len(vnames) != draws.shape[1]:
            raise ValueError('Wrong number of variable names in coda -- check vnames argument')
        Vname = vnames[:16]  # Truncate variable names to 16 characters
    
    # SACF diagnostics
    nlag = 25
    aout = np.zeros((nlag, draws.shape[1]))
    for i in range(draws.shape[1]):
        aout[:, i], _ = sacf(draws[:, i], nlag)
    
    aprt = np.column_stack([aout[0, :], aout[4, :], aout[9, :], aout[24, :]])

    # Raftery-Lewis diagnostics
    rafout = raftery(draws, q, r, s)
    selected_keys = ['kthin', 'nburn', 'n', 'nmin', 'irl']
    rout = np.array([[x[key] for key in selected_keys] for x in rafout.values() if isinstance(x, dict)])
    # Geweke diagnostics
    geweke = momentg(draws)

    # Tested but I do not consider the results in the diagnostics file
    # Split sample and Geweke chi-squared test
    # nobs1 = int(p1 * num_draws)
    # nobs2 = int(p2 * num_draws)
    # draws1 = draws[:nobs1, :]
    # draws2 = draws[-nobs2:, :]
    # res1 = momentg(draws1)
    # res2 = momentg(draws2)
    # resapm = apm(res1, res2)

    result = {
        'q': q,
        'r': r,
        's': s,
        'nvar': draws.shape[1],
        'meth': 'coda',
        'ndraw': num_draws
    }
    
    # Storing Raftery-Lewis diagnostics results
    for i in range(draws.shape[1]):
        result[f'kthin_{i+1}'] = rout[i, 0]
        result[f'nburn_{i+1}'] = rout[i, 1]
        result[f'n_{i+1}'] = rout[i, 2]
        result[f'nmin_{i+1}'] = rout[i, 3]
        result[f'irl_{i+1}'] = rout[i, 4]
    
    # Storing SACF diagnostics results
    for i in range(draws.shape[1]):
        result[f'auto1_{i+1}'] = aprt[i, 0]
        result[f'auto5_{i+1}'] = aprt[i, 1]
        result[f'auto10_{i+1}'] = aprt[i, 2]
        result[f'auto25_{i+1}'] = aprt[i, 3]
    
    # Storing Geweke diagnostics results
    for i in range(draws.shape[1]):
        result[f'nse1_{i+1}'] = geweke[i]['nse1']
        result[f'rne1_{i+1}'] = geweke[i]['rne1']
        result[f'nse2_{i+1}'] = geweke[i]['nse2']
        result[f'rne2_{i+1}'] = geweke[i]['rne2']
        result[f'nse3_{i+1}'] = geweke[i]['nse3']
        result[f'rne3_{i+1}'] = geweke[i]['rne3']
        result[f'pmean_{i+1}'] = geweke[i]['pmean']
        result[f'pstd_{i+1}'] = geweke[i]['pstd']
        result[f'nse_{i+1}'] = geweke[i]['nse']
        result[f'rne_{i+1}'] = geweke[i]['rne']
    
    # Printing results to the console or file
    if pflag == 0:
        # Print diagnostics results (example)
        print(f'MCMC CONVERGENCE diagnostics\nBased on sample size = {num_draws}')
        print(f'Autocorrelation within each parameter chain:')
        # Example for SACF result
        print("SACF diagnostics:")
        print(aprt)
        # Example for Raftery-Lewis diagnostics
        print("Raftery-Lewis diagnostics:")
        print(rout)
        # Example for Geweke diagnostics
        print("Geweke diagnostics:")
        print(geweke)
    
    return result


def prt_coda(results, vnames = None, fid = None):
    ## ################################################################################## ##
    ## Prints convergence diagnostics from DREAM_Suite samples                         ##
    ##                                                                                    ##
    ## SYNOPSIS: rho, bigY = prtcoda(results, vnames, fid)                                ##
    ##  where                                                                             ##
    ##   results   [outpt]  Dictionary returned by coda, raftery, apm, momentg            ##
    ##   vnames    [input]  OPTIONAL: 1xd array with names (strings) of the parameters    ##
    ##   fid       [input]  OPTIONAL: file identified for writing                         ##
    ## ################################################################################## ##

    if not isinstance(results, dict):
        raise ValueError("prt_coda requires a dictionary of results")

    nvar = results['nvar']  # Assuming 'nvar' is the key in the result dictionary

    # Handle variable names
    if vnames is None:
        vnames = [str(i) for i in range(1, nvar + 1)]
    elif len(vnames) != nvar:
        raise ValueError("Wrong # of variable names in coda -- check vnames argument")

    # Depending on the method in results, process different diagnostics
    if results['meth'] == 'raftery':
        # Raftery diagnostics
        print(f"2. RAFTERY-LEWIS DIAGNOSTIC (q = {results['q']:.2f}, r = {results['r']:.2f}, s = {results['s']:.2f})", file=fid)
        print(f"{'Variable':<10}{'Thin':<10.0}{'Burn':<10.0}{'Total(N)':<10.0}{'(Nmin)':<10.0}{'I-stat':<10}", file=fid)
        for i in range(nvar):
            print(f"{vnames[i]:<10}{results[f'kthin_{i+1}']:<10.0}{results[f'nburn_{i+1}']:<10.0}{results[f'n_{i+1}']:<10.0}{results[f'nmin_{i+1}']:<10.0}{results[f'irl_{i+1}']:<10.3f}", file=fid) 

    elif results['meth'] == 'momentg':
        # Geweke NSE, RNE diagnostics
        print("3. GEWEKE DIAGNOSTIC", file=fid)
        print(f"{'Variable':<10}{'Mean':<10}{'std dev':<10}{'NSE iid':<10}{'RNE iid':<10}", file=fid)
        for i in range(nvar):
            print(f"{vnames[i]:<10}{results[f'pmean_{i+1}']:<10.3f}{results[f'pstd_{i+1}']:<10.3f}{results[f'nse_{i+1}']:<10.3f}{results[f'rne_{i+1}']:<10.3f}", file=fid)
        
        print(f"{'Variable':<10}{'NSE 4%':<10}{'RNE 4%':<10}{'NSE 8%':<10}{'RNE 8%':<10}{'NSE 15%':<10}{'RNE 15%':<10}", file=fid)
        for i in range(nvar):
            print(f"{vnames[i]:<10}{results[f'nse1_{i+1}']:<10.3f}{results[f'rne1_{i+1}']:<10.3f}{results[f'nse2_{i+1}']:<10.3f}{results[f'rne2_{i+1}']:<10.3f}{results[f'nse3_{i+1}']:<10.3f}{results[f'rne3_{i+1}']:<10.3f}", file=fid)

    elif results['meth'] == 'apm':
        # Geweke chi-sqr diagnostics
        p1 = results['p1']
        p2 = results['p2']
        print(f"Geweke Chi-squared test for each parameter based on {results['ndraw']} draws", file=fid)
        print(f"First {100*p1:.0f}% versus Last {100*p2:.0f}% of the sample", file=fid)
        print(f"{'Variable':<10}{'Mean':<12}{'N.S.E.':<12}{'Chi-sq Prob':<12}", file=fid)
        for i in range(nvar):
            print(f"{vnames[i]:<12}{results[f'nse_{i+1}']:<12.6f}{results[f'prob_{i+1}']:<12.6f}", file=fid)

    elif results['meth'] == 'coda':
        # CODA diagnostics
        print("-" * 75, file=fid)
        print(f"CHAIN {results['chain_number']} USING {results['ndraw']} SAMPLES", file=fid)
        print("\nA. AUTOCORRELATION", file=fid)
        print(f"{'Variable':<10}{'Lag 1':<10}{'Lag 5':<10}{'Lag 10':<10}{'Lag 25':<10}", file=fid)
        # Create the correct key names by appending the suffixes '_1' or '_2' (or as required)
        for i in range(nvar):
            print(f"{vnames[i]:<10}{results[f'auto1_{i+1}']:<10.3f}{results[f'auto5_{i+1}']:<10.3f}{results[f'auto10_{i+1}']:<10.3f}{results[f'auto25_{i+1}']:<10.3f}", file=fid)

        print(f"B. RAFTERY-LEWIS DIAGNOSTIC (q = {results['q']:.2f}, r = {results['r']:.2f}, s = {results['s']:.2f})", file=fid)
        print(f"{'Variable':<10}{'Thin':<10}{'Burn':<10}{'N':<10}{'Nmin':<10}{'I-stat':<10}", file=fid)
        
        for i in range(nvar):
            print(f"{vnames[i]:<10}{results[f'kthin_{i+1}']:<10.0f}{results[f'nburn_{i+1}']:<10.0f}{results[f'n_{i+1}']:<10.0f}{results[f'nmin_{i+1}']:<10.0f}{results[f'irl_{i+1}']:<10.3f}", file=fid)            

        print("C. GEWEKE DIAGNOSTIC", file=fid)
        print(f"{'Variable':<10}{'Mean':<10}{'std dev':<10}{'NSE iid':<10}{'RNE iid':<10}", file=fid)

        for i in range(nvar):
            print(f"{vnames[i]:<10}{results[f'pmean_{i+1}']:<10.3f}{results[f'pstd_{i+1}']:<10.3f}{results[f'nse_{i+1}']:<10.3f}{results[f'rne_{i+1}']:<10.3f}", file=fid)            

        print(f"{'Variable':<10}{'NSE 4%':<10}{'RNE 4%':<10}{'NSE 8%':<10}{'RNE 8%':<10}{'NSE 15%':<10}{'RNE 15%':<10}", file=fid)
        for i in range(nvar):
            print(f"{vnames[i]:<10}{results[f'nse1_{i+1}']:<10.3f}{results[f'rne1_{i+1}']:<10.3f}{results[f'nse2_{i+1}']:<10.3f}{results[f'rne2_{i+1}']:<10.3f}{results[f'nse3_{i+1}']:<10.3f}{results[f'rne3_{i+1}']:<10.3f}", file=fid)            

    else:
        raise ValueError("results structure not known by prt_coda function")


def sacf(y, m, gflag=0):
    ## ################################################################################## ##
    ## Sample autocorrelation function of MCMC draws                                      ##
    ##                                                                                    ##
    ## SYNOPSIS: rho, bigY = sacf(y, m, gflag)                                            ##
    ##  where                                                                             ##
    ##   y         [input]  REQUIRED: a time-series (need not have mean zero)             ##
    ##   m         [input]  REQUIRED: Number of sample autocorrelations to compute        ##
    ##   gflag     [input]  REQUIRED: Graphing flag (0: graph, 1: no graph, def: 0)       ##
    ##   rho       [outpt]  mx1 vector of sample autocorrelation values                   ##
    ##   bigY      [outpt]  matrix of shifted series for autocorrelation computation      ##
    ## ################################################################################## ##

    n = len(y)
    rho = np.zeros(m)
    npm = n + m
    tmp = np.std(y)
    vary = tmp ** 2
    # Put y in deviations from mean form
    ym = np.mean(y)
    e = y - ym
    # Create the bigY matrix
    bigY = np.tile(np.arange(1, n+1), (2*(m+1), 1)).T
    bigY = bigY.flatten()[:(2*n-1)*(m+1)]
    bigY = bigY.reshape(2*n-1, m+1)
    bigY = bigY[:n, 1:m+1]
    bigY = y[bigY.astype(int)-1]  # adjust for zero-indexing
    # Calculate autocorrelation coefficients
    E = np.tile(e, (m, 1)).T
    rho = np.sum(E * bigY, axis=0) / (n * vary)
    # 2-sigma intervals
    ul = 2 * (1 / np.sqrt(n)) * np.ones(m)
    ll = -2 * (1 / np.sqrt(n)) * np.ones(m)
    # Plot results if gflag is 1
    if gflag == 1:
        plt.bar(np.arange(1, m+1), rho)
        plt.title('Sample Autocorrelation Coefficients')
        plt.xlabel('k-values')
        plt.ylabel('SACF values')
        plt.plot(np.arange(1, m+1), ul, '*r', label = "Upper bound")
        plt.plot(np.arange(1, m+1), ll, '*r', label = "Lower bound")
        plt.legend()
        plt.show()

    return rho, bigY


def raftery(runs, q, r, s):
    ## ################################################################################## ##
    ## Computes the number of draws needed in MCMC to estimate the posterior CDF of the   ##
    ## q-quantile to within +/- r accuracy with probability s                             ##
    ##                                                                                    ##
    ## SYNOPSIS: result = raftery(runs, q, r, s)                                          ##
    ##  where                                                                             ##
    ##   runs      [input]  REQUIRED: A Nxd matrix of draws                               ##
    ##   q         [input]  REQUIRED: Quantile of interest (q ∈ [0,1])                    ##
    ##   r         [input]  REQUIRED: Accuracy                                            ##
    ##   s         [input]  REQUIRED: Probability, s                                      ##
    ##   result    [outpt]  Raftery diagnostics of single chain: nburn, irl, etc.         ##
    ## ################################################################################## ##

    n, nvar = runs.shape
    result = {}
    result['meth'] = 'raftery'
    result['draws'] = n
    result['nvar'] = nvar
    result['q'] = q
    result['r'] = r
    result['s'] = s

    for nv in range(nvar):
        if q > 0.0:
            cutpt = empquant(runs[:, nv], q)
            # work = (runs[:, nv] <= cutpt).astype(int)
            work = (runs[:, nv] <= cutpt).astype(int)
        else:
            q = 0.0
            i1 = np.where(runs[:, nv] == 0)[0]
            i2 = np.where(runs[:, nv] == 1)[0]
            if len(i1) + len(i2) != n:
                raise ValueError('raftery needs 0s and 1s in runs')
            work = runs[:, nv]
            q = np.sum(runs[:, nv]) / n

        kthin = (1)
        bic = 1.0
        epss = 0.001

        while bic > 0:
            tcnt, tmp = thin(work, n, kthin)
            _, bic = mctest(tmp, tcnt)
            kthin += 1

        kthin -= 1
        alpha, beta = mcest(tmp, tcnt)  # this lines can produce if chain is short [returns = nan for alpha]
        kmind = kthin
        _, bic = indtest(tmp, tcnt)

        while bic > 0:
            tcnt, tmp = thin(work, n, kmind)
            _, bic = indtest(tmp, tcnt)
            kmind += 1

        psum = alpha + beta

        # print(kmind,bic,alpha,beta,psum,kthin,tcnt)
        # Pre-caution which is not in MATLAB code
        # Check whether there is a typo somewhere
        if math.isnan(alpha) or math.isnan(beta):            
            nburn = nprec = kind = irl = nmin = nburn_plus_nprec = -(99)
        else:        
            tmp1 = np.log(psum * epss / max(alpha, beta)) / np.log(abs(1.0 - psum))
            nburn = int(np.floor((tmp1 + 1.0) * kthin))
            phi = ppnd((s + 1.0) / 2)
            tmp2 = (2.0 - psum) * alpha * beta * (phi ** 2) / (psum ** 3 * r ** 2)
            nprec = int(np.floor(tmp2 + 1.0) * kthin)
            nmin = int(np.floor(((1.0 - q) * q * phi ** 2) / r ** 2) + 1.0)
            irl = (nburn + nprec) / nmin
            kind = max(int(np.floor(irl + 1.0)), kmind)
            nburn_plus_nprec = nburn + nprec 

        result[nv] = {
            'nburn': nburn,
            'nprec': nprec,
            'kthin': kthin,
            'kind': kind,
            'irl': irl,
            'nmin': nmin,
            'n': nburn_plus_nprec
        }

    return result


def momentg(draws):
    ## ################################################################################## ##
    ## Computes Geweke's convergence diagnostics: numerical standard error and relative   ##
    ## numerical efficiencies NSE and RNE                                                 ##
    ##                                                                                    ##
    ## SYNOPSIS: result = momentg(draws)                                                  ##
    ##  where                                                                             ##
    ##   draws     [input]  REQUIRED: A Nxd matrix of draws                               ##
    ##   result    [outpt]  Dictionary containing convergence diagnostics each parameter  ##
    ## ################################################################################## ##
    
    # Get the number of draws and variables
    ndraw, nvar = draws.shape
    results = {}
    results['ndraw'] = ndraw
    results['meth'] = 'momentg'
    results['nvar'] = nvar
    
    NG = 20
    if ndraw < NG:
        raise ValueError('momentg: needs a larger number of ndraws')
    
    ntaper = [4, 8, 15]
    ns = ndraw // NG
    nuse = ns * NG

    for jf in range(nvar):  # Loop over all variables
        cnt = 0
        cn = np.zeros(NG)
        cd = np.zeros(NG)
        cdn = np.zeros(NG)
        cdd = np.zeros(NG)
        cnn = np.zeros(NG)
        cvar = np.zeros(NG)

        td = 0
        tn = 0
        tdd = 0
        tnn = 0
        tdn = 0
        tvar = 0

        # Form sufficiency statistics needed below
        for ig in range(NG):
            gd = 0
            gn = 0
            gdd = 0
            gdn = 0
            gnn = 0
            gvar = 0
            for is_ in range(ns):
                cnt += 1
                g = draws[cnt - 1, jf]
                ad = 1
                an = ad * g
                gd += ad
                gn += an
                gdn += ad * an
                gdd += ad * ad
                gnn += an * an
                gvar += an * g

            td += gd
            tn += gn
            tdn += gdn
            tdd += gdd
            tnn += gnn
            tvar += gvar

            cn[ig] = gn / ns
            cd[ig] = gd / ns
            cdn[ig] = gdn / ns
            cdd[ig] = gdd / ns
            cnn[ig] = gnn / ns
            cvar[ig] = gvar / ns

        eg = tn / td
        varg = tvar / td - eg**2
        sdg = -1
        if varg > 0:
            sdg = np.sqrt(varg)

        # Save posterior means and std deviations to results
        results[jf] = {
            'pmean': eg,
            'pstd': sdg
        }

        # Numerical standard error assuming no serial correlation
        varnum = (tnn - 2 * eg * tdn + tdd * eg**2) / (td**2)
        sdnum = -1
        if varnum > 0:
            sdnum = np.sqrt(varnum)

        results[jf]['nse'] = sdnum
        results[jf]['rne'] = varg / (nuse * varnum)

        # Get autocovariance of grouped means
        barn = tn / nuse
        bard = td / nuse
        for ig in range(NG):
            cn[ig] -= barn
            cd[ig] -= bard

        rnn = np.zeros(NG)
        rdd = np.zeros(NG)
        rnd = np.zeros(NG)
        rdn = np.zeros(NG)

        for lag in range(NG):
            ann = 0
            add = 0
            and_ = 0
            adn = 0
            for ig in range(lag + 1, NG):
                ann += cn[ig] * cn[ig - lag]
                add += cd[ig] * cd[ig - lag]
                and_ += cn[ig] * cd[ig - lag]
                adn += cd[ig] * cd[ig - lag]
            rnn[lag] = ann / NG
            rdd[lag] = add / NG
            rnd[lag] = and_ / NG
            rdn[lag] = adn / NG

        # Numerical standard error with tapered autocovariance functions
        for mm in range(3):
            m = ntaper[mm]
            am = m
            snn = rnn[0]
            sdd = rdd[0]
            snd = rnd[0]
            for lag in range(1, m - 1):
                att = 1 - lag / am
                snn += 2 * att * rnn[lag]
                sdd += 2 * att * rdd[lag]
                snd += att * (rnd[lag] + rnd[lag])

            varnum = ns * nuse * (snn - 2 * eg * snd + sdd * eg**2) / (td**2)
            sdnum = -1
            if varnum > 0:
                sdnum = np.sqrt(varnum)

            if mm == 0:
                results[jf]['nse1'] = sdnum
                results[jf]['rne1'] = varg / (nuse * varnum)
            elif mm == 1:
                results[jf]['nse2'] = sdnum
                results[jf]['rne2'] = varg / (nuse * varnum)
            elif mm == 2:
                results[jf]['nse3'] = sdnum
                results[jf]['rne3'] = varg / (nuse * varnum)

    return results


def apm(results1, results2):
    ## ################################################################################## ##
    ## Computes Geweke's chi-squared test for two sets of MCMC sample draws               ##
    ##                                                                                    ##
    ## SYNOPSIS: result = apm(results1,results2)                                          ##
    ##  where                                                                             ##
    ##   results1  [input]  REQUIRED: List of results from momentg (first MCMC set)       ##
    ##   results2  [input]  REQUIRED: List of results from momentg (second MCMC set)      ##
    ##   result    [outpt]  Dictionary containing Geweke's chi-squared test results       ##
    ## ################################################################################## ##
    
    if not isinstance(results1, dict) or not isinstance(results2, dict):
        raise ValueError('apm: requires a list from momentg as input')

    nvar  = results1['nvar']
    nvar2 = results2['nvar']
    
    if nvar != nvar2:
        raise ValueError('apm: structure arguments have different number of variables')

    ndraw1 = results1['ndraw']
    ndraw2 = results2['ndraw']
    
    result = {}
    result['p1'] = ndraw1 / (ndraw1 + ndraw2)
    result['p2'] = ndraw2 / (ndraw1 + ndraw2)
    result['ndraw'] = ndraw1 + ndraw2
    result['meth'] = 'apm'
    result['nvar'] = nvar

    ng = nvar
    nf = 2

    # Initialize arrays to hold the relevant data
    g = np.zeros((nf, nvar))
    sdnum1 = np.zeros((nf, nvar))
    sdnum2 = np.zeros((nf, nvar))
    sdnum3 = np.zeros((nf, nvar))
    sdnum4 = np.zeros((nf, nvar))
    # Pull out information from the results
    for i in range(nvar):
        g[0, i] = results1[i]['pmean']
        sdnum1[0, i] = results1[i]['nse']
        sdnum2[0, i] = results1[i]['nse1']
        sdnum3[0, i] = results1[i]['nse2']
        sdnum4[0, i] = results1[i]['nse3']
        
        g[1, i] = results2[i]['pmean']
        sdnum1[1, i] = results2[i]['nse']
        sdnum2[1, i] = results2[i]['nse1']
        sdnum3[1, i] = results2[i]['nse2']
        sdnum4[1, i] = results2[i]['nse3']
    
    # Loop over variables
    for i in range(nvar):
        for k in range(4):
            eg = 0
            nse = 0
            wtsum = 0

            if k == 0:
                sdnum = sdnum1
            elif k == 1:
                sdnum = sdnum2
            elif k == 2:
                sdnum = sdnum3
            elif k == 3:
                sdnum = sdnum4

            gvar = np.zeros((nf - 1, nf - 1))

            # Compute the weighted sum
            for j in range(nf):
                eg += g[j, i] / (sdnum[j, i] ** 2)
                wtsum += 1 / (sdnum[j, i] ** 2)

            eg /= wtsum
            nse = 1 / np.sqrt(wtsum)

            for j in range(nf - 2):
                gvar[j, j] = (sdnum[j, i] ** 2) + (sdnum[j + 1, i] ** 2)
                gvar[j, j + 1] = -(sdnum[j + 1, i] ** 2)
                gvar[j + 1, j] = gvar[j, j + 1]

            gvar[nf - 2, nf - 2] = (sdnum[nf - 2, i] ** 2) + (sdnum[nf - 1, i] ** 2)
            ginv = np.linalg.inv(gvar)
            g1 = g[:nf - 1, i]
            g2 = g[1:nf, i]
            cstat = np.dot((g2 - g1).T, np.dot(ginv, (g2 - g1)))
            df = nf - 1
            p = 1 - stats.chi2.cdf(cstat, df)           
            # Save results
            result[i] = result.get(i, {})
            result[i][f'pmean_{k+1}'] = eg
            result[i][f'nse_{k+1}'] = nse
            result[i][f'prob_{k+1}'] = p

    return result


def empquant(runs, q):
    ## ################################################################################## ##
    ## Empirical quantile function. Calculates the q-th quantile from the empirical       ##
    ## distribution of "runs"                                                             ##
    ##                                                                                    ##
    ## SYNOPSIS: y = empquant(runs, q)                                                    ##
    ##  where                                                                             ##
    ##   runs      [input]  REQUIRED: Nx1 vector of values from which to compute quantile ##
    ##   q         [input]  REQUIRED: Quantile of interest (q ∈ [0,1])                    ##
    ##   y         [outpt]  The q-th quantile of the "runs" array                         ##
    ## ################################################################################## ##

    # Sort the input runs
    work = np.sort(runs)
    # Calculate the order index for the quantile
    n = len(runs)
    order = (n - 1) * q + 1.0    
    # Fractional part of the order
    fract = order % 1.0    
    # Low and high indices based on the order
    # low = max(int(np.floor(order)), 0)  # Index should be within bounds (0 to n-1)
    # high = min(low + 1, n - 1)          # Ensure the high index is within bounds    
    # Compute the quantile using linear interpolation
    # y = (1.0 - fract) * work[low] + fract * work[high]
    low = max(int(np.floor(order)), 1)  # Index should be within bounds (0 to n-1)
    high = min(low + 1, n)              # Ensure the high index is within bounds    
    # Compute the quantile using linear interpolation
    y = (1.0 - fract) * work[low - 1] + fract * work[high - 1]
    
    return y


def thin(run, n, kthin):
    ## ################################################################################## ##
    ## Thinning function: returns every k-th element from "run" array                     ##
    ##                                                                                    ##
    ## SYNOPSIS: x, y = thin(run, n, kthin)                                               ##
    ##  where                                                                             ##
    ##   run       [input]  REQUIRED: Nx1 vector or list of values                        ##
    ##   n         [input]  REQUIRED: Length of the array "run"                           ##
    ##   kthin     [input]  REQUIRED: Thinning interval                                   ##
    ##   x         [outpt]  Length of the thinned array                                   ##
    ##   y         [outpt]  Thinned array                                                 ##
    ## ################################################################################## ##

    # Generate the indices for thinning
    ind = np.arange(0, n, kthin)    
    # Thinned result
    y = run[ind]    
    # Length of the thinned array
    x = len(y)
    
    return x, y


def mctest(d, n):
    ## ################################################################################## ##
    ## Perform the Geweke's G^2 test to evaluate the accuracy of MCMC simulations         ##
    ##                                                                                    ##
    ## SYNOPSIS: g2, BIC = mctest(d, n)                                                   ##
    ##  where                                                                             ##
    ##   d         [input]  REQUIRED: Nx2 matrix of MCMC draws                            ##
    ##   n         [input]  REQUIRED: Number N of draws                                   ##
    ##   g2        [outpt]  Geweke G-squared statistic (= likelihood ratio)               ##
    ##   BIC       [outpt]  Bayes Information Criterion                                   ##
    ## ################################################################################## ##

    m1 = np.zeros((2, 2))  # Initialize 2x2 matrix
    m2 = np.zeros((2, 2))  # Initialize 2x2 matrix
    g2 = 0.0   
    # Count the states
    for i in range(2, n):  # Loop through draws starting from the third one
        i1 = int(d[i-2])
        i2 = int(d[i-1])
        i3 = int(d[i]) + 1
        if i3 == 1:
            m1[i1, i2] += 1
        if i3 == 2:
            m2[i1, i2] += 1

    # Compute G^2 and BIC
    for i1 in range(2):         # Loop over the first dimension
        for i2 in range(2):     # Loop over the second dimension
            for i3 in range(2): # Loop over the third dimension
                if i3 == 0:
                    if m1[i1, i2] != 0:
                        t1 = m1[i1, i2] + m2[i1, i2]
                        t2 = m1[0, i2] + m1[1, i2]
                        t3 = m1[0, i2] + m2[0, i2]
                        t4 = m1[1, i2] + m2[1, i2]
                        fitted = (t1 * t2) / (t3 + t4)
                        focus = m1[i1, i2]
                        g2 += np.log(focus / fitted) * focus
                if i3 == 1:
                    if m2[i1, i2] != 0:
                        t1 = m1[i1, i2] + m2[i1, i2]
                        t2 = m2[0, i2] + m2[1, i2]
                        t3 = m1[0, i2] + m2[0, i2]
                        t4 = m1[1, i2] + m2[1, i2]
                        fitted = (t1 * t2) / (t3 + t4)
                        focus = m2[i1, i2]
                        g2 += np.log(focus / fitted) * focus
    
    g2 *= 2.0
    BIC = g2 - np.log(n - 2.0) * 2.0
    
    return g2, BIC


def mcest(d, n):
    ## ################################################################################## ##
    ## Compute the alpha and beta values for the Markov Chain Monte Carlo (MCMC) draws    ##
    ##                                                                                    ##
    ## SYNOPSIS: alpha, beta = mcest(d, n)                                                ##
    ##  where                                                                             ##
    ##   d         [input]  REQUIRED: Nx2 matrix of MCMC draws                            ##
    ##   n         [input]  REQUIRED: Number N of draws                                   ##
    ##   alpha     [outpt]  Transition probability from state 1 to state 2                ##
    ##   beta      [outpt]  Transition probability from state 2 to state 1                ##
    ## ################################################################################## ##

    t = np.zeros((2, 2))                    # Initialize a 2x2 transition matrix
    # Populate the transition matrix with counts
    for i1 in range(1, n):                  # Loop through the MCMC draws
        t[int(d[i1-1]), int(d[i1])] += 1

    # Calculate alpha and beta
    alpha = t[0, 1] / (t[0, 0] + t[0, 1])   # Transition probability from state 1 to state 2
    beta = t[1, 0] / (t[1, 0] + t[1, 1])    # Transition probability from state 2 to state 1
    
    return alpha, beta


def indtest(d, n):
    ## ################################################################################## ##
    ## Perform independence test for Markov Chain Monte Carlo (MCMC) draws                ##
    ##                                                                                    ##
    ## SYNOPSIS: g2, BIC = indtest(d, n)                                                  ##
    ##  where                                                                             ##
    ##   d         [input]  REQUIRED: Nx1 vector of MCMC draws                            ##
    ##   n         [input]  REQUIRED: Number N of draws                                   ##
    ##   g2        [outpt]  Geweke G-squared statistic (= likelihood ratio)               ##
    ##   BIC       [outpt]  Bayes Information Criterion                                   ##
    ## ################################################################################## ##

    t = np.zeros((2, 2))  # Initialize a 2x2 transition matrix

    # Populate the transition matrix with counts
    for i1 in range(1, n):  # Loop through the MCMC draws
        t[int(d[i1-1]), int(d[i1])] += 1

    dcm1 = n - 1.0  # Adjustment for degrees of freedom
    g2 = 0.0  # G-squared statistic initialization
    # Compute the G-squared statistic (likelihood ratio)
    for i1 in range(2):
        for i2 in range(2):
            if t[i1, i2] != 0:
                t1 = t[i1, 0] + t[i1, 1]                # Row sum
                t2 = t[0, i2] + t[1, i2]                # Column sum
                fitted = (t1 * t2) / dcm1               # Fitted value
                focus = t[i1, i2]                       # Observed count
                g2 += np.log(focus / fitted) * focus    # Update G-squared statistic

    g2 *= 2.0                                           # Multiply by 2
    BIC = g2 - np.log(dcm1)                             # Bayesian Information Criterion

    return g2, BIC


def ppnd(p):
    ## ################################################################################## ##
    ## Computes quantile of cumulative standard normal distribution for given probability ##
    ##                                                                                    ##
    ## SYNOPSIS: y = ppnd(p)                                                              ##
    ##  where                                                                             ##
    ##   p         [input]  REQUIRED: probability, p ∈ [0,1]                              ##
    ##   y         [outpt]  Quantile of standard normal distribution                      ##
    ## ################################################################################## ##

    # Constants
    split1 = 0.425
    split2 = 5.0
    const1 = 0.180625
    const2 = 1.6
    a0 = 3.3871327179e+00
    a1 = 5.0434271938e+01
    a2 = 1.5929113202e+02
    a3 = 5.9109374720e+01
    b1 = 1.7895169469e+01
    b2 = 7.8757757664e+01
    b3 = 6.7187563600e+01
    c0 = 1.4234372777e+00
    c1 = 2.7568153900e+00
    c2 = 1.3067284816e+00
    c3 = 1.7023821103e-01
    d1 = 7.3700164250e-01
    d2 = 1.2021132975e-01
    e0 = 6.6579051150e+00
    e1 = 3.0812263860e+00
    e2 = 4.2868294337e-01
    e3 = 1.7337203997e-02
    f1 = 2.4197894225e-01
    f2 = 1.2258202635e-02

    # Calculate q
    q = p - 0.5
    
    if abs(q) <= split1:
        r = const1 - q * q
        y = q * (((a3 * r + a2) * r + a1) * r + a0) / (((b3 * r + b2) * r + b1) * r + 1.0)
        return y
    elif q < 0.0:
        r = p
    else:
        r = 1 - p

    if r <= 0.0:
        return 0.0
    
    r = math.sqrt(-math.log(r))
    
    if r <= split2:
        r = r - const2
        y = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + 1.0)
    else:
        r = r - split2
        y = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + 1.0)
    
    if q < 0.0:
        y = -y
    
    return y


################################## End Single chain diagnostics ####################################

########################3### Scoring rules for posterior distribution ##############################

def log_score(fcst, obs, mu_F = None, sigma_F = None):
    ## ################################################################################## ##
    ## Logarithmic score of the observations for an ensemble forecast at different times  ##
    ##                                                                                    ##
    ## SYNOPSIS: LS, LS_value, num_zero = log_score(fcst, obs, mu_F, sigma_F)             ##
    ##  where                                                                             ##
    ##   fcst      [input]  REQUIRED: nxm matrix of m emsemble forecasts at n times       ##
    ##   obs       [input]  REQUIRED: nx1 vector of observations at n times               ##
    ##   mu_F      [input]  OPTIONAL: nx1 vector with mean of ensemble forecasts          ##
    ##   sigma_F   [input]  OPTIONAL: nx1 vector with standard dev. of ensemble forecasts ##
    ##   LS        [outpt]  Mean of logarithmic score at n times                          ##
    ##   LS_value  [outpt]  nx1 vector with logarithmic score at each time                ##
    ##   num_zero  [outpt]  Number of logarithmic scores > -inf                           ##
    ## ################################################################################## ##
    
    # Determine if G = 0 (kernel density) or G = 1 (normal PDF)
    G = 1 if mu_F is not None and sigma_F is not None else 0
    # Get the shape of the forecast matrix (n x m)
    n, m = fcst.shape
    if obs.shape[0] != n:
        raise ValueError("log_score: The length of the observation vector does not match the number of rows of forecast matrix.")
    
    if obs.ndim != 1:
        raise ValueError("log_score: The observation vector should have one column only.")
    
    if G == 0:
        # Initialize results
        LS = np.full(n, np.nan)
        num_zero = 0
        # Loop over each entry of the measurement vector
        for t in range(n):
            # Get the predictive pdf of forecast using kernel density estimation
            kde = gaussian_kde(fcst[t, :])      # Kernel density estimation for the t-th observation
            py = kde.evaluate([obs[t]])         # Evaluate the pdf at the observed value
            if py < np.finfo(float).eps:        # Check for values close to zero
                py = np.finfo(float).eps
                num_zero += 1
            # Compute the logarithmic score
            LS[t] = np.log(py)
    elif G == 1:
        # Use normal PDF with provided mean and standard deviation
        # py = norm.pdf(obs, mu_F, sigma_F)
        # Check for values close to zero
        # ii = py < np.finfo(float).eps
        # py[ii] = np.finfo(float).eps
        # num_zero = np.sum(ii)
        py = np.maximum(norm.pdf(obs, mu_F, sigma_F), np.finfo(float).eps)
        num_zero = sum(py == np.finfo(float).eps)
        # Compute the logarithmic score
        LS = np.log(py)
    # Return the mean of non-zero logarithmic scores
    mLS = np.nanmean(LS[LS > -np.inf])          # Mean ignoring -inf values (invalid scores)
    
    return mLS, LS, num_zero


def CRP_score(fcst, obs, type_cdf='ecdf', calc_method = 4):
    ## ################################################################################## ##
    ## Continuous ranked probability score of observations for an ensemble forecast at    ##
    ## different times. The CRPS is a quadratic measure of the difference between the     ##
    ## forecast CDF and the empirical CDF of the observations.                            ##
    ##                                                                                    ##
    ## SYNOPSIS: mCRPS, CRPS, num_zero = CRP_score(fcst, obs, type_cdf, method)           ##
    ##  where                                                                             ##
    ##   fcst      [input]  REQUIRED: nxm matrix of m emsemble forecasts at n times       ##
    ##   obs       [input]  REQUIRED: nx1 vector of observations at n times               ##
    ##   type_cdf  [input]  OPTIONAL: Type of empirical cdf (def: 'ecdf')                 ##
    ##   method    [input]  OPTIONAL: Calculation/computation method (def: 4)             ##
    ##   mCRPS     [outpt]  Mean of continuous ranked probability score                   ##
    ##   CRPS      [outpt]  nx1 vector with CRPS values at each time                      ##
    ##   num_nan   [outpt]  Number of CRP scores equal to NaN                             ##
    ## ################################################################################## ##

    # Validate input type_cdf and set p1, p2 based on the provided type
    type_cdf = type_cdf.lower() if isinstance(type_cdf, str) else type_cdf
    if isinstance(type_cdf, str):
        if type_cdf == 'ecdf':
            p1, p2 = 0, 0
        elif type_cdf == 'weib':
            p1, p2 = 0, 1
        elif type_cdf == 'med':
            p1, p2 = 0.3175, 0.365
        elif type_cdf == 'apl':
            p1, p2 = 0.35, 0
        elif type_cdf == 'blom':
            p1, p2 = 0.375, 0.25
        elif type_cdf == 'cunn':
            p1, p2 = 0.4, 0.2
        elif type_cdf == 'grin':
            p1, p2 = 0.44, 0.12
        elif type_cdf == 'hazen':
            p1, p2 = 0.5, 0
        elif type_cdf == 'bern':
            p1, p2 = 0.3, 0.4
        else:
            raise ValueError(f"CRP_score: Unknown empirical CDF: {type_cdf}")
    elif isinstance(type_cdf, tuple) and len(type_cdf) == 2:
        p1, p2 = type_cdf
        if not (0 <= p1 <= 1 and 0 <= p2 <= 1):
            raise ValueError("CRP_score: Parameters p1 and p2 must be between 0 and 1.")
    else:
        raise ValueError("CRP_score: Invalid type_cdf input.")

    # Validate the forecast and observation shapes
    n, m = fcst.shape
    if obs.shape[0] != n or obs.ndim != 1:
        raise ValueError("CRP_score: The length of the observation vector does not match the number of rows of the forecast matrix.")

    # Initialize CRPS values
    CRPS = np.full(n, np.nan)
    ys = np.sort(fcst, axis = 1)  # Sort the forecast ensemble along each row

    # Compute CRPS based on the selected calc_method
    if calc_method == 1:
        p = ((np.arange(1, m + 1) - p1) / (m + p2))
        r1 = p**2
        r2 = (1.0 - p)**2
        for t in range(n):
            missingFcst = np.any(np.isnan(fcst[t, :]))
            missingObs = np.isnan(obs[t])
            if not missingFcst and not missingObs:
                ind = np.searchsorted(fcst[t, :], obs[t])
                crpsLeft = 0
                if ind > 0:
                    fcstLeft = fcst[t, :ind]
                    dxLeft = np.diff(fcstLeft)
                    pLeft = r1[:ind - 1]
                    crpsLeft = np.dot(pLeft, dxLeft)
                if obs[t] < fcst[t, -1]:
                    fcstRight = fcst[t, ind:]
                    dxRight = np.diff(fcstRight)
                    if dxRight.size > 0:
                        pRight = r2[ind:]
                        crpsRight = np.dot(pRight, dxRight)
                    else:
                        crpsRight = 0
                    crpsCentreLeft = r1[ind] * (obs[t] - fcst[t, ind])
                    crpsCentreRight = r2[ind] * (fcst[t, ind + 1] - obs[t])
                    CRPS[t] = crpsLeft + crpsRight + crpsCentreLeft + crpsCentreRight
                else:
                    crps_right_outside = 1.0**2 * (obs[t] - fcst[t, -1])
                    CRPS[t] = crps_right_outside + crpsLeft
            elif not missingObs:
                CRPS[t] = 1.0**2 * (obs[t] - fcst[t, 0])
    elif calc_method == 2:
        p = ((np.arange(0, m + 1) - p1) / (m + p2))
        r1 = p**2
        r2 = (1 - p)**2
        for t in range(n):
            CRPS[t] = np.dot(r1, np.diff(ys[t, :])) + np.dot(r2, np.diff(ys[t, ::-1]))
    elif calc_method == 3:
        for t in range(n):
            CRPS[t] = 2 * np.sum((ys[t, :m] - obs[t]) * (m * (obs[t] < ys[t, :m]) - np.arange(1, m + 1) + 1 / 2)) / m**2
    elif calc_method == 4:
        CRPS = 2 * np.sum((ys[:, :m] - obs[:, None]) * (m * (obs[:, None] < ys[:, :m]) - np.arange(1, m + 1) + 1 / 2), axis = 1) / m**2
    else:
        raise ValueError(f"CRP_score: Unknown calc_method: {calc_method}")

    CRPS = -CRPS                        # Reversing the sign to get a positively oriented score
    mCRPS = np.nanmean(CRPS)            # Mean of non-missing CRPS values
    num_nan = np.sum(np.isnan(CRPS))    # Number of missing CRPS values

    return mCRPS, CRPS, num_nan


def spherical_score(fcst, obs, zeta = 2, calc_method = 2, N = 1000):
    ## ################################################################################## ##
    ## Spherical score of the observations for an ensemble forecast at different times.   ##
    ##                                                                                    ##
    ## SYNOPSIS: SS, SS_value, num_inf = spherical_score(fcst, obs, zeta, method, N)      ##
    ##  where                                                                             ##
    ##   fcst      [input]  REQUIRED: nxm matrix of m emsemble forecasts at n times       ##
    ##   obs       [input]  REQUIRED: nx1 vector of observations at n times               ##
    ##   zeta      [input]  OPTIONAL: Exponential power (def: 2 = spherical score)        ##
    ##   method    [input]  OPTIONAL: Calculation/computation method (def: 2)             ##
    ##   N         [input]  OPTIONAL: Integer (def: N = 1000)                             ##
    ##   mSS       [outpt]  Mean of spherical score                                       ##
    ##   SS        [outpt]  nx1 vector with spherical scores at each time                 ##
    ##   num_inf   [outpt]  Number of spherical scores equal to infinity                  ##
    ## ################################################################################## ##

    # Ensure input dimensions are valid
    n, m = fcst.shape
    if obs.shape[0] != n:
        raise ValueError('spherical_score: The length of the observation vector does not match'
                         'number of rows of forecast matrix')
    if obs.ndim != 1:
        raise ValueError('spherical_score: The observation vector should have one column only')
    if zeta < 1:
        raise ValueError('spherical_score: The pseudospherical score is not defined for zeta < 1')

    # Determine range, minimum, and maximum of forecast
    r = np.ptp(fcst, axis = 1)
    dr = r / 5
    xf_min = np.min(fcst, axis = 1) - dr
    xf_max = np.max(fcst, axis = 1) + dr

    if calc_method == 0:
        n2 = 100
        Xf = np.full((n, n2), np.nan)
    elif calc_method in [1, 2]:
        n2 = 2**np.ceil(np.log2(N)).astype(int)
        Xf = np.full((n, n2), np.nan)
    
    for t in range(n):
        Xf[t, :] = np.linspace(xf_min[t], xf_max[t], n2)

    SS = np.full(n, np.nan)
    
    # Loop over each observation
    for t in range(n):
        # Get predictive pdf of forecast; evaluate pdf at observed and simulated values
        if calc_method == 0:
            # Evaluate pdf of forecast and return 100 points
            kde = gaussian_kde(fcst[t, :])
            Xf[t, :n2] = np.linspace(xf_min[t], xf_max[t], n2)
            f = kde(Xf[t, :n2])
            f_obs = kde(obs[t])
        elif calc_method == 1:
            # Evaluate pdf at obs and xf at the same time
            kde = gaussian_kde(fcst[t, :])
            f = kde([obs[t], *Xf[t, :n2]])
            f = f / np.sum(f[1:])  # Normalize the PDF
            f_obs = f[0]
            f = f[1:]
        elif calc_method == 2:
            # External kernel density estimator calc_method
            kde = gaussian_kde(fcst[t, :])
            f = kde(Xf[t, :n2])
            f_obs = np.interp(obs[t], Xf[t, :n2], f)

        # Now integrate the distribution
        I2 = np.trapezoid(f**zeta, Xf[t, :n2])
        
        # Compute spherical score: positive orientation [larger = better]
        SS[t] = (f_obs**(zeta - 1)) / (I2**((zeta - 1) / zeta))
    
    ii = SS > -np.inf               # Get index of spherical scores > -infinity
    mSS = np.nanmean(SS[ii])        # Return the mean of the non-missing spherical scores
    num_inf = n - np.sum(ii)        # Number of spherical scores of -infinity

    return mSS, SS, num_inf


def dawid_sebas_score(fcst, obs, mu_F = None, sigma_F = None):
    ## ################################################################################## ##
    ## Dawid Sebastiani score of the observations for an ensemble forecast at different   ##
    ## times.                                                                             ##
    ##                                                                                    ##
    ## SYNOPSIS: DS, DS_value, num_inf = dawid_sebas_score(fcst, obs, mu_F, sigma_F)      ##
    ##  where                                                                             ##
    ##   fcst      [input]  REQUIRED: nxm matrix of m emsemble forecasts at n times       ##
    ##   obs       [input]  REQUIRED: nx1 vector of observations at n times               ##
    ##   mu_F      [input]  OPTIONAL: nx1 vector of mean ensemble distribution forecast F ##
    ##   sigma_F   [input]  OPTIONAL: nx1 vector of standard dev. ensemble forecasts      ##
    ##   mDSS      [outpt]  Mean of non-missing Dawid-Sebastiani scores                   ##
    ##   DSS       [outpt]  nx1 vector of Dawid-Sebastiani scores at each time            ##
    ##   num_nan   [outpt]  Number of NaN values of Dawid-Sebastiani score                ##
    ##   mu_F      [outpt]  nx1 vector of mean ensemble distribution forecast F           ##
    ##   sigma_F   [outpt]  nx1 vector of standard dev. ensemble forecasts                ##
    ## ################################################################################## ##

    # Check input dimensions
    n, m = fcst.shape       # Number of measurement times and ensemble members
    if obs.shape[0] != n:
        raise ValueError('dawid_sebas_score: Length of observation vector does not match number '
                         'of rows of forecast matrix')
    if obs.ndim != 1:
        raise ValueError('dawid_sebas_score: Observation vector should have one column')

    if mu_F == None:
        mu_F = np.mean(fcst, axis = 1)      # Mean along the second dimension (= columns)
    if sigma_F == None:
        sigma_F = np.std(fcst, axis = 1)    # Std along the second dimension (= columns)

    # Compute Dawid-Sebastiani score
    DSS = -2 * np.log(sigma_F) - ((obs - mu_F) / sigma_F) ** 2   
    # Compute the mean of non-missing Dawid-Sebastiani scores
    mDSS = np.nanmean(DSS)
    # Count the number of nan values
    num_nan = np.sum(np.isnan(DSS))    
    
    return mDSS, DSS, num_nan, mu_F, sigma_F


def interval_score(LU, obs, p_alfa = 0.05):
    ## ################################################################################## ##
    ## Interval score of the observations for an ensemble forecast at different times     ##
    ##                                                                                    ##
    ## SYNOPSIS: mIS, IS = interval_score(LU, obs, p_alfa)                                ##
    ##  where                                                                             ##
    ##   LU        [input]  REQUIRED: nx2 matrix of lower and upper prediction quantiles  ##
    ##   obs       [input]  REQUIRED: nx1 vector of observations at n times               ##
    ##   p_alfa    [input]  OPTIONAL: significance level (default = 0.05)                 ##
    ##   mIS       [outpt]  Mean of non-missing Dawid-Sebastiani scores                   ##
    ##   IS        [outpt]  nx1 vector of Dawid-Sebastiani scores at each time            ##
    ## ################################################################################## ##
    
    # Check input dimensions
    n, r = LU.shape  # n: number of measurement times, r: should be 2 (lower and upper quantiles)
    if r != 2:
        raise ValueError('The LU array should have two columns with lower and upper quantiles')
    if obs.shape[0] != n:
        raise ValueError('Length of observation vector does not match number of rows of forecast matrix')
    if obs.ndim != 1:
        raise ValueError('The observation vector should have one column only')
    
    # Extract the lower and upper quantiles
    L = LU[:, 0]  # Lower quantiles
    U = LU[:, 1]  # Upper quantiles   
    # Compute the interval score for all entries
    IS = (U - L) + 2/p_alfa * (L - obs) * (obs < L) + 2/p_alfa * (obs - U) * (obs > U)    
    # Return the mean of all interval scores (positive orientation, larger = better)
    mIS = - np.mean(IS)
    
    return mIS, IS

########################## End Scoring rules for posterior distribution ############################

################################ Performance metrics for posterior #################################
def rlbl(p_val):
    ## ################################################################################## ##
    ## Computes the time-averaged reliability of the forecast distribution of ensemble    ##
    ##                                                                                    ##
    ## SYNOPSIS: RLBL, eCDF, uCDF = rlbl(p_val)                                           ##
    ##  where                                                                             ##
    ##   p_val     [input]  REQUIRED: nx1 vector of p-values of ensemble forecasts        ##
    ##   RLBL      [outpt]  Reliability of time series of ensemble forecasts              ##
    ##   eCDF      [outpt]  nx1 vector with empirical CDF of sorted p-values              ##
    ##   uCDF      [outpt]  nx1 vector with CDF of uniform distribution                   ##
    ## ################################################################################## ##
    
    # Number of p-values (number of observations)
    n = len(p_val)    
    # Empirical CDF (sort p-values in ascending order)
    eCDF = np.sort(p_val)    
    # CDF of uniform distribution
    uCDF = np.arange(1, n + 1) / n    
    # Compute reliability
    RLBL = 1 - (2 / n) * np.sum(np.abs(uCDF - eCDF))
    
    return RLBL, eCDF, uCDF


def p_values(fcst, obs, calc_method = 2):
    ## ################################################################################## ##
    ## Computes the p-values of the observations given the ensemble forecasts at n times  ##
    ##                                                                                    ##
    ## SYNOPSIS: p_val, num_nan = p_values(fcst, obs, method)                             ##
    ##  where                                                                             ##
    ##   fcst      [input]  REQUIRED: nxm matrix of m emsemble forecasts at n times       ##
    ##   obs       [input]  REQUIRED: nx1 vector of observations at n times               ##
    ##   method    [input]  OPTIONAL: Calculation/computation method (def: 2)             ##
    ##   p_val     [outpt]  nx1 vector of p-values of ensemble forecasts                  ##
    ##   num_nan   [outpt]  Number of missing values (NaN) of the p-values                ##
    ## ################################################################################## ##

    # Get the size of the forecast matrix
    n, m = fcst.shape
    
    # Check for consistency of dimensions
    if obs.shape[0] != n:
        raise ValueError('p_values: The length of the observation vector does not match number '
                         'of rows of forecast matrix.')
    if obs.ndim != 1:
        raise ValueError('p_values: The observation vector should have one column only.')
    
    # Calculate p-values using the selected calc_method
    if calc_method == 1:     # Loop over time
        p_val = np.full(n, np.nan)
        for t in range(n):
            p_val[t] = np.sum(fcst[t, :] < obs[t]) / m
    elif calc_method == 2:   # Vectorized implementation
        p_val = np.sum(fcst < obs[:, None], axis=1) / m
    else:
        raise ValueError('Invalid calc_method. Use 1 for loop or 2 for direct (vectorized).')
    
    # Count NaN values
    num_nan = np.sum(np.isnan(p_val))
    
    return p_val, num_nan


############################## End Performance metrics for posterior ###############################

# The following items need to be looked at further: 
# 1. Small numerical differences between GL+/UL functions in MATLAB and Python [do not seem consequential]
# 2. Check mt-dream: [Done] (but under/overflow theoretically possible?)
# 3. Kalman proposal: [Done] (but check all possible configurations) 

## OLD
# Determine how many variables [objects] were received from function handle
# nvar = max(1,sum(isinstance(item, np.ndarray) for sublist in worker_result for item in sublist))
# print("nvar",nvar,"N",N,"worker_result",worker_result)
# if N > DREAMPar['CPU']:     ## must divide by number of X's evaluated by each worker
#     nvar = nvar // (task_ranges[i][1] - task_ranges[i][0])
# nobj = sum(isinstance(item, object) for sublist in worker_result for item in sublist)
# nvar2 = nobj // (n * (task_ranges[i][1] - task_ranges[i][0]))
# i = i + 1
# print("var",nvar,"nvar2",nvar2)
