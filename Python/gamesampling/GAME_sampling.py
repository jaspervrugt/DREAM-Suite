# ####################################################################### #
#                                                                         #
#             GGGGGGGGG     AAAA     MMM      MMM  EEEEEEEE               #
#             GGGGGGGGG    AAAAAA    MMM      MMM  EEEEEEEE               #
#             GGG   GGG   AAA  AAA   MMM      MMM  EEE                    #
#             GGG   GGG  AAA    AAA  MMMM    MMMM  EEE                    #
#             GGGGGGGGG  AAA    AAA  MMMMM  MMMMM  EEEEE                  #
#              GGGGGGGG  AAAAAAAAAA  MMM MMMM MMM  EEEEE                  #
#                   GGG  AAAAAAAAAA  MMM  MM  MMM  EEE                    #
#                   GGG  AAA    AAA  MMM      MMM  EEE                    #
#              GGGGGGGG  AAA    AAA  MMM      MMM  EEEEEEEE               #
#             GGGGGGGGG  AAA    AAA  MMM      MMM  EEEEEEEE               #
#                                                                         #
# SSSSSSS     A     MMM     MMM PPPPPP   LLL      III NNN    NNN GGGGGGGG #
# SSSSSS     AAA    MMM     MMM PPPPPPP  LLL      III NNN    NNN GGGGGGGG #
# SSS       AAAAA   MMM     MMM PPP  PPP LLL      III NNNN   NNN GGG  GGG #
# SSS      AAA AAA  MMMM   MMMM PPP  PPP LLL      III NNNNN  NNN GGG  GGG #
# SSSSSSS AAA   AAA MMMMM MMMMM PPPPPPP  LLL      III NNN NN NNN GGGGGGGG #
# SSSSSSS AAAAAAAAA MMM MMM MMM PPPPPP   LLL      III NNN  NNNNN  GGGGGGG #
#     SSS AAAAAAAAA MMM  M  MMM PPP      LLL      III NNN   NNNN      GGG #
#     SSS AAA   AAA MMM     MMM PPP      LLL      III NNN    NNN      GGG #
#  SSSSSS AAA   AAA MMM     MMM PPP      LLLLLLLL III NNN    NNN  GGGGGGG #
# SSSSSSS AAA   AAA MMM     MMM PPP      LLLLLLLL III NNN    NNN GGGGGGGG #
#                                                                         #
# ####################################################################### #
#                                                                         #
# This function computes the marginal likelihood from a collection of     #
# samples, X, of the target distribution derived from (e)DREAM Package.   #
# Type "help GAME_fit_mixture" for more information about the procedure   #
#                                                                         #
# MAIN IDEA:                                                              #
#  Calculate the marginal likelihood as weighted mean of the ratio of the #
#  samples' target density, p(x_i), and the importance density, q(x_i).   #
#  This works because we know that the importance distribution, a mixture #
#  of normal distributions, integrates to 1. As a result, we yield that   #
#  Z = 1/N * sum_{i=1}^{N} (p(x_i)/q(x_i)) or as we solve herein in       #
#  logarithmic form (numerically more stable)                             #
#                                                                         #
# SYNOPSIS:                                                               #
#  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name);            #
#  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         #
#      GAMEoptions);                                                      #
#  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         #
#      GAMEoptions,Par_info);                                             #
#  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         #
#      GAMEoptions,Par_info,Meas_info);                                   #
#  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         #
#      GAMEoptions,Par_info,Meas_info,options);                           #
#  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         #
#      GAMEoptions,Par_info,Meas_info,options,plugin);                    #
# WHERE                                                                   #
#  X           [input] Rx(d+2) matrix posterior samples DREAM             #
#  method      [input] string (name) marginal likelihood estimator        #
#   = 'ris'        reciprocal importance sampling                         #
#   = 'is'         importance sampling                                    #
#   = 'ob'         optimal bridge sampling                                #
#   = 'gb'         geometric bridge sampling                              #
#  DREAMPar    [input] Structure with algorithmic variables               #
#   .d             Dimensionality (# variables) target distribution       #
#   .N             # of Markov chains                                     #
#   .T             # of generations (= # samples of each Markov chain)    #
#   .lik           Choice of likelihood function                          #
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
#                   → DREAM_ZS and DREAM_DZS                              #
#   .m0            Initial size of external archive, Z   DEF: 10*d        #
#                   → DREAM_ZS and DREAM_DZS                              #
#   .k             Growth rate of external archive       DEF: 10          #
#                   → DREAM_ZS and DREAM_DZS                              #
#   .M             # samples archive Z for Kalman jump   DEF: 20          #
#                   → DREAM_KZS                                           #
#   .mt            Number of multi-try proposals         DEF: 5           #
#                   → MTDREAM_ZS                                          #
#  Func_name   [input] Function (string) returns (log)lik or sim values   #
#   → Func_name must return a likelihood if DREAMPar.lik = 1              #
#   → Func_name must return a log-likelihood if DREAMPar.lik = 2          #
#   → Func_name returns vector of n simulated values: DREAMPar.lik > 2    #
#  GAMEoptions [input] (optional) GAME structure with additional options  #
#   = struct('metric','BIC','K',5,'N',1e4,'M',1,'steps',10);              #
#   .metric        metric for optimal mixture selection                   #
#    = 'bic'       Bayesian information criterion        DEFault          #
#    = 'var'       Variance reduction                                     #
#   .K             maximum # components mixture dist.    DEF: 5           #
#   .N             # importance samples to be used       DEF: 10,000      #
#   .M             # times we repeat each estimator      DEF: 1           #
#   .steps         # steps with gb or ob                 DEF: 10          #
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
#  Meas_info   [input] Structure with measurement information (fitting)   #
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
#  plugin      [input] 2nd input argument Func_name. Class set by user    #
#  Z           [outpt] marginal likelihood                                #
#                  = integral of posterior pdf                            #
#                  = integral of p(x|Y)                                   #
#                  = integral of p(x)L(x|Y)                               #
#  logZ        [outpt] logarithmic value of marginal likelihood           #
#  gmix        [outpt] Structure trained normal mixture importance dist.  #
#  .n              # samples (= R) of X                                   #
#  .d              # dimensions (= parameters) of X                       #
#  .k              # components of mixture distribution                   #
#  .p              # parameters of mixture distribution                   #
#  .w              maximum likelihood weights of mixture components       #
#  .mu             jxd matrix of mean values each component               #
#  .Sigma          dxdxK array of covariance matrices each component      #
#  .I              integral of pdf each component before applying weight  #
#  .loglik         log-likelihood of normal mixture                       #
#  .AIC            Akaike information criterion of normal mixture         #
#  .BIC            Bayesian information criterion of normal mixture       #
#                                                                         #
# REFERENCE:                                                              #
#   Volpi, E., G. Schoups, G. Firmani, and J.A. Vrugt (2017), Sworn       #
#       testimony of the model evidence: Gaussian Mixture Importance      #
#       (GAME) sampling, Water Resources Research, 53, pp. 6133-6158,     #
#       https://doi.org/10.1002/2016WR020167                              #
#                                                                         #
#  ###################################################################### #
#                                                                         #
# NOTE 1:                                                                 #
#  Geometric bridge sampling the code uses:                               #
#      t = linspace(1,0,steps);                                           #
#  Optimal bridge sampling the code uses:                                 #
#      nn1 = size(logQ_Xp(R1),1) * linspace(1,0.001,steps);               #
#  If t = 0; --> RIS whereas if t = 1 --> IS                              #
#      t in [0,1] bridge between both end members                         #
#                                                                         #
# NOTE 2:                                                                 #
#  Use a proper likelihood, that is, normalization constant included!!    #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  Step 1: Run example 2 of DREAM package                                 #
#  Step 2: Then after completion do: P = genparset(chain);                #
#  Step 3: Get posterior samples: X = P(end-25000:end,1:DREAMPar.d+2);    #
#  Step 4: Calculate marginal likelihood via GAME_sampling                #
#             [Z,logZ] = GAME_sampling(X,method,DREAMPar,Func_name, ...   #
#                 GAMEoptions,Par_info);                                  #
#             where method = any option from {'ris','is','ob','gb'}       #
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
#  Version 2.1    Feb 2025                                                #
#                                                                         #
# ####################################################################### #

# This module includes exact copies of the functions Evaluate_target, f_handle, worker, worker_task from DREAM_Suite

import os

# For parallel workers: add the directory containing 'DREAM_Suite_functions.py' to the Python path
# module_path = os.getcwd()
# if module_path not in sys.path:
#     sys.path.append(module_path)

# Get the current working directory
# parent_directory = os.path.join(module_path, 'miscellaneous')
# sys.path.append(parent_directory)

import numpy as np
import random
import importlib
from scipy.stats import multivariate_normal, t, norm, gamma, uniform, genextreme, genpareto
from scipy.linalg import cholesky
from sympy import primerange
import shutil
import importlib
import multiprocess as mp 
from DREAM_Suite_functions import DREAM_Suite_calc_setup, Evaluate_target, worker_task, Eval_prior, \
    Calc_likelihood, setup_lik, check_prior, X_unnormalize, Discrete_space, \
        distribute_tasks, convert_memoryview_to_array, X_unnormalize, \
            X_normalize, copy_model_files, cleanup_worker_directories

# Code has been very well tested and should provide the exact same results as the MATLAB code [= benchmark]
def GAME_sampling(X, method, DREAMPar, Func_name, GAMEoptions=None, Par_info=None, Meas_info=None, options=None, LV=None, MAP_info=None, plugin=None):
    # Initialize necessary structures and warnings
    if plugin is None:
        print('GAME_sampling WARNING: plugin is empty: no additional input to target function')
        plugin = None
    if MAP_info is None:
        print('GAME_sampling WARNING: MAP_info is empty: no sandwich correction')
        MAP_info = None
    if LV is None:
        print('GAME_sampling WARNING: Likelihood variables (LV) is empty: no use of distribution-adaptive likelihood functions')
        LV = None
    if options is None:
        print('GAME_sampling WARNING: Structure options is empty')
        options = {'restart': 'no', 'parallel': 'no'}
    if Meas_info is None:
        print('GAME_sampling WARNING: Structure Meas_info is empty')
        Meas_info = {}
    if Par_info is None or 'min' not in Par_info:
        print('GAME_sampling WARNING: Par_info.min not specified: Minimum values of parameters set to -inf')
        Par_info = {} 
        Par_info['min'] = -np.inf * np.ones(int(DREAMPar['d']))
    if 'max' not in Par_info:
        print('GAME_sampling WARNING: Par_info.max not specified: Maximum values of parameters set to +inf')
        Par_info['max'] = np.inf * np.ones(int(DREAMPar['d']))
    if GAMEoptions is None:
        print('GAME_sampling WARNING: Structure GAMEoptions is empty: Default settings assumed')
        GAMEoptions = {}

    # Initialize a set to track warnings for Evaluate_target
    printed_warnings = set()
    
    # Check input variables
    method, DREAMPar, Par_info, options, GAMEoptions, R, d2 = GAME_check(method, X, Func_name, DREAMPar, GAMEoptions, Par_info, Meas_info, options)

    # Initialize main variables
    DREAMPar, metric, K, N, M, steps, logZ, Par_info, Meas_info, Lik_info, options = GAME_setup(method, DREAMPar, Func_name, GAMEoptions, Par_info, Meas_info, options, LV)

    X_un = X_unnormalize(X[:, :DREAMPar['d']], Par_info)                            # Unnormalized posterior parameter values
    logP = np.sum(X[:R, :DREAMPar['d']:d2], axis = 1)                               # Logarithmic posterior density
    X = X_un[:R, :DREAMPar['d']]                                                    # Extract posterior samples
    gmix = GAME_fit_mixture(X, logP, metric, K, Par_info['min'], Par_info['max'])   # Fit Gaussian mixture
    logQ = GAME_mixture_density(gmix, X)                                            # Log density of importance distribution

    for j in range(M):              ## Dynamic part (M iterations)
        R1 = np.random.choice(np.arange(0, R//2), size = min(R//2, N), replace = False)
        if method == 'ris':         ## Reciprocal Importance Sampling
            logZ[j] = GAME_logevidence((0), [], [], logQ[R1], logP[R1])
            print('GAME_sampling WARNING: RIS samples taken from first half of x_p: make sure you use second half for mixture fitting')
        else:                       ## Importance Sampling, Geometric Bridge, Optimal Bridge Sampling
            DREAMPar['N'] = N
            DREAMPar, T, func_handle, base_dir = DREAM_Suite_calc_setup(DREAMPar, Func_name, options)
            X_gmix = mvgmmrnd(gmix, N)
            if 'steps' in Par_info:
                print('GAME_sampling WARNING: Discrete importance sampling as field steps of structure Par_info is specified')
                Par_info['step_size'] = (Par_info['max'] - Par_info['min']) / Par_info['steps']
                X_gmix = Discrete_space(X_gmix, Par_info)                           # Transform samples into discrete space

            logQ_gmix = GAME_mixture_density(gmix, X_gmix)
            FX_gmix, _ = Evaluate_target(X_gmix, DREAMPar, func_handle, Meas_info, options, base_dir, plugin, printed_warnings)
            logPR_gmix = Eval_prior(X_gmix, [], Par_info, Meas_info, options)
            logL_gmix, _, _, _, _, _ = Calc_likelihood(X_gmix, FX_gmix, DREAMPar, Par_info, Meas_info, Lik_info, options, 6, MAP_info)
            logP_X_gmix = logPR_gmix + logL_gmix                                    # Log post. densty mixture pts.

            if method == 'is':      ## Importance Sampling
                logZ[j] = GAME_logevidence((1), logQ_gmix, logP_X_gmix, logQ[R1], logP[R1])
            elif method == 'gb':    ## Geometric Bridge Sampling
                t = np.linspace(1, 0, steps)
                for k in range(steps):
                    logZ[j, k] = GAME_logevidence(t[k], logQ_gmix, logP_X_gmix, logQ[R1], logP[R1])
            elif method == 'ob':    ## Optimal Bridge Sampling
                logZ[j] = GAME_logevidence((100), logQ_gmix, logP_X_gmix, logQ[R1], logP[R1])
                nn1 = len(logQ[R1]) * np.linspace(1, 0.001, steps)
                for k in range(steps):
                    logZ[j, k] = GAME_logevidence((100), logQ_gmix, logP_X_gmix, logQ[R1], logP[R1], nn1[k])

    Z = np.exp(logZ)                                    # Transform to model evidence
    GAME_end(DREAMPar, options, base_dir)               # Terminate program    

    return Z, logZ, gmix


def GAME_check(method, Xp, Func_name, DREAMPar, GAMEoptions, Par_info, Meas_info, options):
    # ####################################################################### #
    # This function verifies the arguments defined by user for GAME sampler   #
    #                                                                         #
    # SYNOPSIS: [method,DREAMPar,Par_info,options,GAMEoptions,R,d2] = ...     #
    #               GAME_check(method,Xp,Func_name,DREAMPar,GAMEoptions, ...  #
    #               Par_info,Meas_info,options)                               #
    #                                                                         #
    # © Written by Jasper A. Vrugt, Jan 2015                                  #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #

    # Determine the number of samples of Xp
    R, d2 = Xp.shape
    # added to python code
    DREAMPar['d'] = int(DREAMPar['d'])

    if (DREAMPar['d'] + 2 != d2):
        raise ValueError(f"GAME_sampling ERROR: Matrix Xp should be of size R x DREAMPar.d+2")
    
    # <><><><><><><><><><><><><><><><> method <><><><><><><><><><><><><><><><><
    if not isinstance(method, str):
        raise TypeError('GAME_sampling ERROR: method should be a string')
    method = method.lower()

    if method not in ['ris', 'is', 'gb', 'ob']:
        print("GAME_sampling WARNING: Unknown marginal likelihood estimation method. Defaulting to 'ris'.")
        method = 'ris'

    # <><><><><><><><><><><> GAMEoptions & other inputs <><><><><><><><><><><><
    if 'metric' in GAMEoptions:
        if not isinstance(GAMEoptions['metric'], str):
            raise TypeError("GAME_sampling ERROR: Field metric of GAMEoptions should be a string with 'bic' or 'var'.")
        GAMEoptions['metric'] = GAMEoptions['metric'].lower()

        if GAMEoptions['metric'] not in ['bic', 'var']:
            print("GAME_sampling WARNING: Unknown metric in field 'metric' of GAME_options. Defaulting to 'bic'.")
            GAMEoptions['metric'] = 'bic'

    # Validation for other fields in GAMEoptions
    for field in ['J', 'N', 'M', 'steps']:
        if field in GAMEoptions and isinstance(GAMEoptions[field], str):
            raise TypeError(f"GAME_sampling ERROR: The field {field} of GAMEoptions should be an integer.")

    # Check if input data structures are dictionaries (equivalent to MATLAB structs)
    required_structs = [DREAMPar, Par_info, Meas_info, options]
    for struct in required_structs:
        if not isinstance(struct, dict):
            raise TypeError(f"GAME_sampling ERROR: input argument {struct} should be a dictionary.")

    # Normalize field names in DREAMPar, Par_info, and options to lowercase
    def normalize_fields(struct):
        return {k.lower(): v for k, v in struct.items()}
    
    DREAMPar = normalize_fields(DREAMPar)
    Par_info = normalize_fields(Par_info)
    options = normalize_fields(options)

    # <><><><><><><><><><><><><><><><> Func_name <><><><><><><><><><><><><><><><>
    if not Func_name:
        raise ValueError("GAME_sampling ERROR: The variable Func_name has to be defined as string.")
    
    if isinstance(Func_name, (int, float)):
        raise ValueError(f"GAME_sampling ERROR: The variable Func_name is defined as a numerical value ({Func_name}). It should be a string.")

    # <><><><><><><><><><><><><><><><> DREAMPar <><><><><><><><><><><><><><><><><
    if 'd' not in DREAMPar:
        raise KeyError("GAME_sampling ERROR: Field 'd' of structure DREAMPar undefined.")
    
    if DREAMPar['d'] <= 0:
        raise ValueError("GAME_sampling ERROR: Number of parameters should be integer and larger than zero -> Set DREAMPar.d >= 1")

    if DREAMPar['lik'] not in [1, 2, *range(11, 19), 21]:
        raise ValueError("GAME_sampling ERROR: Cannot compute marginal likelihood for this type of (informal) likelihood/posterior density function.")

    if DREAMPar['lik'] in [12, 13, 16] and 'Sigma' not in Meas_info:
        raise KeyError("GAME_sampling ERROR: Meas_info.Sigma needs to be defined either an inline function or a numerical value.")

    # Check if 'Y' field in Meas_info is defined and valid
    if DREAMPar['lik'] in range(11, 19):
        if 'Y' not in Meas_info or len(Meas_info['Y']) == 0:
            raise KeyError("GAME_sampling ERROR: Field 'Y' of structure Meas_info has to be defined (stores training observations).")
        if len(Meas_info['Y'].shape) > 1:
            raise ValueError("GAME_sampling ERROR: Field 'Y' of structure Meas_info has to be a column vector.")

    # More checks for other DREAMPar conditions
    if DREAMPar['lik'] == 18 and 'C' not in Meas_info:
        raise KeyError("GAME_sampling ERROR: Field 'C' of structure Meas_info has to be defined for likelihood 18.")

    # Additional checks for 'S' field in Meas_info
    if DREAMPar['lik'] == 21 and 'S' not in Meas_info:
        raise KeyError("GAME_sampling ERROR: Field 'S' of structure Meas_info has to be defined (stores simulated summary metrics).")

    # <><><><><><><><><><><><><><><><> Par_info <><><><><><><><><><><><><><><><>
    # Check and normalize Par_info
    if 'boundhandling' not in Par_info:
        Par_info['boundhandling'] = 'none'

    # Check: do we sample in normalized [0-1] space, or not?
    if 'norm' in Par_info:
        Par_info['norm'] = int(Par_info['norm'])
        if not isinstance(Par_info['norm'], (int, float)):
            raise ValueError("GAME_sampling ERROR: Par_info.norm should be a scalar -> Define Par_info.norm = 0 or 1")
        if Par_info['norm'] not in [0, 1]:
            raise ValueError("GAME_sampling ERROR: Par_info.norm should be zero or one -> Define Par_info.norm = 0 or 1")
        if Par_info['norm'] == 1:
            if 'min' not in Par_info:
                raise ValueError("GAME_sampling ERROR: Parameter normalization is used but minimum parameter values not defined -> Set Par_info.min!!")
            if 'max' not in Par_info:
                raise ValueError("GAME_sampling ERROR: Parameter normalization is used but maximum parameter values not defined -> Set Par_info.max!!")
    else:
        Par_info['norm'] = (0)

    # Validate Par_info min/max lengths match DREAMPar.d
    if len(Par_info.get('min', [])) != DREAMPar['d'] or len(Par_info.get('max', [])) != DREAMPar['d']:
        raise ValueError(f"GAME_sampling ERROR: Par_info.min and Par_info.max should be of length {DREAMPar['d']}.")

    # <><><><><><><><><><><><><><><><> Meas_info <><><><><><><><><><><><><><><><>
    if 'Sigma' in Meas_info:
        if np.any(Meas_info['Sigma'] < 0):
            print("GAME_sampling WARNING: One or more entries of Meas_info.Sigma is negative - we use absolute values.")
            Meas_info['Sigma'] = np.abs(Meas_info['Sigma'])
        if np.any(Meas_info['Sigma'] == 0):
            print("GAME_sampling WARNING: One or more entries of Meas_info.Sigma are zero - we replace them with the smallest non-zero value.")
            Meas_info['Sigma'][Meas_info['Sigma'] == 0] = 1e-3

    # <><><><><><><><><><><><><><><>< options ><><><><><><><><><><><><><><><><>
    if 'epsilon' in options and options['epsilon'] <= 0:
        raise ValueError("GAME_sampling ERROR: Value of 'epsilon' of structure options should be larger than zero.")

    if 'rho' in options and not callable(options['rho']):
        raise TypeError("GAME_sampling ERROR: Field 'rho' of structure options should be defined as an inline function.")

    # Check content of each field of structure options
    for key, value in options.items():
        if key not in ['epsilon', 'rho', 'burnin'] and value not in ['yes', 'no']:
            raise ValueError(f"GAME_sampling ERROR: Field {key} of structure options should be 'yes' or 'no'.")

    # Parallel computation check
    if 'parallel' in options and options['parallel'] == 'yes':
        if not hasattr(os, 'distributed'):
            print("GAME_sampling WARNING: Distributed Computing toolbox not available. Defaulting to 'no'.")
            options['parallel'] = 'no'

    return method, DREAMPar, Par_info, options, GAMEoptions, R, d2


def GAME_setup(method, DREAMPar, Func_name, user_options, Par_info, Meas_info, options, LV):
    # ####################################################################### #
    # Initializes the main variables used in GAME sampling                    #
    #                                                                         #
    # SYNOPSIS: [DREAMPar,metric,K,N,M,steps,Par_info,Meas_info,Lik_info, ... #
    #               GAMEoptions,slash_dir] = GAME_setup(DREAMPar, ...         #
    #               Func_name,user_options,Par_info,Meas_info,options)        #
    #                                                                         #
    # © Written by Jasper A. Vrugt, Aug. 2015                                 #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #

    random.seed(1 + round(100 * random.random()))

    # Default values for missing options
    name = ['parallel', 'IO', 'modout', 'save', 'restart', 'DB', 'epsilon']
    value = ['no', 'no', 'no', 'no', 'no', 'no', '0.025']
    
    # Set default values
    for i, var_name in enumerate(name):
        if var_name not in options:
            options[var_name] = value[i]

    # Default settings for GAME options
    def_GAMEoptions = {'metric': 'bic', 'K': 5, 'M': 1, 'N': 1e4, 'steps': 10}
    
    # Update GAME options with user options
    for key, value in user_options.items():
        if key in def_GAMEoptions:
            if isinstance(value, str):  # Only apply 'lower()' to strings
                def_GAMEoptions[key] = value.lower()
            else:
                def_GAMEoptions[key] = value

    # Unpack GAME options into individual variables
    metric = def_GAMEoptions['metric']
    K = int(def_GAMEoptions['K'])
    N = int(def_GAMEoptions['N'])
    M = int(def_GAMEoptions['M'])
    steps = int(def_GAMEoptions['steps'])
    if method in ['ris','is']:
        steps = 1

    # Initialize DREAMPar.CPU to be 1
    DREAMPar['CPU'] = 1

    # Approximate Bayesian computation
    if 20 < DREAMPar['lik'] < 24:
        options['ABC'] = 'yes'
        if 'rho' not in options:
            options['rho'] = lambda X, Y: abs(X - Y)
    else:
        options['ABC'] = 'no'

    # Check parameter normalization
    if Par_info['norm'] == 1:
        Par_info['minun'] = Par_info['min']
        Par_info['maxun'] = Par_info['max']
        Par_info['min'] = np.zeros(1, DREAMPar['d'])
        Par_info['max'] = np.ones(1, DREAMPar['d'])

    # Set number of observations
    Meas_info['n'] = len(Meas_info.get('Y', []))

    # Check summary metrics
    if 'S' in Meas_info:
        Meas_info['n_S'] = len(Meas_info['S'])
        options['epsilon'] = np.ones(Meas_info['n_S']) * options['epsilon']
    else:
        Meas_info['n_S'] = int(0)

    # Handling Sigma
    if 'Sigma' in Meas_info:
        if np.isscalar(Meas_info['Sigma']):
            Meas_info['Sigma'] = np.ones(Meas_info['n']) * Meas_info['Sigma']
    else:
        Meas_info['Sigma'] = []

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
        print(f"GAME_sampling: Code analyzed that prior handle returns a {Par_info['pr']}")

    # Setup likelihood function
    Lik_info, DREAMPar, Par_info = setup_lik(Func_name, DREAMPar, Par_info, Meas_info, LV)

    # Order fields
    DREAMPar = dict(sorted(DREAMPar.items()))
    Par_info = dict(sorted(Par_info.items()))
    Meas_info = dict(sorted(Meas_info.items()))
    options = dict(sorted(options.items()))

    # Initialize logZ
    if steps == 1:
        logZ = np.nan * np.ones(M)
    else:
        logZ = np.nan * np.ones((M,steps))
        
    # Print welcome message
    print('  ----------------------------------------------------------------------------------------------------------------------------------------             ')
    print("  GGGGGGGGG     AAAA     MMM      MMM  EEEEEEEE     SSSSSSSSS     AAAA     MMM      MMM  PPPPPPPP   LLL        III  NNN     NNN  GGGGGGGGG             ")
    print("  GGGGGGGGG    AAAAAA    MMM      MMM  EEEEEEEE     SSSSSSSS     AAAAAA    MMM      MMM  PPPPPPPPP  LLL        III  NNN     NNN  GGGGGGGGG             ")
    print("  GGG   GGG   AAA  AAA   MMM      MMM  EEE          SSS         AAA  AAA   MMM      MMM  PPP    PPP LLL        III  NNNN    NNN  GGG   GGG             ")
    print("  GGG   GGG  AAA    AAA  MMMM    MMMM  EEE          SSS        AAA    AAA  MMMM    MMMM  PPP    PPP LLL        III  NNNNN   NNN  GGG   GGG             ")
    print("  GGGGGGGGG  AAA    AAA  MMMMM  MMMMM  EEEEE   ---- SSSSSSSSS  AAA    AAA  MMMMM  MMMMM  PPPPPPPPP  LLL        III  NNN NN  NNN  GGGGGGGGG     /^ ^\   ")
    print("   GGGGGGGG  AAAAAAAAAA  MMM MMMM MMM  EEEEE   ---- SSSSSSSSS  AAAAAAAAAA  MMM MMMM MMM  PPPPPPPP   LLL        III  NNN  NN NNN   GGGGGGGG    / 0 0 \  ")
    print("        GGG  AAAAAAAAAA  MMM  MM  MMM  EEE                SSS  AAAAAAAAAA  MMM  MM  MMM  PPP        LLL        III  NNN   NNNNN        GGG    V\ Y /V  ")
    print("        GGG  AAA    AAA  MMM      MMM  EEE                SSS  AAA    AAA  MMM      MMM  PPP        LLL        III  NNN    NNNN        GGG     / - \   ")
    print("   GGGGGGGG  AAA    AAA  MMM      MMM  EEEEEEEE      SSSSSSSS  AAA    AAA  MMM      MMM  PPP        LLLLLLLLL  III  NNN     NNN   GGGGGGGG    /     |  ")
    print("  GGGGGGGGG  AAA    AAA  MMM      MMM  EEEEEEEE     SSSSSSSSS  AAA    AAA  MMM      MMM  PPP        LLLLLLLLL  III  NNN     NNN  GGGGGGGGG    V__) ||  ")
    print('  ----------------------------------------------------------------------------------------------------------------------------------------             ')
    print('  © Jasper A. Vrugt, University of California Irvine & GPT-4 OpenAI''s language model')
    print('    ________________________________________________________________________')
    print('    Version 2.1, Feb. 2025, Beta-release: MATLAB implementation is benchmark')
    print('\n')

    return DREAMPar, metric, K, N, M, steps, logZ, Par_info, Meas_info, Lik_info, options


def GAME_mixture_density(gmix, X, cov_type = 2, shared = 0):
    # ####################################################################### #
    # Computes log density of normal mixture importance distribtion           #
    #                                                                         #
    # SYNOPSIS:                                                               #
    #  logpdf = GAME_mixture_density(gmix,X,cov_type)                         #
    # WHERE                                                                   #
    #  gmix        [input] Structure trained normal mixture importance dist.  #
    #  .n              # samples (= R) of X                                   #
    #  .d              # dimensions (= parameters) of X                       #
    #  .k              # components of mixture distribution                   #
    #  .p              # parameters of mixture distribution                   #
    #  .w              maximum likelihood weights of mixture components       #
    #  .mu             jxd matrix of mean values each component               #
    #  .Sigma          dxdxj array of covariance matrices each component      #
    #  .I              integral of pdf each component before applying weight  #
    #  .loglik         log-likelihood of normal mixture                       #
    #  .AIC            Akaike information criterion of normal mixture         #
    #  .BIC            Bayesian information criterion of normal mixture       #
    #  X          [input] Nxd matrix of samples                               #
    #  cov_type   [input] OPT: Covariance type                                #
    #   = 1            Diagonal covariance matrix [diagonal entries only]     #
    #   = 2            Full covariance matrix [= all entries are estimated]   #
    #  shared     [input] OPT: shared covariance matrix Σ mixture components  #
    #   = 0            Mixture components share same covariance matrix Σ      #
    #   = 1            Mixture components have different covariance matrices  #
    #  logpdf     [outpt] Nx1 vector logarithmic density mixture at X samples #
    #                                                                         #
    # © Written by Jasper A. Vrugt, Jan 2015                                  #
    #   Adapted from wdensity                                                 #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #

    # mu is d x k matrix with k components
    w = gmix['w']                                                                       # Component weights
    I = gmix['I']                                                                       # Integral of componet PDFs
    mu = gmix['mu'].T                                                                   # Copy of mu
#    Sigma = np.reshape(gmix['Sigma'],(2,2,int(gmix['k'])))                              # Copy of Sigma
    Sigma = gmix['Sigma']                                                               # Copy of Sigma
    known = 0                                                                           # Initialized, Python
    log_pr = np.log(w)                                                                  # Log value component weights
    n, d = X.shape                                                                      # Number of samples n, number of dimensions d
    ell = np.full((n, gmix['k']), np.nan)                                               # Log-likelihood for mixture at X
    for k in range(0, gmix['k']):
        if shared == 1:                                                                 # Components share the same covariance matrix
            z = 0                                                                       # Same covariance matrix for all components
            if k > 0:
                known = 1
        else:
            z = k                                                                       # Different covariance matrix for each component
        if cov_type == 1 and known == 0:                                                # Diagonal Σ matrix 
            L = np.sqrt(Sigma[:, :, z])                                                 # Can take sqrt, Sigma = L'*L    
            if np.any(L < np.finfo(float).eps * d):                                     # Check L entries
                raise ValueError("GAME_mixture_density: Ill-conditioned covariance")
            logdetSigma = np.sum(np.log(Sigma[:, :, z]))                                # Compute log(|Sigma(:,:,z)|)

        elif cov_type == 2 and known == 0:                                              # Full Σ matrix
            L, f = np.linalg.cholesky(Sigma[:, :, z]), False                            # Cholesky fact., Sigma = L'*L
            diagL = np.diagonal(L)                                                      # Diagonal L entries
            if f != 0 or np.any(np.abs(diagL) < np.finfo(float).eps * L.shape[1]):      # Check diagL entries
                raise ValueError("GAME_mixture_density: Ill-conditioned covariance")
            logdetSigma = 2 * np.sum(np.log(diagL))                                     # Compute log(|Sigma(:,:,z - 1)|)

        else:
            raise ValueError("GAME_mixture_density: unknown covariance type")

        Xc = X.copy() - mu[k, :]                                                        # Center samples by subtracting the mean
        if cov_type == 2:
            XcRinv = np.linalg.solve(L, Xc.T).T                                         # Apply the inverse of L (Cholesky factor)
        else:
            XcRinv = Xc * (1. / L )                                                     # Element-wise division for diagonal covariance

        d_M = np.sum(XcRinv ** 2, axis = 1)                                             # nx1 vector of Mahalanobis distances
        ell[:, k] = -0.5 * d_M - 0.5 * logdetSigma + log_pr[k] \
            - 0.5 * d * np.log(2 * np.pi) - np.log(I[k])                                # Loglik each sample for each component

    ell_max = np.max(ell, axis = 1)                                                     # Max log-likelihood for each sample
    p = np.exp(ell - ell_max[:, None])                                                  # Subtract ell_max to avoid underflow
    pdf = np.sum(p, axis = 1)                                                           # Sum the weighted probabilities for each sample
    logpdf = np.log(pdf) + ell_max                                                      # Final log-densities

    return logpdf


def GAME_logevidence(t, logq00, logq10, logq01, logq11, n1 = None):
    # ####################################################################### #
    # Bridge sampling estimation of Bayesian model evidence. This function    #
    # will compute the log normalizing constant of q1, which is the           #
    # logarithmic value of the Bayesian model evidence                        #
    #                                                                         #
    # MAIN IDEA: Take q1 as unnormalized posterior, and q0 as "importance"    #
    # distribution, i.e. a normalized density that ideally approximates the   #
    # shape of q1                                                             #
    #                                                                         #
    # SYNOPSIS: logZ = GAME_logevidence(t,logq00,logq10,logq01,logq11,n1)     #
    #                                                                         #
    # WHERE                                                                   #
    #  t           [input] specifies the type of estimator                    #
    #   = 0            reciprocal imprtnce smplng [= harm est. if q0 = prior] #
    #   = 1            importance sampling                                    #
    #   = [0,1]        geometric bridge sampling                              #
    #   = otherwise    optimal bridge sampling                                #
    #  logq00      [input] n0x1 vector log(q0) n0 samples q0 [not used t=0]   #
    #  logq10      [input] n0x1 vector log(q1) n0 samples q0 [not used t=0]   #
    #  logq01      [input] n1x1 vector log(q0) n1 samples q1 [not used t=1]   #
    #  logq11      [input] n1x1 vector log(q1) n1 samples q1 [not used t=1]   #
    #  n1          [input] OPT: # samples of q1                               #
    #  logZ        [outpt] logarithmic value of Bayesian model evidence       #
    #                                                                         #
    # © Written by Gerrit Schoups, July 2016                                  #
    # Delft Technical University, The Netherlands                             #
    # Modified by Jasper A. Vrugt, Aug. 2016                                  #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #

    if n1 is None:
        n1 = len(logq01)
    
    if t == 0:              ## qhalf = q0 (RIS)
        logZ = -logmeanexp(logq01 - logq11)
    elif t == 1:            ## qhalf = q1 (IS)
        logZ = logmeanexp(logq10 - logq00)
    elif 0 < t < 1:         ## qhalf = geometric bridge = q0^(1-t) * q1^t
        logZ = logmeanexp(t * (logq10 - logq00)) - logmeanexp((1 - t) * (logq01 - logq11))
    else:                   ## qhalf = optimal bridge = 1 / (Z*s0/q1 + s1/q0)
        n0 = len(logq00)
        logn = np.log(n0 + n1)
        logs0 = np.log(n0) - logn
        logs1 = np.log(n1) - logn
        logl0 = logq10 - logq00
        logl1 = logq11 - logq01
        # initial guess
        if t == 10:         ## use RIS as initial guess
            logZ = -logmeanexp(-logl1)
        elif 10 < t < 11:   ## use geometric bridge as initial guess
            t = t - 10
            logZ = logmeanexp(t * logl0) - logmeanexp(- (1 - t) * logl1)
        else:               ## use IS as initial guess
            logZ = logmeanexp(logl0)
        
        # iterate
        for _ in range(10):
            logs0Z0 = (logs0 + logZ) * np.ones((n0))
            logs0Z1 = (logs0 + logZ) * np.ones((len(logq01)))                       # ones(n1)
            logdenom0 = logsumexp(np.column_stack([logs0Z0, logs1 + logl0]), 1)     # log(s0*Z + s1*l0)
            logdenom1 = logsumexp(np.column_stack([logs0Z1, logs1 + logl1]), 1)     # log(s0*Z + s1*l1)
            logZ = logmeanexp(logl0 - logdenom0) - logmeanexp(-logdenom1)
    
    return logZ


def GAME_fit_mixture(X, logP, metric, K, lower, upper):
    # ####################################################################### #
    # This function determines the optimal Gaussian mixture distribution      #
    # for Nxd matrix of posterior samples, X, derived from (e)DREAM Package.  #
    # We try mixture distributions consisting of k = 1,...,K components. Each #
    # of these mixtures is trained using the expectation maximization method. #
    # Then, the optimal mixture distribution is determined using model        #
    # selection criteria such as the AIC, BIC and variance ratio, 'var', of   #
    # the target (= posterior) and Gaussian mixture PDF.                      #
    #                                                                         #
    # The normal mixture PDF of x \in X [= 1xd vector] is written as          #
    #          PDF(x) = Σ_{j=1}^{k} w_{j} f_{N_{d}}(x,µ_{j},Σ_{j})            #
    #                                                                         #
    # where f_{N_{d}}(x,µ_{j},Σ_{j}) is the d-variate normal pdf with mean µ  #
    # [= 1xd vector] and dxd covariance matrix Σ and the weights w_{1},...,   #
    # w_{k} must sum to 1. The weights, w_{j}, mean vector, µ_{j}, and dxd    #
    # covariance matrix, Σ_{j}, of each mixture component, j = 1,...,k, are   #
    # determined using Expectation-Maximization. The default option is that   #
    # each mixture component has its own and full covariance matrix. Then,    #
    # the EM algorithm has to train p1 = d - 1 weights, namely, w_{1},...,    #
    # w_{k-1}, p2 = k*d mean entries, µ_{1},...,µ_{k}, and p3 = k*d*(d+1)/2   #
    # elements of the covariance matrices, Σ_{1},...,Σ_{k}. The total # of    #
    # parameters of the normal mixture equals p = p1+p2+p3. In MATLAB, we     #
    # yield that f_{N_{d}}(x,µ_{j},Σ_{j}) = mvnpdf(x,µ_{j},Σ_{j})             #
    #                                                                         #
    # SYNOPSIS:                                                               #
    #  gmix = GAME_fit_mixture(X,logP,metric,K,lower,upper)                   #
    # WHERE                                                                   #
    #  X           [input] Rxd matrix posterior samples DREAM                 #
    #  logP        [input] Rx1 vector of logarithmic values posterior density #
    #  metric      [input] string (name) metric used for mixture selection    #
    #   = 'bic'        Bayesian information criterion        DEFault          #
    #   = 'var'        Variance reduction                                     #
    #  K           [input] Maximum # components mixture      DEF: 5           #
    #                      [= importance] distribution                        #
    #  lower       [input] 1xd vector lower bounds parameters                 #
    #  upper       [input] 1xd vector upper bounds parameters                 #
    #  gmix        [outpt] Structure trained normal mixture importance dist.  #
    #  .n              # samples (= R) of X                                   #
    #  .d              # dimensions (= parameters) of X                       #
    #  .K              # components of mixture distribution                   #
    #  .p              # parameters of mixture distribution                   #
    #  .w              maximum likelihood weights of mixture components       #
    #  .mu             jxd matrix of mean values each component               #
    #  .Sigma          dxdxj array of covariance matrices each component      #
    #  .I              integral of pdf each component before applying weight  #
    #  .loglik         log-likelihood of normal mixture                       #
    #  .AIC            Akaike information criterion of normal mixture         #
    #  .BIC            Bayesian information criterion of normal mixture       #
    #                                                                         #
    # © Written by Jasper A. Vrugt, Jan 2015                                  #
    # University of California Irvine                                         #
    #                                                                         #
    # Version 1:    June 2012       Initial setup and definition              #
    # Version 1.1:  Jan. 2015       Major overhaul of code                    #
    # Version 1.2:  Feb. 2015       VAR method: min. variance ratio q and p   #
    #				Oct. 2015       Add bounded Gaussian mixture distribution #
    #                                                                         #
    # ####################################################################### #

    var_ratios = np.inf * np.ones(K)  # Variance of the ratio between q and p
    bic = np.inf * np.ones(K)         # Initialize mixture BIC
    gmix = []
    
    for k in range(1, K+1):  # Loop over the number of components in the mixture
        #try:
            # Fit a Gaussian mixture using custom EM algorithm
        g = emgm(X.T.copy(), k, lower, upper)[0]
        logQ = GAME_mixture_density(g, X)
        var_ratios[k-1] = np.var(logP - logQ)
        bic[k-1] = g['BIC']
        gmix.append(g)
        #except Exception as e:
            # If fitting fails, print the error and continue with the next iteration
        #    print(f"GAME_fit_mixture: Error fitting Gaussian mixture for k={k}: {e}")
        #    continue

    # Sort BIC in ascending order to find the best mixture [= lowest BIC]
    j1 = np.argmin(bic[:K])         # Index of best mixture by BIC
    j2 = np.argmin(var_ratios[:K])  # Index of best mixture by variance ratio

    # Return the mixture based on the selected metric
    if metric == 'bic':
        gmix = gmix[j1]
    elif metric == 'var':
        gmix = gmix[j2]

    # Check if BIC and VAR metrics yield the same result
    if j1 == j2:
        print("GAME_fit_mixture: BIC and VAR metric lead to the same optimal mixture distribution.")
    else:
        print("GAME_fit_mixture: BIC and VAR metric do not yield the same optimal mixture distribution. Check results with GAMEoptions.metric = 'var'.")

    return gmix


def GAME_end(DREAMPar, options, base_dir):
    # ####################################################################### #
    # This function close the MATLAB pool (if CPU > 1) and/or removes files   #
    #                                                                         #
    # SYNOPSIS: GAME_end(DREAMPar,options)                                    #
    #                                                                         #
    # © Written by Jasper A. Vrugt, Jan 2015                                  #
    # University of California Irvine                                         #
    #                                                                         #
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
        fid.write('----------- End of GAME_sampling warning file ----------\n')

    # Optionally open the warning file in the default editor for the platform
    # if platform.system() in ['Windows', 'Darwin']:  # Windows or macOS
    #     os.system('notepad warning_file.txt' if platform.system() == 'Windows' else 'open warning_file.txt')

    return



def mvgmmrnd(gmix, N):
    # ####################################################################### #
    #MVGMMRND Multivariate Gaussian mixture model random. This function draws #
    # N samples from the normal mixture distribution defined in the structure #
    # gmix. The normal mixture PDF at x [= 1xd vector] is written as          #
    #          PDF(x) = Σ_{j=1}^{k} w_{j} f_{N_{d}}(x,µ_{j},Σ_{j})            #
    # where f_{N_{d}}(x,µ_{j},Σ_{j}) is the d-variate normal pdf with mean µ  #
    # [= 1xd vector] and dxd covariance matrix Σ. In MATLAB, we yield that    #
    # f_{N_{d}}(x,µ_{j},Σ_{j}) = mvnpdf(x,µ_{j},Σ_{j})                        #
    #                                                                         #
    # SYNOPSIS:                                                               #
    #  X = mvgmmrnd(gmix,N);                                                  #
    # WHERE                                                                   #
    #  gmix        [outpt] Structure trained normal mixture importance dist.  #
    #  .d              # dimensions (= parameters) of X                       #
    #  .k              # components of mixture distribution                   #
    #  .w              maximum likelihood weights of mixture components       #
    #  .mu             kxd matrix of mean values each component               #
    #  .Sigma          dxdxk array of covariance matrices each component      #
    #  N           [input] # samples drawn from normal mixture ditribution    #
    #  X           [outpt] Nxd matrix of samples drawn from normal mixture    #
    #                                                                         #
    # © Written by Jasper A. Vrugt, Jan. 2015                                 #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #

    nw = np.histogram(np.random.rand(N), bins = np.append([0], np.cumsum(gmix['w']) / np.sum(gmix['w'])))[0]
    cum_nw = np.append([0], np.cumsum(nw))                                              # Cumulative count
    X = np.nan * np.ones((N, gmix['d']))                                                # Initialize Nxd matrix X
    for j in range(0, gmix['k']):                                                       # Draw samples for each component
        id = np.arange(cum_nw[j], cum_nw[j + 1])                                        # Start at zeroth sample
        # X[id, :] = np.random.normal(0, 1, (nw[j], gmix['d'])) @ \
        #     np.linalg.cholesky(gmix['Sigma'][:, :, j]) + gmix['mu'][:, j].T           # Store the samples
        # Reminder: Next sentence is correct - took a long time to find this bug!!
        X[id, :] = np.random.randn(nw[j], gmix['d']) @ \
            cholesky(gmix['Sigma'][:, :, j]) + np.tile(gmix['mu'][:, j].T, (nw[j], 1))  # Store the samples
        
    return X


def emgm(X, k, lower, upper):
    # ####################################################################### #
    #EMGM Training of normal mixture model using the Expectation-Maximization #
    # algorithm                                                               #
    #                                                                         #
    # SYNOPSIS:                                                               #
    #  [gmix,loglik] = emgm(X,k,lower,upper)                                  #
    # WHERE                                                                   #
    #  X           [input] Rxd matrix posterior samples DREAM                 #
    #  k           [input] scalar with # components normal mixture dist.      #
    #  lower       [input] 1xd vector lower bounds parameters                 #
    #  upper       [input] 1xd vector upper bounds parameters                 #
    #  gmix        [outpt] Structure trained normal mixture importance dist.  #
    #  .d              # dimensions (= parameters) of X                       #
    #  .n              # samples (= R) of X                                   #
    #  .k              # components of mixture distribution                   #
    #  .p              # parameters of mixture distribution                   #
    #  .w              maximum likelihood weights of mixture components       #
    #  .mu             jxd matrix of mean values each component               #
    #  .Sigma          dxdxj array of covariance matrices each component      #
    #  .I              integral of each component pdf before applying weight  #
    #  .loglik         log-likelihood of normal mixture                       #
    #  .AIC            Akaike information criterion of normal mixture         #
    #  .BIC            Bayesian information criterion of normal mixture       #
    #  m_loglik    [outpt] 1xit vector of mean log-likelihoods mixture dist.  #
    #                                                                         #
    # © Modified by Jasper A. Vrugt, Jan 2015                                 #
    #   based on code of Michael Chen (sth4nth@gmail.com)                     #
    #                                                                         #
    # ####################################################################### #

    print("GAME_fit_mixture:EMGM: Fitting Gaussian mixture distribution...")

    # Initialize with random values for means and covariances
    R = initialization(X, k)        # X is the data matrix, k is the number of components
    label = np.argmax(R, axis=1)    # Find the labels by selecting the maximum responsibility
    R = R[:, np.unique(label)]      # Remove empty components (if any)

    tol = 1e-10
    maxiter = 2500
    m_loglik = -np.inf * np.ones(maxiter)
    converged = False
    it = 1

    while not converged and it < maxiter:
        
        gmix = maximization(X, R)
        R, m_loglik[it] = expectation(X, gmix)
        label = np.argmax(R, axis = 1)
        unique_labels = np.unique(label)
        if R.shape[1] != len(unique_labels):
            R = R[:, unique_labels]
        else:
            converged = np.abs(m_loglik[it] - m_loglik[it - 1]) < tol * np.abs(m_loglik[it])

        it += 1
        
    m_loglik = m_loglik[0:it]
    
    if converged:
        print(f"GAME_fit_mixture:EMGM: Parameters of k = {k}-component mixture distribution converged in {it-1} steps.")
    else:
        print(f"GAME_fit_mixture:EMGM: Parameters of k = {k}-component mixture distribution did not converge in {maxiter} steps.")

    # Collect final mixture parameters
    gmix['k'] = int(k)
    gmix['d'], gmix['n'] = X.shape
    gmix['loglik'] = gmix['n'] * m_loglik[it - 1]
    gmix['p'] = gmix['k'] - (1) + gmix['k'] * gmix['d'] + gmix['k'] * gmix['d'] * (gmix['d'] + (1)) // 2
    gmix['AIC'] = -2 * gmix['loglik'] + 2 * gmix['p']
    gmix['BIC'] = -2 * gmix['loglik'] + gmix['p'] * np.log(gmix['n'])

    # Compute the integral of components' pdfs for the bounded case (lower/upper bounds)
    gmix['I'] = np.nan * np.ones(k)
    precision = int(1e4)
    if gmix['d'] > 1:  # If more than 1 dimension
        low = up = np.zeros((gmix['d']))
        for k in range(gmix['k']):
            low = lower - gmix['mu'][:, k]
            up = upper - gmix['mu'][:, k]
            gmix['I'][k] = qsimvn(precision, gmix['Sigma'][:, :, k], low, up)[0]
    else:  # For 1D case
        for k in range(gmix['k']):
            gmix['I'][k] = multivariate_normal.cdf(upper, mean = gmix['mu'][:, k], cov = gmix['Sigma'][:, :, k]) - \
                           multivariate_normal.cdf(lower, mean = gmix['mu'][:, k], cov = gmix['Sigma'][:, :, k])

    # gmix['mu'] is vertical -> d x k matrix, where k is the number of components

    return gmix, m_loglik


def initialization(X, init):
    # ####################################################################### #
    #INITIALIZATION Initialization of the mixture                             #
    #                                                                         #
    # SYNOPSIS:                                                               #
    #  R = initialization(X,init)                                             #
    #                                                                         #
    # © Written by Michael Chen (sth4nth@gmail.com)                           #
    #   Modified by Jasper A. Vrugt                                           #
    #                                                                         #
    # ####################################################################### #

    d, n = X.shape                                      # d = number of features, n = number of samples

    if isinstance(init, dict):                          # Initialize with a model (assuming `gmix` structure)
        R = expectation(X, init)[0]
    
    elif isinstance(init, int):                         # Random initialization
        k = init
        idx = np.random.choice(n, k, replace=False)     # Randomly pick k indices from n
        m = X[:, idx]                                   # Choose initial means from random samples
        
        # Compute the labels based on the distance to the means
        # label = np.argmax(np.dot(m.T, X) - np.sum(m**2, axis=0).T / 2, axis=0)
        A = np.dot(m.T, X)              # k x n matrix
        B = np.sum(m**2, axis=0) / 2    # 1 x k matrix
        label = np.argmax(A - B[:,None], axis=0)
        # Ensure labels are in range [0, k-1]
        label = np.clip(label, 0, k-1)

        # Compute the labels based on the distance to the means 
        # Compute the squared Euclidean distances from each sample to each mean
        # diff = X[:, :, np.newaxis] - m  # Shape (d, n, k)
        # dist = np.sum(diff**2, axis=0)  # Shape (n, k), squared distances
        # label = np.argmin(dist, axis=1)  # Find the index of the closest mean for each sample

        # Ensure labels are in range [0, k-1]    
        unique_labels = np.unique(label)

        while k != len(unique_labels):
            idx = np.random.choice(n, k, replace=False)
            m = X[:, idx]
            # Separated the following line in two lines and added [':, None]
            # label = np.argmax(np.dot(m.T, X) - np.sum(m**2, axis=0) / 2, axis=0)
            z = np.sum(m**2, axis = 0) / 2
            label = np.argmax(np.dot(m.T, X) - z[:, None], axis = 0)
            # Clip between 0 and k - 1 (= largest component of k-components with 0 index in Python)
            label = np.clip(label, 0, k-1)
            unique_labels = np.unique(label)

        # Create responsibility matrix (n x k), with 1 for the corresponding label (= component)
        R = np.zeros((n, k), dtype=int)
        R[np.arange(n), label] = 1

    elif len(init.shape) == 1 and init.shape[0] == n:  # Initialize with labels
        label = init
        k = np.max(label)
        R = np.zeros((n, k), dtype=int)
        R[np.arange(n), label] = 1

    elif len(init.shape) == 2 and init.shape[0] == d:  # Initialize with centers
        k = init.shape[1]
        m = init
        # Separated the following line in two lines and added [':, None]          
        # label = np.argmax(np.dot(m.T, X) - np.sum(m**2, axis=0) / 2, axis=0)
        z = np.sum(m**2, axis = 0) / 2
        label = np.argmax(np.dot(m.T, X) - z[:, None], axis = 0)
        # Clip between 0 and k - 1 (= largest component of k-components with 0 index in Python)
        label = np.clip(label, 0, k-1)
        R = np.zeros((n, k), dtype=int)
        R[np.arange(n), label] = 1

    else:
        raise ValueError("GAME_fit_mixture:EMGM:Initialization: Invalid initialization.")

    return R


def expectation(X, gmix):
    # ####################################################################### #
    #EXPECTATION Expectation step of EM algorithm                             #
    #                                                                         #
    # SYNOPSIS:                                                               #
    #  [R,loglik] = expectation(X,gmix)                                       #
    #                                                                         #
    # © Written by Michael Chen (sth4nth@gmail.com)                           #
    #   Modified by Jasper A. Vrugt                                           #
    #                                                                         #
    # ####################################################################### #

    n = X.shape[1]                                  # Number of samples
    k = gmix['mu'].shape[1]                         # Number of components
    logrho = np.zeros((n, k))                       # Initialize log responsibilities

    for j in range(k):                              # Compute pdf of each mixture component
        mu = gmix['mu'][:, j]
        Sigma = gmix['Sigma'][:, :, j]
        logrho[:, j] = loggausspdf(X, mu, Sigma)    # Function to compute the log PDF

    logrho += np.log(gmix['w'])                     # Add log(weights)
    T = logsumexp(logrho, 1)                        # Log-sum-exp for normalization [ sum according to rows, dim = 2 in MATLAB]
    loglik = np.sum(T) / n                          # Log-likelihood of the mixture model

    # Normalize to get responsibilities (scaled densities)
    R = np.exp(logrho - T[:, None])                 # Broadcasting to normalize

    return R, loglik


def maximization(X, R):
    # ####################################################################### #
    #MAXIMIZATION Maximization step of EM algorithm                           #
    #                                                                         #
    # SYNOPSIS:                                                               #
    #  gmix = maximization(X,R)                                               #
    #                                                                         #
    # © Written by Michael Chen (sth4nth@gmail.com)                           #
    #   Modified by Jasper A. Vrugt                                           #
    #                                                                         #
    # ####################################################################### #
    
    d, n = X.shape
    k = R.shape[1]
    nk = np.sum(R, axis=0)

    w = nk / n
    # mu = np.dot(X, R) / nk[:, None]
    mu = np.dot(X, R) / nk
    Sigma = np.zeros((d, d, k))
    sqrtR = np.sqrt(R)

    for j in range(k):
        Xo = X - mu[:, j][:, None]
        Xo = Xo * sqrtR[:, j]
        Sigma[:, :, j] = np.dot(Xo, Xo.T) / nk[j]
        Sigma[:, :, j] += 1e-6 * np.eye(d)          # Regularization for numerical stability

    return {'mu': mu, 'Sigma': Sigma, 'w': w}


def loggausspdf(X, mu, Sigma):
    # ####################################################################### #
    #LOGGAUSSPDF Logarithmic density of Nd(µ,Σ) at N samples of X             #
    # = log(mvnpdf(X,mu,Sigma))'                                              #
    # but with protection against numerical underflow                         #
    #                                                                         #
    # SYNOPSIS:                                                               #
    #  logpdf = loggausspdf(X,mu,Sigma)                                       #
    #                                                                         #
    # © Written by Michael Chen (sth4nth@gmail.com)                           #
    #   Modified by Jasper A. Vrugt                                           #
    #                                                                         #
    # ####################################################################### #

    d = X.shape[0]                      # Number of features
    diff = X - mu[:, np.newaxis]        # Shape (d, n)
    inv_Sigma = np.linalg.inv(Sigma)    # Inverse of the covariance matrix
    norm_factor = 0.5 * (d * np.log(2 * np.pi) + np.log(np.linalg.det(Sigma)))
    exponent = -0.5 * np.sum(diff.T @ inv_Sigma * diff.T, axis=1)  # Quadratic form

    return -norm_factor + exponent      # Log of the Gaussian PDF


def logmeanexp(x, dim = 0):
    # ####################################################################### #
    #LOGMEANEXP Returns log(mean(exp(a),dim)) while avoiding numerical        #
    # underflow                                                               #
    #                                                                         #
    # dim = 1 for taking the mean along each column of a (default)            #
    # dim = 2 for taking the mean along each row of a                         #
    #                                                                         #
    # © Modified by Jasper A. Vrugt, Aug. 2015                                #
    # University of California Irvine                                         #
    #                                                                         #
    # ####################################################################### #

    # Use logsumexp to avoid numerical underflow
    logsum_exp = logsumexp(x, dim)
    s = logsum_exp - np.log(x.shape[dim])
    
    return s


def logsumexp(x, dim = None):
    # ####################################################################### #
    #LOGSUMEXP Returns log(sum(exp(x),dim)) while avoiding numerical          #
    # underflow                                                               #
    #                                                                         #
    # dim = 1 for taking the mean along each column of x (default)            #
    # dim = 2 for taking the mean along each row of x                         #
    #                                                                         #
    # © Written by Michael Chen (sth4nth@gmail.com)                           #
    #                                                                         #
    # ####################################################################### #

    if dim is None:
        # Determine the dimension along which to sum (default behavior)
        dim = np.argmax(np.array(x.shape) != 1)
        if dim is None:
            dim = 0             # Summing along columns (= default)
    else:
        dim = int(dim)

    # Subtract the largest value along the specified dimension to prevent overflow
    y = np.max(x, axis = dim, keepdims = True)
    x = x - y                   # Subtract the max

    # Calculate the sum of the exponentials and then take the log
    s = y + np.log(np.sum(np.exp(x), axis = dim, keepdims = True))

    # Handle the case where the values in y are not finite
    s[np.isnan(y) | np.isinf(y)] = y[np.isnan(y) | np.isinf(y)]

    return s.squeeze()  # Return the result with the same shape as the input (squeeze the singleton dimensions)


def qsimvn(m, r, a, b, set_generator=2):
    # ####################################################################### #
    #QSIMVN Uses a randomized quasi-random rule with m points to estimate an  #
    # multivariate normal probability distribution for positive definite      #
    # covariance matrix r, with lower integration limits a and upper          #
    # integration limits b. Probability p is output with error estimate e     #
    #                                                                         #
    #  Example usage:                                                         #
    #     >> r = [4 3 2 1;3 5 -1 1;2 -1 4 2;1 1 2 5];                         #
    #     >> a = -inf*[1 1 1 1 ]'; b = [ 1 2 3 4 ]';                          #
    #     >> [p e] = qsimvn(5000,r,a,b); disp([p e])                          #
    #                                                                         #
    #  This function uses an algorithm given in the paper                     #
    #   Genz, Alan (1992), Numerical Computation of Multivariate Normal       #
    #       Probabilities, Journal of Computational and Graphical Statistics, #
    #       1, pp. 141-149. WSU Math, PO Box 643113, Pullman, WA 99164-3113   #
    #                                                                         #
    #  The primary references for the numerical integration are               #
    #   Niederreiter, H. (1972), On a number-theoretical integration method,  #
    #       Aequationes Mathematicae, 8, pp. 304-11                           #
    #   Cranley, R. and T.N.L. Patterson (1976), Randomization of number      #
    #       theoretic methods for multiple integration, SIAM Journal of       #
    #       Numerical Analysis, 13, pp. 904-14                                #
    #                                                                         #
    # Copyright (C) 2013, Alan Genz,  All rights reserved                     #
    #                                                                         #
    # Redistribution and use in source and binary forms, with or without      #
    # modification, are permitted provided the following conditions are met:  #
    #   1. Redistributions of source code must retain the above copyright     #
    #      notice, this list of conditions and the following disclaimer       #
    #   2. Redistributions in binary form must reproduce the above copyright  #
    #      notice, this list of conditions and the following disclaimer in    #
    #      the documentation and/or other materials provided with the         #
    #      distribution                                                       #
    #   3. The contributor name(s) may not be used to endorse or promote      #
    #      products derived from this software without specific prior         #
    #      written permission                                                 #
    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     #
    # "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT       #
    # LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS       #
    # FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE          #
    # COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,     #
    # INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,    #
    # BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS   #
    # OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  #
    # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR   #
    # TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE  #
    # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE     #
    #                                                                         #
    # Alan Genz is the author of this function and following Matlab functions #
    # Email: AlanGenz@wsu.edu                                                 #
    # Location: WSU Math, PO Box 643113, Pullman, WA 99164-3113               #
    #                                                                         #
    # ####################################################################### #

    n = r.shape[1]
    ch, as_, bs_ = chlrdr(r, a, b)          # Code benchmarked against MATLAB
    ct = ch[0, 0]
    ai = as_[0]
    bi = bs_[0]
    
    if np.abs(ai) < 9 * ct:
        c = norm.cdf(ai / ct)
    else:
        c = (1 + np.sign(ai)) / 2
    if np.abs(bi) < 9 * ct:
        d = norm.cdf(bi / ct)
    else:
        d = (1 + np.sign(bi)) / 2
    
    ci = c
    dci = d - ci
    p = 0
    e = 0

    ns = int(12)
    nv = max([m // ns, 1])

    # Set generator choice for quasi-random points [benchmarked]
    if set_generator == 1:
        q = 2 ** (np.arange(1, n) / n)
    elif set_generator == 2:
        # ps = np.sqrt(np.array([p for p in primes(int(5 * n * np.log(n + 1) / 4))]))       # Generate primes (adjusted for Python)
        # Calculate the upper limit for primes
        primes_up_to = int(1 + 5 * n * np.log(n + 1) / 4)
        primes_list = list(primerange(1, primes_up_to))
        ps = np.sqrt(np.array(primes_list))
        # original code
        q = ps[:n-1]

    # Randomization loop for ns samples
    for i in range(0, ns):
        vi = 0
        xr = np.random.rand(n - 1)
        for j in range(0, nv):
            x = np.abs(2 * np.mod((j+1) * q + xr, 1) - (1))                               # Periodizing transformation
            vp = mvndns(n, ch, ci, dci, x, as_, bs_)
            vi = vi + (vp - vi) / (j + 1)
        
        d = (vi - p) / (i + 1)
        p = p + d
        
        if np.abs(d) > 0:
            e = np.abs(d) * np.sqrt(1 + (e / d) ** 2 * (i - 1) / (i + 1))
        else: 
            if i > 0:
                e = e * np.sqrt((i - 1) / (i + 1))

    e = 3 * e  # Error estimate is 3 x standard error with ns samples

    return p, e


def mvndns(n, ch, ci, dci, x, a, b):
    # ####################################################################### #
    #MVNDNS Transformed integrand for computation of MVN probabilities        #
    # ####################################################################### #

    y = np.zeros(n - 1)
    s = 0
    c = ci
    dc = dci
    p = dc
    
    for i in range(1, n):
        y[i - 1] = norm.ppf(c + x[i - 1] * dc)
        # s = ch(i,1:i-1)*y(1:i-1);
        # s = np.dot(ch[i, 0:i-1], y[0:i-1])
        # Python: i_max loop = n-1 --> then 0:n-2 is length y of n-1  
        s = np.dot(ch[i, np.arange(0,i)], y[np.arange(0,i)])
        ct = ch[i, i]
        ai = a[i] - s
        bi = b[i] - s
        if np.abs(ai) < 9 * ct:
            c = norm.cdf(ai / ct)
        else:
            c = (1 + np.sign(ai)) / 2
        if np.abs(bi) < 9 * ct:
            d = norm.cdf(bi / ct)
        else:
            d = (1 + np.sign(bi)) / 2
        dc = d - c
        p *= dc
    
    return p


def chlrdr(R, a, b):
    # ####################################################################### #
    #CHLRDR Computes the permuted lower Cholesky factor c for R, which may be #
    # singular, and also permutes the integration limit vectors a and b.      #
    # ####################################################################### #

    ep = 1e-10              # Singularity tolerance
    n = R.shape[1]
    c = np.copy(R)
    ap = np.copy(a)
    bp = np.copy(b)
    d = np.sqrt(np.maximum(np.diag(c), 0))
    for i in range(0, n):
        if d[i] > 0:
            c[:, i] /= d[i]
            c[i, :] /= d[i]
            ap[i] /= d[i]
            bp[i] /= d[i]

    y = np.zeros(n) 
    sqtp = np.sqrt(2 * np.pi)
    for k in range(0, n):
        im = k
        ckk = 0
        dem = 1
        s = 0
        for i in range(k, n):
            if c[i, i] > np.finfo(float).eps:
                cii = np.sqrt(max([c[i, i], 0]))
                if i > 0:
                    # s = np.dot(c[i, 0:k - 1], y[0:k - 1])
                    s = c[i, 0:k - 1] @ y[0:k - 1]
                ai = (ap[i] - s) / cii
                bi = (bp[i] - s) / cii
                de = norm.cdf(bi) - norm.cdf(ai)
                if de <= dem:
                    ckk = cii.copy()
                    dem = de
                    am = ai
                    bm = bi
                    im = i

        if im > k:
            # Swap values between ap(im) and ap(k)
            tv = ap[im]
            ap[im] = ap[k]
            ap[k] = tv
            # Swap values between bp(im) and bp(k)
            tv = bp[im]
            bp[im] = bp[k]
            bp[k] = tv
            # Swap c(im, im) and c(k, k)
            c[im, im] = c[k, k].copy()
            # Swap rows and columns for the submatrix in c
            t = c[im, np.arange(0,k)].copy()
            #c[im, 0:k-1] = c[k, 0:k-1]
            c[im, np.arange(0,k)] = c[k, np.arange(0,k)]
            #c[k, 0:k-1] = t
            c[k, np.arange(0,k)] = t

            t = c[im+1:n, im].copy()
            c[im+1:n, im] = c[im+1:n, k]
            c[im+1:n, k] = t
            #t = c[k+1:im-1, k]
            t = c[k+1:im, k].copy()
            #c[k+1:im-1, k] = c[im, k+1:im-1].T
            c[k+1:im, k] = c[im, k+1:im].T
            #c[im, k+1:im-1] = t.T
            c[im, k+1:im] = t.T

        if ckk > (ep * (k+1)):
            c[k, k] = ckk
            c[k, np.arange(k+1,n)] = (0)
            for i in range(k + 1, n):
                c[i, k] = c[i, k] / ckk
                # c[i, k+1:i] = c[i, k+1:i] - c[i, k] * c[k+1:i, k].T
                c[i, np.arange(k + 1,i+1)] = c[i, np.arange(k + 1,i+1)] - c[i, k] * c[np.arange(k+1,i+1), k].T

            if np.abs(dem) > ep:
                y[k] = (np.exp(-am ** 2 / 2) - np.exp(-bm ** 2 / 2)) / (sqtp * dem)
            else:
                if am < -10:
                    y[k] = bm
                elif bm > 10:
                    y[k] = am
                else:
                    y[k] = (am + bm) / 2
        else:
            c[k:n, k] = 0
            y[k] = 0
    
    return c, ap, bp
