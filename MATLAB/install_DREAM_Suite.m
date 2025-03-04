% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%             DDDDD     RRRRR     EEEEEEE     AAA     MM    MM            %
%             DDDDDD    RRRRRR    EEEEEEE    AAAAA    MM    MM            %
%             DD  DD    RR   RR   EE        AA   AA   MMM  MMM            %
%             DD   DD   RR  RR    EEEE      AA   AA   MMMMMMMM            %
%             DD   DD   RRRRR     EEEE      AAAAAAA   MMM  MMM            %
%             DD  DD    RR RR     EE        AAAAAAA   MM    MM            %
%             DDDDDD    RR  RR    EEEEEEE   AA   AA   MM    MM            %
%             DDDDD     RR   RR   EEEEEEE   AA   AA   MM    MM            %
%                                                                         %
%              SSSSSSSS  UU    UU   II   TTTTTTTTTT   EEEEEEE             %
%              SSSSSSS   UU    UU   II   TTTTTTTTTT   EEEEEEE             %
%              SS        UU    UU   II       TT       EE                  %
%              SSSS      UU    UU   II       TT       EEEE                %
%                 SSSS   UU    UU   II       TT       EEEE                %
%                   SS   UU    UU   II       TT       EE                  %
%               SSSSSS   UUUUUUUU   II       TT       EEEEEEE             %
%              SSSSSSS   UUUUUUUU   II       TT       EEEEEEE             %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% DREAM-Suite: DiffeRential Evolution Adaptive Metropolis algorithm       %
% with continuous and/or discrete variables and single- or multi-try      %
% sampling from an archive of current or past states using parallel       %
% direction, snooker and/or Kalman candidate points. This toolbox         %
% implements in one package the DREAM, DREAM_D, DREAM_ZS, DREAM_ZS,       %
% MTDREAM_ZS and DREAM_KZS algorithms.                                    %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% SYNOPSIS:                                                               %
%  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          %
%      DREAMPar,Par_info)                                                 %
%  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          %
%      DREAMPar,Par_info,Meas_info)                                       %
%  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          %
%      DREAMPar,Par_info,Meas_info,options)                               %
%  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          %
%      DREAMPar,Par_info,Meas_info,options,MAP_info)                      %
%  [chain,output,FX,Z,loglik] = DREAM_Suite(method,Func_name,...          %
%      DREAMPar,Par_info,Meas_info,options,MAP_info,plugin)               %
% WHERE                                                                   %
%  method      [input] Name (string) of MCMC method                       %
%   = 'dream'                                                             %
%   = 'dream_zs'                                                          %
%   = 'dream_d'                                                           %
%   = 'dream_dzs'                                                         %
%  Func_name   [input] Function (string) returns (log)lik or sim values   %
%   → Func_name must return a likelihood if DREAMPar.lik = 1              %
%   → Func_name must return a log-likelihood if DREAMPar.lik = 2          %
%   → Func_name returns vector of n simulated values: DREAMPar.lik > 2    %
%  DREAMPar    [input] Structure with algorithmic variables               %
%   .d             Dimensionality (# variables) target distribution       %
%   .N             # of Markov chains                                     %
%   .T             # of generations (= # samples of each Markov chain)    %
%   .lik           Choice of likelihood function                          %
%     = 1          User computes likelihood in Func_name                  %
%     = 2          User computes log-likelihood in Func_name              %
%     = 11         Normal likelihood (= post. density) [Box and Tiao]     %
%     = 12         Normal likelihood & (non)const err Meas_info.Sigma     %
%     = 13         Normal likelihood & (non)const var. + AR(2) process    %
%     = 14         Genrlzd likelihood Schoups & Vrugt (2010) [OBSOLETE]   %
%     = 15         Whittle likelihood function (see Whittle, 1953)        %
%     = 16         Laplace likelihood & (non)const var. + AR(1) proc.     %
%     = 17         Student t likelihood & (non)const var. + AR(2) proc.   %
%     = 21         Approximate Bayesian Computation: normal likelihood    %
%     = 22         Approximate Bayesian Computation: boxcar function      %
%     = 23         Limits of Acceptability: log-liklhood is # pts LOA     %
%     = 31         GLUE: log-likelihood: Table 1a Beven & Freer 2001      %
%     = 32         GLUE: log-likelihood: Table 1b Beven & Freer 2001      %
%     = 33         GLUE: log-likelihood: Table 1c Beven & Freer 2001      %
%     = 34         GLUE: log-likelihood: Page 284 Beven & Binley 1992     %
%     = 44         Generalized likelihood ++: Vrugt et al. 2022           %
%     = 45         Universal likelihood: Vrugt et al. 2022                %
%     = 52         Log-likelihood: Generalized Least Squares (GLS) form   %
%     = 61         Laplace power likelihood: unit intgrl, lambda estmtd   %
%     = 62         Normal power likelihood: unit intgrl, lambda estmtd    %
%   .nCR           # crossover values                    DEF: 3           %
%   .delta         # chain pairs for proposal            DEF: 3           %
%   .lambda        Random error for ergodicity           DEF: 0.05        %
%   .zeta          Randomization                         DEF: 0.05        %
%   .p_unit_gamma  Probability unit jumprate (gamma)     DEF: 0.2         %
%   .adapt_pCR     Adapt crossover probabilities?        DEF: 'yes'       %
%   .thinning      Each thinning(th) chain sample stored DEF: 1           %
%   .GLUE          GLUE likelihood parameter             DEF: 10          %
%   .beta0         Scaling factor built-in jump rate     DEF: 1           %
%   .outlier       Outlier chain detection test          DEF: 'iqr'       %
%                   → DREAM and DREAM_D                                   %
%   .psnooker      Selection probability snooker jump    DEF: 0.1         %
%                   → DREAM_ZS, DREAM_DZS and MTDREAM_ZS                  %
%   .m0            Initial size of external archive, Z   DEF: 10*d        %
%                   → DREAM_ZS, DREAM_DZS and MTDREAM_ZS                  %
%   .k             Growth rate of external archive       DEF: 10          %
%                   → DREAM_ZS, DREAM_DZS and MTDREAM_ZS                  %
%   .M             # samples archive Z for Kalman jump   DEF: 20          %
%                   → DREAM_KZS                                           %
%   .mt            Number of multi-try proposals         DEF: 5           %
%                   → MTDREAM_ZS                                          %
%  Par_info    [input] Parameter structure: Ranges, initial/prior & bnd   %
%   .names         1xd-cell array with parameter names   DEF: []          %
%   .min           1xd-vector of min parameter values    DEF: -inf(1,d)   %
%   .max           1xd-vector of max parameter values    DEF: inf(1,d)    %
%   .norm          Work in normlzed parameter space      DEF: 0           %
%     = 0          Work in unnormlzed parameter space    DEFault          %
%     = 1          Work in normalized [0-1] parameter space               %
%   .boundhandling Treat the parameter bounds or not?                     %
%     = 'reflect'  Reflection method                                      %
%     = 'bound'    Set to bound                                           %
%     = 'fold'     Folding [Vrugt&Braak: doi:10.5194/hess-15-3701-2011]   %
%     = 'reject'   Reject out of bound proposals                          %
%     = 'none'     No boundary handling                  DEFault          %
%   .initial       Method to draw initial chain states                    %
%     = 'uniform'  Uniform: U(Par_info.min,Par_info.max)                  %
%     = 'latin'    Latin hypercube: LH(Par_info.min,Par_info.max)         %
%     = 'normal'   Normal:  N(Par_info.mu,Par_info.cov)                   %
%     = 'prior'    User specified prior distribution                      %
%     = 'user'     Initial chain states taken from Par_info.x0            %
%   .mu            1xd-mean vector: µ if .initial = 'normal'              %
%   .cov           dxd-covariance matrix: Σ if .initial = 'normal'        %
%   .x0            N x d matrix of initial states if .initial = 'user'    %
%   .prior         Prior distribution (manual) if .initial = 'prior'      %
%                  Ex 1: Par_info.prior = @(x,a,b) mvnpdf(x,a,b);         %
%                            Par_info.a = [-2 -2];                        %
%                            Par_info.b = eye(2);                         %
%                  Ex 2: Par_info.prior = {'normpdf(x,0,1)',...           %
%                                          'unifpdf(x,-2,2)'}             %
%                  Note: prior handle can return log(pdf): Code checks    %
%   .steps         d-vector with # intervals for each parameter           %
%                   → DREAM_D/DREAM_DZS                                   %
%  Meas_info   [input] Structure with measurement information (fitting)   %
%   .Y             nx1 vector against which model output is compared      %
%   .Sigma         Measurement error standard deviation of Y              %
%     = scalar     Constant measurement error                             %
%     = vector     Nonconstant measurement error (nx1 vector)             %
%   .sigma2        Treatment meas. error var. likelihood 13/16/17/44/45   %
%     = 'constant      Homoscedastic measurement error variance           %
%     = 'nonconstant'  Heteroscedastic measurement error variance         %
%   .S             Scalar/vector with summary metrics                     %
%   .C             dxd meas. err. cov. matrix GLS likelihood [= 52]       %
%   .R             Measurement error covariance matrix                    %
%                   → DREAM_KZS                                           %
%  options     [input] Structure with computational settings/options      %
%   .parallel      Multi-core computation chains?        DEF: 'yes'       %
%   .IO            If parallel, IO writing model?        DEF: 'no'        %
%   .rho           ABC dist. func. (anonymous handle)  DEF: @(X,Y)|X-Y|   %
%   .epsilon       ABC epsilon value (scalar/vector)     DEF: 0.025       %
%   .DB            Diagnostic Bayes?                     DEF: 'no'        %
%   .modout        Return model simulations?             DEF: 'no'        %
%   .save          Save DREAM output during the run?     DEF: 'no'        %
%   .restart       Restart run? (only with "save")       DEF: 'no'        %
%   .diagnostics   Compute within-chain diagnostics?     DEF: 'yes'       %
%   .print         Output writing screen (tables/figs)   DEF: 'yes'       %
%   .burnin        Burn-in % chain for conv. diagnstcs   DEF: 50          %
%  MAP_info    [input] Structure with information about MAP solution      %
%   .map           1xd vector with MAP solution                           %
%   .An            dxd sensitivity matrix MAP solution                    %
%   .Betan         dxd variability matrix MAP solution = Bn if unbiased   %
%  plugin      [input] 2nd input argument Func_name. Class set by user    %
%                                                                         %
%  chain       [outpt] T x (d+2) x N array N chains: pars+logpr+loglik    %
%                      If thinning is used then length ~ T/thinning       %
%  output      [outpt] Structure summarizes algorithmic performance       %
%   .R_stat        Univariate \hat{R} convergence diagnostic              %
%   .MR_stat       Multivariate \hat{R} convergence diagnostic            %
%   .AR            Acceptance rate (%)                                    %
%   .CR            Crossover selection probabilities                      %
%   .outlier       Information about outlier chains                       %
%   .RunTime       CPU time in seconds                                    %
%  FX          [outpt] T/thinning x n matrix model simuls chain samples   %
%  Z           [outpt] External archive used by DREAM_ZS and DREAM_DZS    %
%  loglik      [outpt] Log-likelihood sampled chains [diagnostics only]   %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %
%                                                                         %
% DREAM PACKAGE IS DESCRIBED IN                                           %
%   Vrugt, J.A., R. de Punder, and P. Grünwald, A sandwich with water:    %
%       Bayesian/Frequentist uncertainty quantification under model       %
%       misspecification, Submitted to Water Resources Research,          %
%       May 2024, https://essopenarchive.org/users/597576/articles/...    %
%           937008-a-sandwich-with-water-bayesian-frequentist-...         %
%           uncertainty-quantification-under-model-misspecification       %
%   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and           %
%       Diagnostics: Elicitability, Propriety, and Scoring Rules for      %
%       Hydrograph Functionals, Water Resources Research, 60,             %
%       e2023WR036710, https://doi.org/10.1029/2023WR036710               %
%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022),    %
%       On the use of distribution-adaptive likelihood functions:         %
%       Generalized and universal likelihood functions, scoring rules     %
%       and multi-criteria ranking, Journal of Hydrology, 615, Part B,    %
%       2022, https://doi.org/10.1016/j.jhydrol.2022.128542               %
%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the     %
%       DREAM software package: Theory, concepts, and MATLAB              %
%       implementation, Environmental Modeling and Software, 75,          %
%       pp. 273-316, https://doi.org/10.1016/j.envsoft.2015.08.013        %
%   Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation   %
%       using Markov chain Monte Carlo simulation: DREAM_(ABC), Water     %
%       Resources Research, https://doi.org/10.1002/2014WR015386          %
%   Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model            %
%       calibration and evaluation: Approximate Bayesian computation,     %
%       Water Resources Research, 49, pp. 4335–4345,                      %
%           https://doi.org/10.1002/wrcr.20354                            %
%   Laloy, E., and J.A. Vrugt (2012), High-dimensional posterior          %
%       exploration of hydrologic models using multiple-try DREAM_(ZS)    %
%       and high-performance computing, Water Resources Research, 48,     %
%       W01526, https://doi.org/10.1029/2011WR010608                      %
%   Vrugt, J.A., and C.J.F. ter Braak (2011), DREAM_(D): An adaptive      %
%       Markov chain Monte Carlo simulation algorithm to solve            %
%       discrete, noncontinuous, and combinatorial posterior parameter    %
%       estimation problems, Hydrology and Earth System Sciences, 15,     %
%       pp. 3701-3713, https://doi.org/10.5194/hess-15-3701-2011          %
%   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and                        %
%       B.A. Robinson (2009), Equifinality of formal (DREAM) and          %
%       informal (GLUE) Bayesian approaches in                            %
%       hydrologic modeling? Stochastic Environmental Research and Risk   %
%       Assessment, 23(7), pp. 1011-1026,                                 %
%           https://doi.org/10.1007/s00477-008-0274-y                     %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon,                %
%       B.A. Robinson, and J.M. Hyman (2009), Accelerating Markov chain   %
%       Monte Carlo simulation by differential evolution with             %
%       self-adaptive randomized subspace sampling, International         %
%       Journal of Nonlinear Sciences and Numerical Simulation, 10(3),    %
%       pp. 271-288                                                       %
%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and            %
%       B.A. Robinson (2008), Treatment of input uncertainty in           %
%       hydrologic modeling: Doing hydrology backward with Markov chain   %
%       Monte Carlo simulation, Water Resources Research, 44, W00B09,     %
%       https://doi.org/10.1029/2007WR006720                              %
%   Ter Braak, C.J.F., and J.A. Vrugt (2008), Differential Evolution      %
%       Markov Chain with snooker updater and fewer chains, Statistics    %
%       and Computing, https://doi.org/10.1007/s11222-008-9104-9          %
%   Ter Braak, C.J.F. (2006), A Markov Chain Monte Carlo version of the   %
%       genetic algorithm differential evolution: easy Bayesian           %
%       computing for real parameter spaces, Statistics and Computing,    %
%       16, pp. 239-249, doi:10.1007/s11222-006-8769-1                    %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% DIFFERENT TEST EXAMPLES                                                 %
%  example 1: d-dimensional banana shaped Gaussian distribution           %
%  example 2: d-dimensional Gaussian distribution                         %
%  example 3: d-dimensional multimodal normal mixture distribution        %
%  example 4: real-world example rainfall-runoff (hymod in C++/MATLAB)    %
%  example 5: rainfall-runoff (hymod as external executable)              %
%  example 6: hmodel with distribution-adaptive likelihood functions      %
%  example 7: HYDRUS-1D soil hydraulic model: multiplicative prior        %
%  example 8: Approximate Bayesian Computation: Benchmark function        %
%  example 9: Spectral likelihood function in watershed modeling          %
%  example 10: Gaussian mixture distibution: multivariate prior           %
%  example 11: d-variate t-distribution: df ° freedom & corr. matrix R    %
%  example 12: pedometrics problem involving variogram fitting            %
%  example 13: Nash-Cascade hydrograph                                    %
%  example 14: Approx. Bayesian Comp. watershed signatures                %
%  example 15: Approx. Bayesian Comp. bivariate normal benchmark test     %
%  example 16: Hydrogeophysical inversion                                 %
%  example 17: Watershed model, normal, AR(1) and heteroscedastic lik.    %
%  example 18: Lotka-Volterra model: informal likelihood (GLUE)           %
%  example 19: Bayesian Model Averaging: I recommend MODELAVG toolbox!    %
%  example 20: Limits of acceptability: Soil temperature modeling         %
%  example 21: Limits of acceptability: Soil moisture model HYDRUS-1D     %
%  example 22: Limits of acceptability: Nash-Cascade hydrograph           %
%  example 23: Limits of acceptability: SAC-SMA (old C-code Euler int.)   %
%  example 24: Flow duration curve fitting                                %
%  example 25: Bedrock depth from high-res topo data & geomorph model     %
%  example 26: Data assimilation Lorenz model (SODA: Vrugt et al. 2005)   %
%  example 27: Data assimilation interception model (Vrugt et al. 2003)   %
%  example 28: Rainfall and hmodel parameter estimation from streamflow   %
%  example 29: Gaussian mixture distribution & multiplicative prior       %
%  example 30: Predator prey interactions                                 %
%  example 31: AR(2)-parameter estimation: Test of likelihood functions   %
%  example 32: Distribution-adaptive likelihood functions                 %
%  example 33: 2-dimensional rectangular target distribution              %
%  example 34: Haverkamp infiltration equation using HYDRUS-1D data       %
%  example 35: Haverkamp infiltration equation using SWIG database        %
%  example 36: Sediment transport modeling                                %
%  example 99: Bayesian inference & M-estimation: Sandwich correction     %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory                                          %
%                                                                         %
% FURTHER CHECKING                                                        %
%  Website:  http://faculty.sites.uci.edu/jasper                          %
%  Papers: http://faculty.sites.uci.edu/jasper/publications/              %
%  Google Scholar: https://scholar.google.com/citations?user=...          %
%                      zkNXecUAAAAJ&hl=nl                                 %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% HOW TO INSTALL TOOLBOX AND EXECUTE EXAMPLES                             %
%  Step 1: Go to main DREAM-Suite directory                               %
%  Step 2: Type: install_DREAM_Suite                                      %
%  Step 3: Go to an example directory, type: cd example_1                 %
%  Step 4: Execute this example in the command prompt: example_1          %
%    Visual results are printed to screen unless options.print = 'no'     %
%    You can create a 2d matrix from the 3d-chain array by typing:        %
%        parset = genparset(chain);                                       %
%    This gives you all sampled parameter vectors by the chains. You      %
%    will need to apply burn-in to obtain the posterior samples           %
%                                                                         %
% HOW TO USE GAME SAMPLING                                                %
%  Step 1: Unpack all chain samples:                                      %
%              P = genparset(chain);                                      %
%  Step 2: Burn-in for posterior samples, say:                            %
%              Post = P(end-9999:end,1:DREAMPar.d+2);                     %
%  Step 3: Now calculate marginal likelihood                              %
%              Z = GAME_sampling(Post,method,DREAMPar,Func_name);         %
%          where "method" is name (string) marginal likelihood estmtr     %
%                     = 'ris'   reciprocal importance sampling            %
%                     = 'is'    importance sampling                       %
%                     = 'ob'    optimal bridge sampling                   %
%                     = 'gb'    geometric bridge sampling                 %
%  Type "edit readme.txt" MATLAB prompt for more about GAME_sampling      %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% In root directory of eDREAM PACKAGE 
addpath(pwd,[pwd '/postprocessing'],[pwd '/diagnostics'], ...
    [pwd '/gamesampling'], [pwd '/gamesampling/miscellaneous'], ...
    [pwd '/export'], [pwd '/num2words'], [pwd '/miscellaneous']);
% Now go to example directory; say example_1
cd example_1
% Now execute this example by typing in command prompt: "example_1" 

% After DREAM-Suite terminates, create a 2d matrix from 3D-chain array
% P = genparset(chain);
% This matrix will have all samples of the joint chains (thinning active)
% Burn-in of P is required to get posterior samples

% Figures should be printed to the screen and text/binary files written