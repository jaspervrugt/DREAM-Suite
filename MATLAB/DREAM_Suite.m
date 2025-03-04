function varargout = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,varargin)
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
%   = 'mtdream_zs'                                                        %
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
%   .mt            Number of multi-try proposals         DEF: 5           %
%                   → MTDREAM_ZS                                          %
%   .M             # samples archive Z for Kalman jump   DEF: 20          %
%                   → DREAM_KZS                                           %
%   .a_1           # Kalman jump begins at a_1 *.T gens  DEF: 0.1         %
%                   → DREAM_KZS                                           %
%   .a_2           # Kalman jump ends at a_2 *.T gens    DEF: 0.25        %
%                   → DREAM_KZS                                           %
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
% The different components/algorithms of DREAM Package are described in   %
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
%  example 99: Bayesian inference & M-estimation: Sandwich correction     %
%                                                                         %
% © Written by Jasper A. Vrugt, Feb 2007                                  %
% Los Alamos National Laboratory                                          %
% University of California Irvine                                         %
%                                                                         %
% FURTHER CHECKING                                                        %
%  Website:  http://faculty.sites.uci.edu/jasper                          %
%  Papers: http://faculty.sites.uci.edu/jasper/publications/              %
%  Google Scholar: https://scholar.google.com/citations?user=...          %
%                      zkNXecUAAAAJ&hl=nl                                 %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% <><><><><><><><><><><><>< CHECK INPUT ARGS ><><><><><><><><><><><><><><>
plugin = []; [MAP_info,options] = deal(struct); method = lower(method);
if nargin == 8, [Meas_info,options,MAP_info,plugin] = varargin{1:4}; end
if nargin == 7, [Meas_info,options,MAP_info] = varargin{1:3}; end
if nargin == 6, [Meas_info,options] = varargin{1:2}; end
if nargin == 5, Meas_info = varargin{1}; end
if nargin < 5 || isempty(Meas_info), Meas_info = struct; end
if nargin < 4, error(['DREAM_Suite ERROR:TooFewInputs: Requires at '...
        'least four input arguments.']); end
if ~isfield(options,'restart'), options.restart = 'no'; end
tot_accept = 0; start_time = cputime; count = 0;
% <><><><><><><><><><><>< END CHECK INPUT ARGS ><><><><><><><><><><><><><>

% <><><><><><><><><><><><><> INITIALIZATION <><><><><><><><><><><><><><><>
switch options.restart
    case 'no' % Start from scratch: check algorithmic variables
        [DREAMPar,Par_info,options] = DREAM_Suite_check(method,...
            Func_name,DREAMPar,Par_info,Meas_info,options);
        % Setup of the main variables used in DREAM
        [DREAMPar,Par_info,Meas_info,Lik_info,options] = ...
            DREAM_Suite_setup(method,Func_name,DREAMPar,Par_info,...
            Meas_info,options);
        % Computing environment (multi-core?)
        [DREAMPar,f_handle,T] = DREAM_Suite_calc_setup(DREAMPar,...
            Func_name,options,plugin);
        % Sample initial chain states (initial population) and so forth
        [chain,output,X,FX,S,Table_gamma,CR,pCR,cCR,TdCR,sdX_Xp,loglik,...
            EX,std_mX,id_r,id_c,iloc,it,g,LLX,MAP_info] = ...
            DREAM_Suite_initialize(method,DREAMPar,Par_info,...
            Meas_info,Lik_info,options,f_handle,MAP_info);
    case 'yes' % Restart: Load variables from last run
        [DREAMPar,Par_info,Meas_info,Lik_info,options,f_handle,...
            MAP_info,chain,output,X,FX,S,Table_gamma,CR,pCR,cCR,...
            TdCR,sdX_Xp,loglik,EX,std_mX,id_r,id_c,iloc,it,g,LLX,T] = ...
            DREAM_Suite_restart('DREAM_Suite.mat');
end
% <><><><><><><><><><><><> END INITIALIZATION <><><><><><><><><><><><><><>
N = DREAMPar.N * DREAMPar.mt;
save temp.mat
% <><><><><><><><><><><><><><> DYNAMIC PART <><><><><><><><><><><><><><><>
for t = T : DREAMPar.T
    % Save current population to be X_old
    X_old = X(1:DREAMPar.N,1:DREAMPar.d);
    switch method
        case 'mtdream_zs'
            % Determine whether to do parallel direction or snooker update
            if rand <= DREAMPar.psnooker, j_mthd = 2; else, j_mthd = 1; end
            % Create DREAMPar.mt proposals in each chain
            [Xp,log_alfa_sn_Xp,vp,Z,CR_par,gamma] = Calc_MTproposal( ...
                X(1:DREAMPar.N,1:DREAMPar.d),CR(:,g),DREAMPar, ...
                Table_gamma,Par_info,j_mthd,1); 
        otherwise
            % Create candidate points in each chain
            [Xp,log_alfa_sn,CR(:,g),vp] = Calc_proposal(method,...
                X(1:DREAMPar.N,1:DREAMPar.d),EX,std_mX,CR(:,g),...
                DREAMPar,Table_gamma,Par_info,Meas_info);
    end
    % Compute unnormalized parameter values
    Xp_un = X_unnormalize(Xp(:,1:DREAMPar.d),Par_info);
    % Now evaluate the model ( = pdf ) and return fx
    [FXp,SXp] = Evaluate_target(Xp_un,DREAMPar,Meas_info,options,f_handle);
    FXp
    pause
    % Compute the log(prior) of candidate points
    logPR_Xp = Eval_prior(Xp_un,SXp,Par_info,Meas_info,options);
    % Compute the log-likelihood of candidate points
    [logL_Xp,EXp,std_mXp,~,~,LLp] = Calc_likelihood(Xp_un,FXp,DREAMPar,...
        Par_info,Meas_info,Lik_info,options,MAP_info);
    % Append log(prior) and log-likelihood to Nxd matrix Xp
    Xp(1:N,DREAMPar.d+1:DREAMPar.d+2) = [logPR_Xp logL_Xp];
    % prop. loglik if range violation: Par_info.boundhandling = 'reject'
    % then, Xp(vp==1,:) [= Xp(vp,:)] will be declined
    Xp(vp,DREAMPar.d+2) = -inf;

    switch method
        case 'mtdream_zs'
            % Sample proposal chain with probability proportional to weight
            [log_wXp,id_pop] = select_Xp(DREAMPar,Xp(1:N, ...
                DREAMPar.d+1:DREAMPar.d+2),log_alfa_sn_Xp);
            % Make sure that reference set uses same crossover values as Xp
            CR(:,g) = CR_par(id_pop);
            % Sample reference set, Xr, centered on Xp [= pref cand. points]
            [Xr,log_alfa_sn_Xr,vr] = Calc_MTproposal(Xp(id_pop, ...
                1:DREAMPar.d),CR(:,g),DREAMPar,Table_gamma,Par_info, ...
                j_mthd,2,gamma,Z);
            % Compute unnormalized reference parameter values
            Xr_un = X_unnormalize(Xr(:,1:DREAMPar.d),Par_info);
            % Now evaluate the model ( = pdf ) and return fx
            [FXr,SXr] = Evaluate_target(Xr_un(id_r,1:DREAMPar.d), ...
                DREAMPar,Meas_info,options,f_handle);
            % Compute the log(prior) of reference points
            % Must fix with SXr as this is not an empty cell
            logPR_Xr = Eval_prior(Xr_un(id_r,1:DREAMPar.d),SXr, ...
                Par_info,Meas_info,options);
            % Compute the log-likelihood of the reference points
            [logL_Xr,~,~,~,~,LLr] = Calc_likelihood(Xr_un(id_r, ...
                1:DREAMPar.d),FXr,DREAMPar,Par_info,Meas_info, ...
                Lik_info,options);                                  %#ok
            % Append log(prior) and log-likelihood to Nxd matrix Xr
            Xr(id_r,DREAMPar.d+1:DREAMPar.d+2) = [logPR_Xr logL_Xr];
            % ref loglik range violation: Par_info.boundhandling = 'reject'
            % then, Xr(vr==1,:) [= Xr(vr,:)] will be declined
            Xr(id_r(vr),DREAMPar.d+2) = -inf;
            % Augment Xr with current position of each chain
            Xr(id_c,1:DREAMPar.d+2) = X(1:DREAMPar.N,1:DREAMPar.d+2); 
            % And set snooker alfa of current chain positions to zero
            log_alfa_sn_Xr(id_c) = 0;
            % Compute log_wXr, Eq. 10 of https://arxiv.org/pdf/1801.09065 
            log_wXr = reshape(sum(Xr(1:N,DREAMPar.d+1:DREAMPar.d+2),2) - ...
                log_alfa_sn_Xr,DREAMPar.mt,DREAMPar.N);
            % Accept/not cand. points (JAV: extended for informtve priors)
            [accept,id,ii] = mtMetropolis_rule(log_wXp,log_wXr, ...
                DREAMPar,id_pop);
        otherwise
            % Calculate Metropolis ratio
            [accept,id] = Metropolis_rule(DREAMPar,Meas_info, ...
                log_alfa_sn,Xp(1:N,DREAMPar.d+1:DREAMPar.d+2), ...
                X(1:N,DREAMPar.d+1:DREAMPar.d+2),options); ii = id;
    end
    % Accept candidate points
    X(ii,1:DREAMPar.d+2) = Xp(id,1:DREAMPar.d+2); 
    % Do same with simulated outputs and likelihood vector (if time series)
    FX(:,ii) = FXp(:,id); LLX(:,ii) = LLp(:,id);
    % And summary metrics
    if ~isempty(S), S(:,ii) = SXp(:,id); end

    % Update residuals and measurement errors for Kalman proposal
    if strcmp(method,'dream_kzs')
        std_mX(:,ii) = std_mXp(:,id); EX(:,ii) = EXp(:,id);
    end
    % Store current states in chains?
    if mod(t,DREAMPar.thinning) == 0
        % Store current state in chain
        iloc = iloc + 1; chain(iloc,1:DREAMPar.d+2,1:DREAMPar.N) = ...
            reshape(X',1,DREAMPar.d+2,DREAMPar.N);
        % Store model simulations (if appropriate)
        DREAM_Suite_store_results(options,[FX;S],Meas_info,'a+');
    end
    % Update crossover values?
    if strcmp(DREAMPar.adapt_pCR,'yes')
        % Standardized Euclidean distance new and old population
        sdX_Xp(:,g) = sum(bsxfun(@rdivide,X(1:DREAMPar.N,1:DREAMPar.d) ...
            - X_old,max(std(X_old),eps)).^2,2);
    end
    % Append X to external archive, Z?
    if any(strcmpi(method,{'dream_zs','dream_dzs','dream_kzs', ...
            'mtdream_zs'}))
        if mod(t,DREAMPar.k) == 0
            % Append X to Z
            fid_Z = fopen('Z.bin','a+','n');
            fwrite(fid_Z,X','double'); fclose(fid_Z);
            % Append FX to FXZ
            fid_FXZ = fopen('FXZ.bin','a+','n');
            fwrite(fid_FXZ,FX,'double'); fclose(fid_FXZ);
            % Update DREAMPar.m
            DREAMPar.m = DREAMPar.m + DREAMPar.N;
        end
        % Adapt Kalman jump probability?
        if strcmp(method,'dream_kzs')
            if ( t >= round(DREAMPar.a_1 * DREAMPar.T) ) && ...
                    ( t < round(DREAMPar.a_2 * DREAMPar.T) )
                % No more use of Kalman jump distribution
                DREAMPar.pkalman = DREAMPar.oldpkalman;
                DREAMPar.pparallel = 1 - DREAMPar.psnooker ...
                    - DREAMPar.pkalman;
            else
                DREAMPar.pkalman = 0; DREAMPar.pparallel = 1 - ...
                    DREAMPar.psnooker - DREAMPar.pkalman;
            end
        end
    end
    % Update gen and # candidate points accepted
    g = g + 1; tot_accept = tot_accept + sum(accept);
    % Update log likelihood of chains [= for outlier chain detection only]
    loglik(t,1:DREAMPar.N+1) = [t*DREAMPar.N X(:,DREAMPar.d+2)'];
    % If t is equal to MCMC.steps then compute convergence diagnostics
    if mod(t,DREAMPar.steps) == 0
        % Store acceptance Rate
        output.AR(it,1:2) = [ t * DREAMPar.N 100 * sum(tot_accept) / ...
            (DREAMPar.N * DREAMPar.steps)];
        % Update pCR values?
        if (t < DREAMPar.T/5)
            if strcmp(DREAMPar.adapt_pCR,'yes')
                [TdCR,cCR,pCR] = Calc_pCR(DREAMPar,sdX_Xp,...
                    TdCR,cCR,CR);
            end
        end
        if any(strcmp(method,{'dream','dream_d'}))
            % Outlier chains? If so, reset them to current best X of chains
            [X,loglik(floor(t/2):t,2:DREAMPar.N+1),outlier] = ...
                Remove_outlier(method,DREAMPar,X,t,...
                loglik(floor(t/2):t,2:DREAMPar.N+1),options);
            % Append ID of outlier chains
            output.outlier = [ output.outlier ; outlier ];
        end
        % Store individual crossover values
        output.CR(it,1:DREAMPar.nCR+1) = [ t * DREAMPar.N pCR ];
        % Sample randomly crossover values of 1:DREAMPar.nCR
        CR = reshape(randsample(1:DREAMPar.nCR,DREAMPar.N * ...
            DREAMPar.steps,'true',pCR),DREAMPar.N,DREAMPar.steps);
        % Now it is time for convergence diagnostics
        start_id = max(1,floor(options.burnin/100 * iloc)); end_id = iloc;
        % hatR and hatRd diagnostics using chain burn-in (DEF = 50%)
        [hatR,hatRd] = Gelman(chain(start_id:end_id,1:DREAMPar.d,...
            1:DREAMPar.N),t,method);
        % Append univariate \hat{R} convergence diagnostic
        output.R_stat(it,1:DREAMPar.d+1) = [ t*DREAMPar.N hatR ];
        % Append multivariate \hat{R}d convergence diagnostic
        output.MR_stat(it,1:2) = [ t*DREAMPar.N hatRd ];
        % Print progress
        if t > 1
            fprintf(1, repmat('\b',1,count)); % delete line before
            str = strcat('count = fprintf(''',upper(method),{' '},...
                'CALCULATING:',{' '},...
                ['%3.2f%% done and \\hat{R}^{d} convergence ' ...
                'diagnostic: '],{' '},'%3.2f''', ...
                ',100*(t/DREAMPar.T)',...
                ',output.MR_stat(it,2));'); eval(char(str));
        end
        % Update iter, set generation back to 1 and totaccept to zero
        it = it + 1; g = 1; tot_accept = 0;
        % Save memory or not?
        if strcmp(options.save,'yes'), save DREAM_Suite.mat; end
    end
end
% <><><><><><><><><><><><> END DYNAMIC PART <><><><><><><><><><><><><><><>

% <><><><><><><><><><><><>< POSTPROCESSING ><><><><><><><><><><><><><><><>
% Determine total run time
output.RunTime = cputime - start_time; fprintf('\n');
% Run postprocessor (tables and figures)
[chain,output,FX,Z] = DREAM_Suite_postproc(method,DREAMPar,f_handle,...
    Par_info,Meas_info,Lik_info,options,chain,output,it,iloc,MAP_info);
% Close the cluster - if used
DREAM_Suite_end(DREAMPar,options);
% Now define return arguments
varargout(1:5) = {chain,output,FX,Z,loglik};
% <><><><><><><><><><><>< END POSTPROCESSING <><><><><><><><><><><><><><><

end