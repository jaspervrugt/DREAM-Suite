function [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...
    GAMEoptions,Par_info,Meas_info,options,plugin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% GGGG   AA   M   M  EEEE  SSSS   AA   M   M  PPPP  L     I  N    N  GGGG %
% G  G  A  A  M   M  E     S     A  A  M   M  P  P  L     I  NN   N  G  G %
% G  G  A  A  MM MM  EEE   S     A  A  MM MM  P  P  L     I  N N  N  G  G %
% GGGG  AAAA  M M M  EEE - SSSS  A  A  M M M  PPPP  L     I  N  N N  GGGG %
%    G  A  A  M   M  E        S  AAAA  M   M  P     L     I  N   NN     G %
% GGGG  A  A  M   M  EEEE  SSSS  A  A  M   M  P     LLLL  I  N    N  GGGG %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% This function computes the marginal likelihood from a collection of     %
% samples, X, of the target distribution derived from (e)DREAM Package.   %
% Type "help GAME_fit_mixture" for more information about the procedure   %
%                                                                         %
% MAIN IDEA:                                                              %
%  Calculate the marginal likelihood as weighted mean of the ratio of the %
%  samples' target density, p(x_i), and the importance density, q(x_i).   %
%  This works because we know that the importance distribution, a mixture %
%  of normal distributions, integrates to 1. As a result, we yield that   %
%  Z = 1/N * sum_{i=1}^{N} (p(x_i)/q(x_i)) or as we solve herein in       %
%  logarithmic form (numerically more stable)                             %
%                                                                         %
% SYNOPSIS:                                                               %
%  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name);            %
%  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         %
%      GAMEoptions);                                                      %
%  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         %
%      GAMEoptions,Par_info);                                             %
%  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         %
%      GAMEoptions,Par_info,Meas_info);                                   %
%  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         %
%      GAMEoptions,Par_info,Meas_info,options);                           %
%  [Z,logZ,gmix] = GAME_sampling(X,method,DREAMPar,Func_name, ...         %
%      GAMEoptions,Par_info,Meas_info,options,plugin);                    %
% WHERE                                                                   %
%  X           [input] Rx(d+2) matrix posterior samples DREAM             %
%  method      [input] string (name) marginal likelihood estimator        %
%   = 'ris'        reciprocal importance sampling                         %
%   = 'is'         importance sampling                                    %
%   = 'ob'         optimal bridge sampling                                %
%   = 'gb'         geometric bridge sampling                              %
%  DREAMPar    [input] Structure with algorithmic variables               %
%   .d             Dimensionality (# variables) target distribution       %
%   .N             # of Markov chains                                     %
%   .T             # of generations (= # samples of each Markov chain)    %
%   .lik           Choice of likelihood function                          %
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
%                   → DREAM_ZS and DREAM_DZS                              %
%   .m0            Initial size of external archive, Z   DEF: 10*d        %
%                   → DREAM_ZS and DREAM_DZS                              %
%   .k             Growth rate of external archive       DEF: 10          %
%                   → DREAM_ZS and DREAM_DZS                              %
%   .M             # samples archive Z for Kalman jump   DEF: 20          %
%                   → DREAM_KZS                                           %
%   .mt            Number of multi-try proposals         DEF: 5           %
%                   → MTDREAM_ZS                                          %
%  Func_name   [input] Function (string) returns (log)lik or sim values   %
%   → Func_name must return a likelihood if DREAMPar.lik = 1              %
%   → Func_name must return a log-likelihood if DREAMPar.lik = 2          %
%   → Func_name returns vector of n simulated values: DREAMPar.lik > 2    %
%  GAMEoptions [input] (optional) GAME structure with additional options  %
%   = struct('metric','BIC','K',5,'N',1e4,'M',1,'steps',10);              %
%   .metric        metric for optimal mixture selection                   %
%    = 'bic'       Bayesian information criterion        DEFault          %
%    = 'var'       Variance reduction                                     %
%   .K             maximum # components mixture dist.    DEF: 5           %
%   .N             # importance samples to be used       DEF: 10,000      %
%   .M             # times we repeat each estimator      DEF: 1           %
%   .steps         # steps with gb or ob                 DEF: 10          %
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
%  plugin      [input] 2nd input argument Func_name. Class set by user    %
%  Z           [outpt] marginal likelihood                                %
%                  = integral of posterior pdf                            %
%                  = integral of p(x|Y)                                   %
%                  = integral of p(x)L(x|Y)                               %
%  logZ        [outpt] logarithmic value of marginal likelihood           %
%  gmix        [outpt] Structure trained normal mixture importance dist.  %
%  .n              # samples (= R) of X                                   %
%  .d              # dimensions (= parameters) of X                       %
%  .k              # components of mixture distribution                   %
%  .p              # parameters of mixture distribution                   %
%  .w              maximum likelihood weights of mixture components       %
%  .mu             jxd matrix of mean values each component               %
%  .Sigma          dxdxK array of covariance matrices each component      %
%  .I              integral of pdf each component before applying weight  %
%  .loglik         log-likelihood of normal mixture                       %
%  .AIC            Akaike information criterion of normal mixture         %
%  .BIC            Bayesian information criterion of normal mixture       %
%                                                                         %
% REFERENCE:                                                              %
%   Volpi, E., G. Schoups, G. Firmani, and J.A. Vrugt (2017), Sworn       %
%       testimony of the model evidence: Gaussian Mixture Importance      %
%       (GAME) sampling, Water Resources Research, 53, pp. 6133-6158,     %
%       https://doi.org/10.1002/2016WR020167                              %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% Written by Jasper A. Vrugt & Elena Volpi & Gerrit Schoups               %
%                                                                         %
% Version 1:    June 2012       Initial setup and definition              %
% Version 1.1:  Jan. 2015       Major overhaul of code                    %
% Version 1.2:  Oct. 2017       Final update                              %
% Version 1.3:  Aug. 2021       Return logZ as well                       %
% Version 2.0:  Aug. 2024       Major overhaul of code                    %
%                                                                         %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
% NOTE 1:                                                                 %
%  Geometric bridge sampling the code uses:                               %
%      t = linspace(1,0,steps);                                           %
%  Optimal bridge sampling the code uses:                                 %
%      nn1 = size(logQ_Xp(R1),1) * linspace(1,0.001,steps);               %
%  If t = 0; --> RIS whereas if t = 1 --> IS                              %
%      t in [0,1] bridge between both end members                         %
%                                                                         %
% NOTE 2:                                                                 %
%  Use a proper likelihood, that is, normalization constant included!!    %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%                                                                         %
%  Step 1: Run example 2 of DREAM package                                 %
%  Step 2: Then after completion do: P = genparset(chain);                %
%  Step 3: Get posterior samples: X = P(end-25000:end,1:DREAMPar.d+2);    %
%  Step 4: Calculate marginal likelihood via GAME_sampling                %
%             [Z,logZ] = GAME_sampling(X,method,DREAMPar,Func_name, ...   %
%                 GAMEoptions,Par_info);                                  %
%             where method = any option from {'ris','is','ob','gb'}       %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

MAP_info = struct;                          % Needs to be treated later
% Check input arguments
if nargin < 9
    fprintf(['GAME_sampling WARNING: plug is empty: no additional ' ...
        'input to target function\n']); plugin = [];
end
if nargin < 8 || isempty(options)
    fprintf('GAME_sampling WARNING: Structure options is empty\n');
    options = struct('restart','no','parallel','no');
end
if nargin < 7 || isempty(Meas_info)
    fprintf('GAME_sampling WARNING: Structure Meas_info is empty\n');
    Meas_info = struct;
end
if nargin < 6 || isempty(Par_info) || ~isfield(Par_info,'min')
    fprintf(['GAME_sampling WARNING: Par_info.min not specified: ' ...
        'Minimum values of parameters set to -inf\n']);
    Par_info.min = -inf(1,DREAMPar.d);
end
if nargin < 6 || isempty(Par_info) || ~isfield(Par_info,'max')
    fprintf(['GAME_sampling WARNING: Par_info.max not specified: ' ...
        'Maximum values of parameters set to inf\n']);
    Par_info.max = +inf(1,DREAMPar.d);
end
if nargin < 5 || isempty(GAMEoptions)
    fprintf(['GAME_sampling WARNING: Structure GAMEoptions is empty: ' ...
        'Default settings assumed\n']);
    GAMEoptions = struct;
end

[method,DREAMPar,Par_info,options, ...      % Check input variables
    GAMEoptions,R,d2] = GAME_check( ...
    method,X,Func_name,DREAMPar, ...
    GAMEoptions,Par_info,Meas_info, ...
    options);
[DREAMPar,metric,K,N,M,steps,Par_info, ...  % Init. main variables of DREAM
    Meas_info,Lik_info,options] = ...
    GAME_setup(DREAMPar,Func_name, ...
    GAMEoptions,Par_info,Meas_info, ...
    options);
X_un = X_unnormalize(X(:,1:DREAMPar.d), ... % Unnormalized post. par values
    Par_info);
logP = sum(X(1:R,d2-1:d2),2);               % Logarithmic posterior density
X = X_un(1:R,1:DREAMPar.d);                 % Extract posterior samples
gmix = GAME_fit_mixture(X,logP, ...         % Train Gaussian mixture up to 
    metric,K,Par_info.min,Par_info.max);    % J components
logQ = GAME_mixture_density(gmix,X);        % Log density importance dist.
logZ = nan(M,1);                            % Initialize logZ
 
for j = 1 : M   % DYNAMIC PART
    R1 = randsample(1:round(R/2), ...       % RIS/GB: draw DREAM samples
        min(round(R/2),N));                 % without those ID estimation!)
    if strcmp(method,'ris') % → Reciprocal Importance Sampling
        logZ(j,1) = ...                     % qhalf = q0 (RIS)
            GAME_logevidence(0,[], ...
            [],logQ(R1),logP(R1));
        fprintf(['GAME_sampling ' ...       % print warning
            'WARNING: RIS samples ' ...
            'taken from first half ' ...
            'of x_p: make sure you ' ...
            'use second half for ' ...
            'mixture fitting\n']);
    else                    % → Importance/Geometric/Optimal Bridge Sampling
        DREAMPar.N = N;                     % # chains is equal to N
        [DREAMPar,f_handle] = ...           % Function handle and do setup
            eDREAM_package_calc_setup( ...
            DREAMPar,Func_name,options, ...
            plugin);
        X_gmix = mvgmmrnd(gmix,N);          % Draw unnrmlzd samples mixture
        % <> do boundary handling <> ?      % Satisfy boundary constraints
   
        % Now check whether we used DREAM_D or DREAM_DZS
        if isfield(Par_info,'steps')
            fprintf(['GAME_sampling'...     % Warning discrte imprtnce sam.
                ' WARNING: Discrete ' ...
                'importance sampling ' ...
                'as field steps of ' ...
                'structure Par_info ' ...
                'is specified\n']);
            Par_info.step_size = ...        % Now lets define step size
                (Par_info.max - ...
                Par_info.min) ./ ...
                Par_info.steps;
            X_gmix = ...                    % Trans. mixture discrete space
                Discrete_space( ...
                X_gmix,Par_info);
        end
        logQ_gmix = ...                     % log(density) imprtnce dist.
            GAME_mixture_density( ...   
            gmix,X_gmix);
        FX_gmix = Evaluate_target( ...      % Evaluate model and return fx
            X_gmix,DREAMPar, ...
            Meas_info,options, ...
            f_handle);
        logPR_gmix = Eval_prior(...         % Compute log(prior) mix. pts
            X_gmix,[],Par_info, ...
            Meas_info,options);
        logL_gmix = ...                     % Compute loglik mixture points
            Calc_likelihood(X_gmix, ...
            FX_gmix,DREAMPar,Par_info, ...
            Meas_info,Lik_info, ...
            options,MAP_info);
        logP_X_gmix = logPR_gmix ...        % Log post. densty mixture pts.
            + logL_gmix;

        switch method
            case {'is'}  % → Importance Sampling
                logZ(j,1) = ...
                    GAME_logevidence( ...
                    1,logQ_gmix, ...
                    logP_X_gmix, ...
                    logQ(R1),logP(R1));
            case {'gb'}  % → Geometric Bridge Sampling
                t = linspace(1,0,steps);
                for k = 1 : steps
                    logZ(j,k) = ...
                        GAME_logevidence( ...
                        t(k),logQ_gmix, ...
                        logP_X_gmix, ...
                        logQ(R1),logP(R1));
                end
            case {'ob'}  % → Optimal Bridge Sampling
                logZ(j,1) = ...
                    GAME_logevidence( ...
                    100,logQ_gmix, ...
                    logP_X_gmix, ...
                    logQ(R1),logP(R1));
                nn1 = size(logQ(R1),1) * ... % Optimal Bridge Sampling
                    linspace(1,0.001,steps);
                for k = 2 : steps
                    logZ(j,k) = ...
                        GAME_logevidence( ...
                        100,logQ_gmix, ...
                        logP_X_gmix,logQ(R1), ...
                        logP(R1),nn1(k));
                end
        end
    end             % END: GAME SAMPLING
end             % END: DYNAMIC PART

Z = exp(logZ);                              % Transform to model evidence

end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% SECONDARY FUNCTIONS
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
function X = mvgmmrnd(gmix,N)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%MVGMMRND Multivariate Gaussian mixture model random. This function draws %
% N samples from the normal mixture distribution defined in the structure %
% gmix. The normal mixture PDF at x [= 1xd vector] is written as          %
%          PDF(x) = Σ_{j=1}^{k} w_{j} f_{N_{d}}(x,µ_{j},Σ_{j})            %
% where f_{N_{d}}(x,µ_{j},Σ_{j}) is the d-variate normal pdf with mean µ  %
% [= 1xd vector] and dxd covariance matrix Σ. In MATLAB, we yield that    %
% f_{N_{d}}(x,µ_{j},Σ_{j}) = mvnpdf(x,µ_{j},Σ_{j})                        %
%                                                                         %
% SYNOPSIS:                                                               %
%  X = mvgmmrnd(gmix,N);                                                  %
% WHERE                                                                   %
%  gmix        [outpt] Structure trained normal mixture importance dist.  %
%  .d              # dimensions (= parameters) of X                       %
%  .k              # components of mixture distribution                   %
%  .w              maximum likelihood weights of mixture components       %
%  .mu             kxd matrix of mean values each component               %
%  .Sigma          dxdxk array of covariance matrices each component      %
%  N           [input] # samples drawn from normal mixture ditribution    %
%  X           [outpt] Nxd matrix of samples drawn from normal mixture    %
%                                                                         %
% © Written by Jasper A. Vrugt, Jan. 2015                                 %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

nw = histcounts(rand(N,1), ...              % # draws each component
    [0; cumsum(gmix.w(:))./sum(gmix.w)]);
cum_nw = [0 cumsum(nw)];                    % Cumulative count
X = nan(N,gmix.d);                          % Initialize Nxd matrix X
for j = 1:gmix.k                            % Draw jth comp. based on w(j) 
    id = cum_nw(j) + 1 : cum_nw(j+1);       % Row id of samples jth comp.
    X(id,1:gmix.d) = randn(nw(j), ...       % Draw nw(j) samples jth comp.
        gmix.d) * chol(gmix.Sigma(:,:,j)) ...
        + repmat(gmix.mu(j,:),nw(j),1);
end

end