function [loglik,std_e,eps_n,f_eps_n,Y_r] = UL(iflag,nuisvar,method,...
    y_sim,y_meas,sigma,N)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Universal likelihood (UL) function: correlated, heteroscedastic and non-Gaussian   %%
%% errors                                                                             %%
%%  1. Partial stand. residuals expected to follow skew generalized Student's t (SGT) %%
%%     distribution                                                                   %%
%%  2. Serial correlation described using an autoregressive model - up to order 2     %%
%% UL function performs standardization before treatment serial correlation           %%
%%                                                                                    %%
%%  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r] = GL_plus(iflag,nuisvar,method,...     %%
%%                 y_sim,Y_meas,sigma,N)                                              %%
%%  where                                                                             %%
%%   iflag      [input] Estimation ('est') or simulation ('sim')                      %%
%%   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         %%
%%    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       %%
%%    s1:  nuisvar(2) Slope of linear heteroscedastic model                           %%
%%    lbd: nuisvar(3) Skewness (-1,1): (-1,0) neg. skew, 0 symmetric, (0,1) pos. skew %%
%%    p:   nuisvar(4) Kurtosis (> 0) The larger p is, the more uniform distribution   %%
%%    q:   nuisvar(5) Kurtosis (> 0) The larger q is, the peakier the distribution    %%
%%    fi1: nuisvar(5) First-order AR coefficient (0,1: check book)                    %%
%%    fi2: nuisvar(7) Second-order AR coefficient (0,1: check book)                   %%
%%   method     [input] Treatment of s0 and/or s1: [1-8]                              %%
%%   y_sim      [input] n x 1 vector with simulated values                            %%
%%   Y_meas     [input] n x 1 vector with observed values                             %%
%%   sigma      [input] Optional: Measurement sigma defined by user                   %%
%%   N          [input] Optional: Number of replicates - resampling                   %%
%%   loglik     [outpt] n x 1 vector of log-likelihood values                         %%
%%   std_e      [outpt] n x 1 vector with standard deviation of raw residuals         %%
%%   eps_n      [outpt] n x 1 vector with standardized decorrelated residuals         %%
%%   f_eps_n    [outpt] n x 1 vector with density stand. decorrelated residuals       %%
%%   Y_r        [outpt] n x N matrix of replicate simulations (for Bayes_pdf)         %%
%%                                                                                    %%
%%  Reference:                                                                        %%
%%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of %%
%%       distribution-adaptive likelihood functions: Generalized and universal        %%
%%       likelihood functions, scoring rules and multi-criteria ranking, Journal of   %%
%%       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             %%
%%       https://www.sciencedirect.com/science/article/pii/S002216942201112X          %%
%%                                                                                    %%
%%  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  %%
%%         2. Measurement error computed from simulated data                          %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Dec. 2016                                            %%
%% DREAM Package                                                                      %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% If number of replicates, N, is not defined
if nargin < 7, N = 1; end

% Initialize return arguments
loglik = -inf*ones(size(y_meas)); [std_e,eps_n,f_eps_n,Y_r] = deal([]);

% Unpack nuisance variables
nuisvar = num2cell(nuisvar); [s0,s1,lbd,p,q,fi1,fi2] = deal(nuisvar{:});

% Determine mean and variance of SGT distribution
if p*q < 2  % Do not admit a zero-mean and unit variance SGT distribution
    return
else        % we expect standardized partial residuals to be distributed
    % according to the standardized skewed generalized t
    % distribution, SGT(.|0,1,lambda,p,q) if
    kappa_lpq = beta(1/p,q/p) / sqrt((1+3*lbd^2).*...
        beta(1/p,q/p) * beta(3./p,(q-2)/p) - ...
        4*lbd^2 * beta(2/p,(q-1)/p).^2);
    mu_lpq = 2*kappa_lpq * lbd * 1 * beta(2/p,(q-1)/p) / ...
        beta(1/p,q/p);
end

n = size(y_sim,1);      % # elements of model simulation
fi_p = [1 -fi1 -fi2];   % Coefficients of AR filter
e = y_meas - y_sim;     % nx1 vector of raw residuals

% Studentize/standardize the raw residuals: measurement error std.
switch method
    case 0 % CONSTANT/NONCONSTANT: specifed by user: s0 and s1 irrelevant
        std_e = sigma;
        % --> s0 removed, s1 removed
    case 1 % CONSTANT: s0 not selected, s1 not selected
        std_e = s0 * ones(n,1);
        % --> s0 fixed, s1 removed
    case 2 % CONSTANT: s0 not selected, s1 selected
        std_e = std(e) * ones(n,1);
        % --> s0 phantom variable, s1 removed
    case 3 % CONSTANT: s0 selected, s1 not selected
        std_e = s0 * ones(n,1);
        % --> s0 estimated, s1 removed
    case 4 % CONSTANT: s0 selected, s1 is selected
        std_e = s0 * ones(n,1);
        % --> s0 estimated, s1 removed
    case 5 % NONCSTNT: s0 not selected, s1 not selected
        [std_e,~,exitflag] = s1_phantom(s0,e,y_sim);
        if exitflag ~= 1, return; end
        % --> s0 fixed, s1 phantom variable
    case 6 % NONCSTNT: s0 not selected, s1 selected
        std_e = s0 + s1 * y_sim;
        % --> s0 fixed, s1 estimated
    case 7 % NONCSTNT: s0 selected, s1 not selected
        [std_e,~,exitflag] = s1_phantom(s0,e,y_sim);
        if exitflag ~= 1, return; end
        % --> s0 estimated, s1 phantom variable
    case 8 % NONCSTNT: s0 selected, s1 selected
        std_e = s0 + s1 * y_sim;
        % --> s0 estimated, s1 estimated
end

% Compute theoretical std. of partial residuals
std_eps = sqrt( ( 1 + fi1^2 - fi2^2 - 2*fi1^2/(1-fi2) ) * 1 );

% Multiplication with one should be replaced by multiplication with the
% variance of the studentized raw residuals. Only cases 2,5,7 guarantee a
% unit variance of the studentized raw residuals, e_n = e/std_e,
% hence * 1 is fine only for these three cases. Practical experience
% suggests that free estimation of s0/s1 leads to equivalent results

%% Now check flag
switch iflag
    case 'est'  % Likelihood estimation
        Y_r = nan(n,1);
        if ~isreal(std_eps), return; end    % abs(fi1) + abs(fi2) <= 1!!
        e_n = e./std_e;                     % Studentized raw residuals; var(e_n) = 1
        eps = filter(fi_p,1,e_n);           % Partial residuals
        eps_n = eps./std_eps;               % Standardized partial residuals
        eps_sgt = eps_n + mu_lpq;           % Shift for SGT distribution
        loglik = - log(std_eps) - ...       % Individual log-likelihoods
            log(std_e) + log(p) - log(2) ...
            - log(kappa_lpq) - betaln(1/p,q/p) ...
            - ((q + 1)/p) * log (1 + ...
            abs( eps_sgt ./ (kappa_lpq * ...
            (1 + lbd*sign(eps_sgt))) ).^p );
        if nargout > 3
            f_eps_n = f_SGT(eps_n,lbd,p,q); % SGT density of eps_n
        end

    case 'sim'  % Simulation: Replicates, Y_r, for prediction uncertainty
        eps_n = SGTrnd(0,1,lbd,p,q,[n N]);  % nxN matrix samples from SGT(0,1,lbd,p,q)
        dl_eps = zeros(n,N);                % Correction to mean of eps
                                            % Scharnagl 2015: dl_eps = fi1 * mean_of_e_n
        eps = std_eps * eps_n + dl_eps;     % nxN matrix residuals
        e_n = filter(1,fi_p,eps);           % nxN matrix std. raw residuals, var(e_n) = ones(1,N)
        E_r = bsxfun(@times,std_e,e_n);     % nxN matrix raw residuals
        Y_r = bsxfun(@plus,y_sim,E_r);      % N replicates of simulation

end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% NOTES: 1. This is an approximate log-likelihood: The preceding two observations    %%
%%           y_(-1) and y_0 are assumed zero in AR model                              %%
%%        2. Crucial aspect is how we compute std_res; from simulated values (as in   %%
%%           Schoups and Vrugt, 2010) or from measured data. Latter approach may      %%
%%           converge faster.                                                         %%
%%        3. Should the studentized raw residuals and partial residuals have the      %%
%%           same mean value? See Scharnagl et al., 2015 (HESS-D paper). Tested,      %%
%%           but maybe worthwhile to revisit. The suggested correction poses a        %%
%%           problem during resampling, what correction would one apply to the        %%
%%           resampled time series?                                                   %%
%% Re: 2     Do we express heteroscedasticity using measured or simulated data?       %%
%%        a. If we use simulated data then prediction intervals may look strange at   %%
%%           times as heteroscedasticity follows simulation and not measured data     %%
%%        b. If we use measured data, and model is "bad", say we significantly        %%
%%           underpredict the current observation, then meas_error will be            %%
%%           unexpectedly large. Vice, versa if model output is much larger than      %%
%%           observation, the meas. error may be unexpectedly small.                  %%
%%           If we use measured data than phi --> 1 as then model/meas error is       %%
%%           relatively small (you can see this for French Broad --> s1 = 0.5 is      %%
%%           unrealistically large)                                                   %%
%%                                                                                    %%
%% The above issues are not that important if model can track data relatively well,   %%
%% otherwise, inference may produce undesirable credible/prediction intervals.        %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

