function [loglik,std_e,eps_n,f_eps_n,Y_r] = SL(iflag,nuisvar,method,...
    y_sim,y_meas,sigma,N)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Student Likelihood (SL) function: correlated, non-constant and non-Gaussian errors %%
%%  1. Partial stand. residuals expected to follow skew Student t (SST) distribution  %%
%%  2. Serial correlation treated using an AR(1) or AR(2) model                       %%
%%                                                                                    %%
%% This SST function performs standardization before treatment correlation            %%
%%                                                                                    %%
%%  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r] = SL(iflag,nuisvar,method,y_sim,...    %%
%%                 y_meas,sigma,N)                                                    %%
%%  where                                                                             %%
%%   iflag      [input] Estimation ('est') or simulation ('sim')                      %%
%%   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         %%
%%    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       %%
%%    s1:  nuisvar(2) Slope of linear heteroscedastic model                           %%
%%    nu:  nuisvar(3) Degrees of freedom (q > 0)                                      %%
%%    xi:  nuisvar(4) Skewness parameter (xi > 0)                                     %%
%%    fi1: nuisvar(5) First-order AR coefficient (0,1: check book)                    %%
%%    fi2: nuisvar(6) Second-order AR coefficient (0,1: check book)                   %%
%%   method     [input] Treatment of s0 and/or s1: [1-8]                              %%
%%   y_sim      [input] n x 1 vector with simulated values                            %%
%%   y_meas     [input] n x 1 vector with observed values                             %%
%%   sigma      [input] Optional: Measurement sigma defined by user                   %%
%%   N          [input] Optional: Number of replicates - resampling                   %%
%%   loglik     [outpt] n x 1 vector of log-likelihood values                         %%
%%   std_e      [outpt] n x 1 vector with standard deviation of raw residuals         %%
%%   eps_n      [outpt] n x 1 vector with standardized decorrelated residuals         %%
%%   f_eps_n    [outpt] n x 1 vector with density stand. decorrelated residuals       %%
%%   Y_r        [outpt] n x N matrix of replicate simulations (for Bayes_pdf)         %%
%%                                                                                    %%
%%  Reference:                                                                        %%
%%   Scharnagl, B., S.C. Iden, W. Durner, H. Vereecken, and M. Herbst (2015), Inverse %%
%%       modelling of in situ soil water dynamics: accounting for heteroscedastic,    %%
%%       autocorrelated, and non-Gaussian distributed residuals, Hydrology and Earth  %%
%%       System Sciences Discussions, vol. 12, pp. 2155-2199, 2015.                   %%
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
if nargin < 7, N = 1; end; calc_method = 1;

% Initialize return arguments
loglik = -inf*ones(size(y_meas)); [std_e,eps_n,f_eps_n,Y_r] = deal([]);

% Unpack nuisance variables
nuisvar = num2cell(nuisvar); [s0,s1,nu,xi,fi1,fi2] = deal(nuisvar{:});

% Compute first and second moment of SST distribution
switch calc_method
    case 0 % OK for small nu; if nu > 100 or so, M1 and M2 not defined
        M = @(j,nu) gamma((j+1)/2)*gamma((nu-j)/2)*(nu-2)^(j/2) / ...
            (sqrt(pi)*gamma(nu/2));
        M1 = M(1,nu); M2 = M(2,nu);
    case 1 % OK for all nu --> log formulation to avoid numerical overflow
        logM = @(j,nu) gammaln((j+1)/2) + gammaln((nu-j)/2) + ...
            (j/2)*log(nu-2) - 1/2*log(pi) - gammaln(nu/2);
        M1 = exp(logM(1,nu)); M2 = exp(logM(2,nu));
end
% Now compute mu_xi and sig_xi
mu_xi = M1*(xi - 1/xi); % Matches Scharnagl paper
sig_xi = sqrt((M2-M1^2)*(xi^2 + xi^-2) + 2*M1^2 - M2);

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

% Now check flag
switch iflag
    case 'est'  % Likelihood estimation
        Y_r = nan(n,1);                         % Replicate simulation
        if ~isreal(std_eps), return; end        % abs(fi1) + abs(fi2) <= 1!!
        e_n = e./std_e;                         % Studentized raw residuals; var(e_n) = 1
        eps = filter(fi_p,1,e_n);               % Partial residuals
        eps_n = eps./std_eps;                   % Standardized partial residuals
        a_sst = (mu_xi + sig_xi.*eps_n) ./ ...  % Adjust to SST distribution
            (xi.^sign(mu_xi + sig_xi.*eps_n));
        % OK for all degrees of freedom, nu > 2
        % gammaln(x) avoids gamma(x) to inf for large x
        loglik = - log(std_e) - log(std_eps)... % Individual log-likelihoods
            + log(2) + log(sig_xi) + gammaln((nu+1)/2) ...
            - log(xi+1/xi) - gammaln(nu/2) - 1/2*log(pi) ...
            - 1/2*log(nu-2) - (nu+1)/2*log( 1 ...
            + 1/(nu-2) * a_sst.^2 );
        % Compute SST density
        if nargout > 3
            f_eps_n = f_SST(eps_n,nu,xi,calc_method);
        end

    case 'sim'  % Simulation: Replicates, Y_r, for prediction uncertainty
        eps_n = SSTrnd(nu,xi,[n N]);    % nxN matrix samples from SST(0,1,nu,xi)
        dl_eps = zeros(n,N);            % Correction to mean of eps
                                        % Scharnagl 2015: dl_eps = fi1 * mean_of_e_n
        eps = std_eps * eps_n + dl_eps; % nxN matrix residuals
        e_n = filter(1,fi_p,eps);       % nxN matrix std. raw residuals, var(e_n) = ones(1,N)
        E_r = bsxfun(@times,std_e,e_n); % nxN matrix raw residuals
        Y_r = bsxfun(@plus,y_sim,E_r);  % N replicates of simulation

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES
% 1: https://www.vosesoftware.com/riskwiki/Laplacedistribution.php
%    f(x|a,sigma) = 1/(sqrt(2)*sigma ) * exp(-sqrt(2)|x-a|/sigma)
%                 = prod_t=1^n 1/(sqrt(2)*sigma ) * exp(-sqrt(2)|x-a|/sigma)
%                 = (1/(sqrt(2)*sigma ))^n * prod_t=1^n exp(-sqrt(2)|x-a|/sigma)
% log(f(x|a,sigma)) = -n*log(sqrt(2)*sigma) + log(exp( sum_t=1^n -sqrt(2)|x-a|/sigma ))
%                   = -n*log(sqrt(2)*sigma) - sqrt(2)/sigma * sum_t=1^n |x-a|
%                   = -n/2*log(2) - n*log(sigma) - sqrt(2)*sum_t=1^n |x-a|/sigma
% 2: 1/std(eps) * lapdf(eps_n,0,1) = lapdf(eps,0,std(eps))
% 3: First observation has a variance: sigma2_eps/(1 - phi1^2);
% 4: Later observations have variance: sigma2_eps --> from normal likelihood
% 5: Combine: - n/2*log(2) - log(sqrt(std_e(1)^2/(1 - phi1^2))) ...
%             - (n-1)*log(sigma_eps) - sqrt(2) * abs(e_n(1)) / ...
%               sqrt(sigma2_eps/(1 - phi1^2)) ...
%             - sqrt(2)*sum (abs(eps(2:n))/sigma_eps) - sum(log(std_e));
% 6: GL+/UL/SL large sample approximation
%    logL_Xp(ii,1) = - n/2 * log(2) - n*log(sigma_eps) ...
%                    - sqrt(2) * sum(abs(eps_n(1:n))) ...
%                    - sum(log(std_e),'omitnan');
% This matches if we treat first observation as zeros
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
