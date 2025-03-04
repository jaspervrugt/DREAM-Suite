function [loglik,std_e,eps_n,f_eps_n,Y_r] = NL(iflag,nuisvar,method,...
    y_sim,y_meas,sigma,N)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Normal Likelihood (NL) function: correlated and non-constant errors                %%
%%  1. Partial stand. residuals expected to follow normal distribution                %%
%%  2. Serial correlation treated using an AR(1) or AR(2) model                       %%
%%                                                                                    %%
%% NL function performs standardization before treatment correlation                  %%
%%                                                                                    %%
%%  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r] = NL(iflag,nuisvar,method,y_sim,...    %%
%%                 y_meas,sigma,N)                                                    %%
%%  where                                                                             %%
%%   iflag      [input] Estimation ('est') or simulation ('sim')                      %%
%%   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         %%
%%    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       %%
%%    s1:  nuisvar(2) Slope of linear heteroscedastic model                           %%
%%    fi1: nuisvar(3) First-order AR coefficient (0,1)                                %%
%%    fi2: nuisvar(4) Second-order AR coefficient (0,1)                               %%
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
%%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of %%
%%       distribution-adaptive likelihood functions: Generalized and universal        %%
%%       likelihood functions, scoring rules and multi-criteria ranking, Journal of   %%
%%       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             %%
%%       https://www.sciencedirect.com/science/article/pii/S002216942201112X          %%
%%                                                                                    %%
%%  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  %%
%%         2. Measurement error computed from simulated data                          %%
%%                                                                                    %%
%% © Written by Jasper A. Vrugt, Oct. 2013                                            %%
%% DREAM Package                                                                      %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% If number of replicates, N, is not defined
if nargin < 7, N = 1; end; lik_calc = 'cond';

% Initialize return arguments
loglik = -inf*ones(size(y_meas)); [std_e,eps_n,f_eps_n,Y_r] = deal([]);

% Unpack nuisance variables
nuisvar = num2cell(nuisvar); [s0,s1,fi1,fi2] = deal(nuisvar{:});

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
        eps = filter(fi_p,1,e_n);           % AR(2) filter
        eps_n = eps./std_eps;               % Standardized partial residuals
        switch lik_calc
            case 'exact'    % Exact log-likelihood function
                % sum of log-likelihood
                loglik = -n/2 * log(2*pi) - n*log(std_eps) ...
                    + 1/2*log( (1+fi2)^2 * ((1-fi2)^2 - fi1^2) ) ...
                    - (1+fi2)/(2*std_eps^2) * ( (1-fi2)*eps(1)^2 ...
                    - 2*fi1*eps(1)*eps(2) + (1-fi2)*eps(2)^2 ) ...
                    - 1/2*sum(eps_n(3:n).^2) - sum(log(std_e),'omitnan');
                % Return individual log-likelihoods, TO DO
                % loglik = ...
            case 'cond'     % Conditional log-likelihood function
                            % unobserved observations are assumed zero
                % Return individual log-likelihoods
                loglik = - (1/2) * log(2*pi) - log(std_eps) ...
                    - 1/2 * eps_n(1:n).^2 - log(std_e);
                % Sum of log-likelihoods
                % loglik = ...
        end
        % Normal density of eps_n
        if nargout > 3
            f_eps_n = normpdf(eps_n,0,1);
        end

    case 'sim'  % Simulation: Replicates, Y_r, for prediction uncertainty
        eps_n = normrnd(0,1,[n N]);     % nxN matrix samples from N(0,1)
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

% OLD CASE 13
% fi1 = Xp(ii,DREAMPar.d);                               % Extract fi1
% e_n = E(1:n,ii)./std_e;                                 % Standardize residuals
% eps = filter ( [1 -fi1 ] , 1 , e_n );                  % Compute partial residuals
% % Note: the normalized residuals, e_n, do not necessarily
% % have a variance of unity!
% sigma2_eps = 1 - fi1^2; eps_n = eps./sqrt(sigma2_eps); % Standardize partial residuals
% % Now determine log-likelihood (Check Page 4: http://faculty.washington.edu/ezivot/econ584/notes/armaestimation.pdf)
% % ~= sum(log(normpdf(a,0,std_e))); --> cannot compute a(1) as we do not have zeroth observation
% % == sum(log(normpdf(a(2:end),0,std_e(2:end)))) + log(1/(sqrt(2*pi)*sqrt(c)) * exp ( -1/2 * a(1)^2 / c))
% % with c = std_e(1)^2/(1 - rho^2);
% logL_Xp(ii,1) = - ( n/2 ) * log(2*pi) ...
%     - 1/2*log( std_e(1)^2/(1-fi1^2)) ...
%     - (1-fi1^2)/(2*sigma2_eps) * e_n(1)^2 ...
%     - (n-1)/2*log(sigma2_eps) ...
%     - 1/(2*sigma2_eps) * sum ( eps ( 2 : n ).^2 ) ...
%     - sum(log(std_e)); % For transformation of e to e_n;
% % Note: 1/std(eps) * normpdf(eps_n,0,1) = normpdf(eps,0,std(eps))
% % NOTE: sum(log(normpdf(z,0,std(z)))) == sum(log(1/std(z) * normpdf(z/std(z),0,1)))
% % NOTE: identical to likelihood 53
% % NOTE: for fi1 = 0, logL_Xp is identical to logL_Xp of likelihood 12

%  case 43 % Gaussian likelihood function with homoscedastic/heteroscedastic error and AR-1 model of residuals
%             % Same as likelihood 13, yet, normalization precedes treatment of serial correlation
%             fi1 = Xp(ii,DREAMPar.d);                             % Extract fi1
%             e_n = E(1:n,ii)./std_e(1:n);                          % Standardize residuals first
%             eps = filter ( [1 -fi1 ] , 1 , e_n );                % Partial residuals
%             sigma_eps = sqrt(1-fi1^2);                           % Variance partial residuals
%             eps_n(2:n,1) = eps(2:n,1)./sigma_eps;                 % Standardize partial residuals
%             logL_Xp(ii,1) = - ( n/2 ) * log(2*pi) ...             % Compute likelihood
%                 - 1/2 * sum( log ( std_e.^2 ) ) - 1/2*eps(1)^2 ...
%                 - (n - 1) * log ( sigma_eps ) ...
%                 - 1/2 * sum ( eps_n(2:n).^2 );
%             % NOTE: sum(log(normpdf(z,0,std(z)))) == sum(log(1/std(z) * normpdf(z/std(z),0,1)))
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  case 53 % Matrix implementation of Gaussian likelihood using LU decomposition
%             fi1 = Xp(ii,DREAMPar.d);                             % Extract fi1
%             e_n = E(:,ii)./std_e;                                % Studentize residuals
%             sigma2_eps = 1 - fi1^2;                              % Variance partial residuals
%             L = diag(ones(n,1)) + diag(-fi1 * ones(n-1,1),-1);   % Lower triangular matrix
%             U = 1/(1-fi1^2) * [ diag(ones(n,1)) + ...
%                 diag(-fi1 * ones(n-1,1),1) ];                    % Upper triangular matrix
%             U(n,n) = 1;
%             logL_Xp(ii,1) = - ( n/2 ) * log(2 * pi) - sum(log(std_e)) ...     % Log-likelihood
%                 - (n-1)/2 * log(sigma2_eps) - 1/2 * e_n' * L * U * e_n;
