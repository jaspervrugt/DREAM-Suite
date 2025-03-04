function [loglik,std_e,eps_n,f_eps_n,Y_r,y_bc] = GL(iflag,nuisvar,...
    method,y_sim,y_meas,sigma,N)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Generalized likelihood function: correlated, heteroscedastic and                   %%
%%  non-Gaussian errors                                                               %%
%%  1. Partial standardized residuals are expected to follow a skew exponential       %%
%%     power (SEP) distribution                                                       %%
%%  2. Serial correlation described using an autoregressive model - up to order 2     %%
%%                                                                                    %%
%% This GL function treats serial correlation before standardization (= not ideal)    %%
%% This work has been published in Schoups and Vrugt (2010)                           %%
%%                                                                                    %%
%%  SYNOPSIS: [loglik,std_e,eps_n,f_eps_n,Y_r,y_bc] = GL(iflag,nuisvar,method,...     %%
%%                 y_sim,y_meas,N)                                                    %%
%%  where                                                                             %%
%%   iflag      [input] Estimation ('est') or simulation ('sim')                      %%
%%   nuisvar    [input] Column vector of nuisance variables (fixed/estimated)         %%
%%    s0:  nuisvar(1) Intercept of linear heteroscedastic model                       %%
%%    s1:  nuisvar(2) Slope of linear heteroscedastic model                           %%
%%    ba:  nuisvar(3) Kurtosis (-1: uniform, 0: normal; 1: Laplace)                   %%
%%    xi:  nuisvar(4) Skewness (1: symmetric; <1: negative skew; >1: positive skew)   %%
%%    mu1: nuisvar(5) Bias correction parameter                                       %%
%%    fi1: nuisvar(6) First-order AR coefficient (0,1)                                %%
%%    fi2: nuisvar(7) Second-order AR coefficient (0,1)                               %%
%%    fi3: nuisvar(8) Third-order AR coefficient (0,1)                                %%
%%    fi4: nuisvar(9) Fourth-order AR coefficient (0,1)                               %%
%%    K:   nuisvar(10) Box-Cox transformation parameter (skewness)                    %%
%%    lba: nuisvar(11) Box-Cox transformation parameter (heteroscedasticity)          %%
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
%%   y_bc       [outpt] n x 1 vector of bias corrected simulated values               %%
%%                                                                                    %%
%%  Reference:                                                                        %%
%%  Schoups, G., and J. A. Vrugt (2010), A formal likelihood function for parameter   %%
%%       and predictive inference of hydrologic models with correlated,               %%
%%       heteroscedastic, and non-Gaussian errors, Water Resources Research, 46,      %%
%%       W10531, doi:10.1029/2009WR008933                                             %%
%%                                                                                    %%
%%  Notes: 1. This is a conditional likelihood function: y_(-1) and y_0 assumed zero  %%
%%         2. Measurement error computed from simulated data                          %%
%%         3. This function is OBSOLETE/ERRONEOUS, please use GL_plus                 %%
%%                                                                                    %%
%% Written by Gerrit Schoups (TU Delft) and modified by Jasper A. Vrugt               %%
%% DREAM Package                                                                      %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% To avoid the model to become inactive and let the AR model do the job, we 
% can impose a constraint that simulated discharge of model should be 
% within a certain distance to the measured value. This avoids the model to 
% deactivate with AR inference.

% If number of replicates, N, is not defined
if nargin < 6, N = 1; end

% Initialize return arguments
loglik = -inf*ones(size(y_meas)); [std_e,eps_n,f_eps_n,Y_r] = deal([]);

% Unpack nuisance variables
nuisvar = num2cell(nuisvar); 
[s0,s1,ba,xi,mu1,fi1,fi2,fi3,fi4,K,lba] = deal(nuisvar{:}); 

n = size(y_sim,1);                      % # entries of model simulation
y_bc = y_sim.*min(10,exp(mu1.*y_sim));  % Bias-corrected simulation
y_bc2 = ((y_bc+K).^lba - 1)/lba;        % Box-Cox bias corrected simulation
e_bc = ((y_meas+K).^lba - 1)/lba - ...  % Box-Cox raw residuals
    y_bc2;

% Now we have to studentize/standardize the residuals
switch method
    case 0 % CONSTANT/NONCONSTANT: specifed by user: s0 and s1 irrelevant 
        std_e = sigma;
        % --> s0 removed, s1 removed 
    case 1 % CONSTANT: s0 not selected, s1 not selected 
        std_e = s0 * ones(n,1);
        % --> s0 fixed, s1 removed
    case 2 % CONSTANT: s0 not selected, s1 selected
        std_e = std(e_bc) * ones(n,1);
        % --> s0 phantom variable, s1 removed
    case 3 % CONSTANT: s0 selected, s1 not selected        
        std_e = s0 * ones(n,1);
        % --> s0 estimated, s1 removed
    case 4 % CONSTANT: s0 selected, s1 is selected
        std_e = s0 * ones(n,1);
        % --> s0 estimated, s1 removed
    case 5 % NONCSTNT: s0 not selected, s1 not selected
        [std_e,~,exitflag] = s1_phantom(s0,e_bc,y_bc);
        if exitflag ~= 1, return; end
        % --> s0 fixed, s1 phantom variable
    case 6 % NONCSTNT: s0 not selected, s1 selected       
        std_e = s0 + s1 * y_bc;
        % --> s0 fixed, s1 estimated
    case 7 % NONCSTNT: s0 selected, s1 not selected        
        [std_e,~,exitflag] = s1_phantom(s0,e_bc,y_bc);
        if exitflag ~= 1, return; end
        % --> s0 estimated, s1 phantom variable
    case 8 % NONCSTNT: s0 selected, s1 selected
        std_e = s0 + s1 * y_bc;
        % --> s0 estimated, s1 estimated
end
% std_e = max(1e-6,std1.*Y_bc + std0);  % OLD: Compute standard deviation
fi_p = [ 1 -fi1 -fi2 -fi3 -fi4 ];       % Coefficients of AR filter

%% Now check flag
switch iflag
    case 'est'  % Likelihood estimation
        Y_r = nan(n,1);
        eps = filter(fi_p,1,e_bc);                  % nx1 vector partial residuals
        eps_n = eps./std_e;                         % nx1 vector studentized partial residuals
        [f_eps_n,logf_eps_n] = f_SEP(eps_n,ba,xi);  % nx1 vector SEP densities
        loglik = logf_eps_n - log(std_e) + ...      % nx1 vector log-likelihoods
            (lba-1) * log(y_meas+K);
        % loglik = sum(logf_eps_n) - ...            % Total log-likelihood
        %     sum(log(std_e),'omitnan') + ...
        %     (lba-1)*sum(log(y_meas+K),'omitnan');
        % NOTES
        % 1. If beta --> -1, the SEP distribution becomes a uniform
        %    distribution on the interval -1.732, 1.732 with pdf = Wb.
        %    Any point outside this interval receives a pdf of zero and
        %    as a result the log(pdf) equals -infinity. This creates 
        %    problems in the default formulation below
        %
        % loglik = n*log(2*sig_xi*Wb) - n*log(xi+1/xi) ...      
        %     - sum(log(std_e)) + (lambda-1)*sum(log(Y_meas+K),'omitnan') ...
        %     - sum(Cb * abs(a_skew).^(2/(1+beta)),'omitnan');
        %    where
        % a_skew = (mu_xi + sig_xi.*eps_n) ./ ...   
        %     (xi.^sign(mu_xi + sig_xi.*eps_n)); 
        % 
        %    With a pdf of zero, the first factor should not equal n, 
        %    but n-1. The default formulation returns -inf if a single 
        %    pdf is zero. That is fine; the normalization constant is 
        %    just not 100% right. This will hardly matter in most cases
        %    as -inf is -inf anyway. 
        %    One could think of using "realmin" instead of 0 for pdf. 

    case 'sim'  % Simulation: Replicates, Y_r, for prediction uncertainty
        eps_n = SEPrnd(ba,xi,[n N]);        % nxN matrix samples from SEP(0,1,beta,xi)
        eps = bsxfun(@times,std_e,eps_n);   % Step 6: nxN matrix partial Box-Cox residuals
        e_r = filter(1,fi_p,eps);           % Step 7: nxN matrix raw Box-Cox residuals
        Y_r = (lba*bsxfun(@plus,y_bc2,e_r) + 1) ...
            .^(1/lba) - K;                  % Step 9: Replicate simulations
                                            % Assumption: E[g^-1(Y)] = g^-1(E[Y])
                                            % where g = Box-Cox transformation
end

%% OLD CODE
% % logL = N.*log(Wb*2*sig_xi/(xi+1/xi)) - sum(log(std_e)) - ...
% %            Cb.*(sum(abs(a_xi).^(2./(1+beta))));
% % logL = logL + (lambda-1)*sum ( log ( ObsY + K ) );