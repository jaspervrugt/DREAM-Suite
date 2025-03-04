function logPR_X = Eval_prior(X,SX,Par_info,Meas_info,options)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function evaluates the prior distribution and returns the        %%
%%  log-prior of the N parameter vectors of X                            %%
%%                                                                       %%
%% SYNOPSIS: logPR_X = Eval_prior(X,SX,DREAMPar,Par_info,Meas_info,...   %%
%%                         options)                                      %%
%%  where                                                                %%
%%   X         [input] N x d matrix of parameter vectors                 %%
%%   SX        [input] m x N matrix of m summary metrics paramtr vectors %%
%%   Par_info  [input] Parameter structure                               %%
%%   Meas_info [input] Measurement data structure                        %%
%%   options   [input] Structure with computational settings/options     %%
%%   logPR_X   [outpt] Nx1 vector with logarithmic prior density of X    %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Feb. 2007                               %%
%% Los Alamos National Laboratory                                        %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

[N,d] = size(X);                        % # parameter vectors & dimensions
logPR_X = zeros(N,1);                   % Initialize log prior

if isfield(Par_info,'prior')            % Evaluate prior and then take log
    switch Par_info.u 
        case 'yes'                      % Univariate case
            PR_X = nan(N,d);            % Initialize prior densities
            for ii = 1 : d              % Eval prior density each parameter
                PR_X(1:N,ii) = ...
                    Par_info.prior{ii}(X(1:N,ii));
            end
            if strcmp(Par_info.pr,'pdf')
                logPR_X = sum( ...      % Prior = pdf --> take log density
                    log(PR_X),2); 
            elseif strcmp(Par_info.pr,'logpdf')
                logPR_X = sum(PR_X,2);  % Prior = logpdf --> no log needed
            end
        case 'no'                       % Multivariate case
            PR_X = nan(N,1);            % Initialize prior densities
            for jj = 1 : N              % Eval joint prior density at once
                PR_X(jj,1) = ...           
                    Par_info.prior(X(jj,1:d));
            end
            if strcmp(Par_info.pr,'pdf')
                logPR_X = log(PR_X);    % Prior = pdf --> take log density
            elseif strcmp(Par_info.pr,'logpdf')
                logPR_X = PR_X;         % Prior = logpdf --> no log needed
            end
    end
end

if strcmp(options.DB,'yes')             % Check if we do diagnostic Bayes
    for ii = 1 : N                      % Distance obs/sim summary metrics
        logPR_X(ii,1) = ...             % Log-prior densty ≠ true logpdf
            min(options.epsilon - ...             
            abs(Meas_info.S(:) - SX(:,ii)));
    end
end