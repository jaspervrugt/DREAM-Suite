function loglik = BMA_lik(par,plugin)
% This function computes the BMA log likelihood for weights and sigma's
% Recommend using the MODELAVG toolbox instead!

D = plugin.BMA.D; y = plugin.BMA.y; % Ensemble forecasts & verifying data
[n,K] = size(D);                    % # forecasts, # ensmble mmbers
w = par(1:K);                       % Unpack weights
w = w./sum(w,2);                    % Normalize weights: quick fix, but not recommended
                                    % Please use MODELAVG package

switch plugin.BMA.VAR   % VARIANCE OPTION 
                        % nxK matrix standard deviation forecasts
    case {'1'} % 1: common constant variance
        sigma = par(K+1) * ones(n,K);                                    
    case {'2'} % 2: individual constant variance
        sigma = bsxfun(@times,par(K+1:2*K),ones(n,K));                               
    case {'3'} % 3: common non-constant variance
        c = par(K+1); sigma = c * D;                               
    case {'4'} % 4: individual non-constant variance
        c = par(K+1:2*K); sigma = bsxfun(@times,c,D);
end
sigma = max(sigma,eps);             % sigma >= 2.22e-16

switch plugin.BMA.PDF   % CONDITIONAL DISTRIBUTION 
                        % nxK matrices A and B
    case 'normal'       % Gaussian with µ = D and standard deviation sigma
        A = D; B = sigma; 
    case 'gamma'    % Gamma with shape A and scale B  
        mu = abs(D); var = sigma.^2; A = mu.^2./var; B = var./mu;
end
Y = repmat(y,1,K);                  % Make K copies verifying data
L = pdf(plugin.BMA.PDF,Y,A,B);      % nxK matrix likelihoods forecasts
lik = L*w' + realmin;                % nx1 vector likelihoods BMA model
loglik = sum(log(lik));             % log-likelihood (= scalar)

end
