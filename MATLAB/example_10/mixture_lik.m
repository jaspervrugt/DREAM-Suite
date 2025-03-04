function loglik = mixture_lik(x)
% d-variate normal mixture log-likelihood
% Target distribution is w1*N(µ1,Σ1) + w2*N(µ2,Σ2) with µ1 = -5, µ2 = 5,
%    Σ1/Σ2 = identity matrix and w1 = 1/3 and w2 = 2/3
% with informative prior

comp_method = 3;    % 1: Exact using built-in functions (slow)
                    % 2: Fast but without normalization constant [= OK] 
                    % 3: Refined/fast [= some protection against underflow] 
% What is underflow? L(x) → zero when x is far removed from target mean
%                         → log(L(x)) goes to -∞.
% For example, if 
% x = [11.8 17.0 19.8 18.1 24.7 23.5 11.0 4.8 29.4 0.02]
% then comp_method = 1 yield loglik = -Inf, whereas comp_method = 3 yields
% loglik = -993.3101. 

% Store local variables in memory
persistent d k lik mu log_Z w1 C Cc
% Determine constant variables
if isempty(k)
    % How many dimensions?
    d = size(x,2);
    % How many Gaussians?
    k = 2;
    % Determine weight of each respective Gaussian
    w1 = 1/3; w2 = 1 - w1;
    log_prior = log([w1 w2]);
    % Mean of both modes
    mu = [ -5*ones(1,d) ; 5*ones(1,d) ];
    % Covariance matrix
    C = eye(d);
    % Cholesky factorization
    Cc = chol(C); logDetSigma = 2 * sum ( log ( diag( Cc ) ) );
    % Calculate ~normalization constant of each normal
    log_Z = -0.5 * logDetSigma + log_prior - d*log(2*pi)/2;
end
ll_max = 0;

switch comp_method
    case 1  % Built-in function: no protection against underflow
        lik = w1 * mvnpdf(x,mu(1,1:d),C) + (1-w1) * mvnpdf(x,mu(2,1:d),C);
    case 2  % Faster but no protection against underflow
        lik = w1 * exp(-(x+5)*(x+5)') + (1-w1) * exp(-(x-5)*(x-5)');
    case 3  % Some protection against underflow 
        ll_j = nan(1,k);
        % Loop over each component of mixture
        for j = 1:k
            % Calculate log_likelihood of jth component
            ll_j(j) = -.5*sum((bsxfun(@minus,x,mu(j,:))/Cc).^2,2) + log_Z(j);
        end
        % Now determine max of log-likelihoods
        ll_max = max(ll_j); 
        % Compute likelihood (minus maxll to avoid underflow)
        lik = sum(exp(bsxfun(@minus,ll_j,ll_max)),2);
        % p(x_i) is \sum_j \alpha_j P(x_i| \theta_j)/ exp(maxll(i))
end
% Log-likelihood of normal mixture
loglik = log(lik) + ll_max;

end