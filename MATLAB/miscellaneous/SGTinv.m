function X = SGTinv(P,mu,sigma,lambda,p,q)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantile function of the Skewed Generalized Student's t (SGT)         %%
%% distribution (inverse CDF). X = SGTinv(P,mu,sigma,lambda,p,q) returns %%
%% the inverse of SGT's cdf with mu, sigma, lambda, p and q, at the      %%
%% values in P                                                           %%   
%%                                                                       %%
%%  SYNOPSIS: X = SGTinv(P,mu,sigma,lambda,p,q)                          %%
%%   where                                                               %%
%%    P      [input] array of size A-by-B-by...                          %%
%%    mu     [input] mean of SGT distribution                            %%
%%    sigma  [input] std. dev. of SGT distribution                       %%
%%    lambda [input] skewness of SGT distribution, -1 < lambda < 1       %%
%%    p      [input] controls kurtosis, p > 0                            %%
%%    q      [input] controls kurtosis, q > 2                            %%
%%    X      [outpt] an A-by-B-by... array with inverse of SGT cdf       %%
%%                                                                       %%
%% If p = 2 this function becomes equivalent to the Standardized         %%
%% Student t-distribution with q degrees of freedom. If q = inf, then    %%
%% this function is equal to the skewed generalized exponential          %%
%% distribution, GED(p). When p = 2 and q = inf, then we yield the       %%
%% skewed normal distribution                                            %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2020 based on existing code        %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    error('Requires at least six input arguments');
end
% Initialize z to zero
if isa(P,'single') || isa(p,'single') || isa(q,'single') || ...
        isa(lambda,'single') || isa(sigma,'single')
    X = zeros(size(P),'single');
else
    X = zeros(size(P));
end

% For q > 200, the GST becomes equivalent to GED
q_max = 200;
% The inverse cdf of 0 is -Inf, and the inverse cdf of 1 is Inf
X(P == 1) = Inf; X(P == 0) = -Inf;
% Find values between zero and one
ii = find(P > 0 & P < 1);
% Special case when lambda = 0 and x = 0.5, SGTinv(1/2) = 0
k = find(P == 1/2 & lambda == 0);
if any(k)
    X(k) = 0; ii(ii == k) = [];
end
% Check whether SGT or GED
if q <= q_max % SGT distribution
    theta = beta(1/p,q/p) / sqrt((1+3*lambda^2) * ...
        beta(1/p,q/p) * beta(3/p,(q-2)/p) - ...
        4*lambda^2 * beta(2/p,(q-1)/p)^2);
    mu = mu - 2*theta * lambda * sigma * beta(2/p,(q-1)/p) / ...
        beta(1/p,q/p);
    w = betaincinv((2*P(ii) - (1-lambda)) ./ (lambda+sign(P(ii) - ...
        0.5*(1-lambda))) , 1/p , q/p);
    X(ii) = theta * sigma * (sign(P(ii) - 0.5*(1-lambda)) + lambda)./...
        (w.^(-1)-1).^(1/p) + mu;
else % GED distribution
    theta = gamma(1/p) / sqrt((1+3*lambda^2) *...
        gamma(1/p) * gamma(3./p) - 4*lambda^2 *...
        gamma(2/p)^2);
    mu = mu - 2*theta * lambda * sigma * gamma(2/p) / ...
        gamma(1/p);
    w = gammaincinv((2*P(ii) - (1-lambda)) ./ (lambda+sign(P(ii) - ...
        0.5*(1-lambda))) , 1/p);
    X(ii) = theta * sigma * (lambda + sign(P(ii) - 0.5*(1-lambda))) .* ...
        w.^(1/p) + mu;
end

end