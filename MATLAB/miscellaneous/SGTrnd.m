function R = SGTrnd(mu,sigma,lambda,p,q,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draws samples, R, from the skewed generalized Student's t (SGT)       %%
%% distribution, R ~ SGT(mu,sigma,lambda,p,q), where R is A-by-B-by...   %%
%%                                                                       %%
%%  SYNOPSIS: R = SGTrnd(mu,sigma,lambda,p,q,[A B ...])                  %%
%%   where                                                               %%
%%    mu     [input] mean of SGT distribution                            %%
%%    sigma  [input] std. dev. of SGT distribution                       %%
%%    lambda [input] skewness of SGT distribution, -1 < lambda < 1       %%
%%    p      [input] controls kurtosis, p > 0                            %%
%%    q      [input] controls kurtosis, q > 2                            %%
%%    R      [outpt] an A-by-B-by... array, R ~ SGT(mu,sigma,lambda,p,q) %%
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

if nargin < 5
    error('Requires at least five input arguments');
end

% Extract distribution parameters from input arguments
[err, sizeOut] = sizechk(5,p,q,lambda,mu,sigma,varargin{:});
if err > 0
    error('Size information is inconsistent');
end

% Return NaN for invalid parameter values
if (sigma <= 0 || lambda <= -1 || lambda >= 1 || p <= 0 || q <= 2 || ...
    isnan(mu) || isnan(sigma) || isnan(lambda) || isnan(p) || isnan(q))
    R = nan(sizeOut);
else
    u = rand(sizeOut);
    % Express SGT random numbers as a function of inverse cdf
    R = SGTinv(u,mu,sigma,lambda,p,q);
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, commonSize, numElements] = sizechk(npar,varargin)
% SIZECHK Check for compatible array sizes.
%   [err, commonSize, numElements] = sizechk(npar,A,B,...,M,N,...) or
%   [err, commonSize, numElements] = sizechk(npar,A,B,...,[M,N,...])
%   in effect computes size( A + B + ... + zeros(M,N,...) ), and catches
%   any size mismatches; npar is the number of SGTrnd input arguments

try
    tmp = 0;
    for argnum = 1:npar
        tmp = tmp + varargin{argnum};
    end
    if nargin > npar + 1
        tmp = tmp + zeros(varargin{npar+1:end});
    end
    err = 0;
    commonSize = size(tmp);
    numElements = numel(tmp);
catch
    err = 1;
    commonSize = [];
    numElements = 0;
end

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%