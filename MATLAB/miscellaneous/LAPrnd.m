function R = LAPrnd(mu,sigma,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draws samples, R, from the Laplacian (LAP) distribution with mean, mu %%
%% and standard dev., sigma, R ~ LAP(mu,sigma), where R is A-by-B-by...  %%
%%                                                                       %%
%%  SYNOPSIS: R = LAPrnd(mu,sigma,[A B ...])                             %%
%%   where                                                               %%
%%    mu     [input] mean of LAP distribution             mu in R        %%
%%    sigma  [input] std. dev. of LAP distribution        sigma > 0      %%
%%    R      [outpt] an A-by-B-by... array of draws, R ~ LAP(mu,sigma)   %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2020                               %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('LAPrnd:Requires at least two input arguments');
end

% Extract distribution parameters from input arguments
[err, sizeOut] = sizechk(2,mu,sigma,varargin{:});
if err > 0
    error('LAPrnd:Size information is inconsistent');
end

% Return NaN for invalid parameter values
if (sigma <= 0 || isnan(sigma))
    R = nan(sizeOut);
else
    u = rand(sizeOut) - .5;                     % Generate Laplacian noise
    b = sigma/sqrt(2);                          % Compute scale parameter
    R = mu - b * sign(u).* log(1 - 2*abs(u));   % Generate draws
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, commonSize, numElements] = sizechk(npar,varargin)
% SIZECHK Check for compatible array sizes.
%   [err, commonSize, numElements] = sizechk(npar,A,B,...,M,N,...) or
%   [err, commonSize, numElements] = sizechk(npar,A,B,...,[M,N,...])
%   in effect computes size( A + B + ... + zeros(M,N,...) ), and catches
%   any size mismatches; npar is the number of LAPrnd input arguments

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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%