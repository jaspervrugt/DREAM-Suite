function R = SEPrnd(beta,xi,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draws samples, R, from the standardized skewed generalized            %%
%% exponential power (SEP) distribution, R ~ SEP(0,1,beta,xi), where R   %%
%% is A-by-B-by...                                                       %%
%%                                                                       %%
%%  SYNOPSIS: R = SEPrnd(0,1,beta,xi,[A B ...])                          %%
%%   where                                                               %%
%%    beta   [input] skewness -1 < beta <= 1                             %%
%%    xi     [input] kurtosis xi > 0                                     %%
%%    R      [outpt] an A-by-B-by... array, R ~ SEP(0,1,beta,xi)         %%
%%                                                                       %%
%% Written by Gerrit Schoups, modified by Jasper A. Vrugt                %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('Requires at least two input arguments');
end

% Extract distribution parameters from input arguments
[err, sizeOut] = sizechk(2,beta,xi,varargin{:});
if err > 0
    error('Size information is inconsistent');
end

% Return NaN for invalid parameter values
if (beta <= -1 || beta > 1 || xi <= 0 || isnan(beta) || isnan(xi))
    R = nan(sizeOut);
else
    % Compute SST variables so that mu = 0 and sigma = 1
    A1 = gamma(3*(1+beta)/2);
    A2 = gamma((1+beta)/2);
    M1 = gamma(1+beta)/sqrt(A1*A2);
    M2 = 1; 
    mu_xi = M1*(xi-1/xi);
    sig_xi = sqrt((M2-M1^2)*(xi^2 + 1/xi^2) + 2*M1^2 - M2);
%     fr = @(r,xi) (xi.^(r+1) + (-1)^r/xi^(r+1))/(xi + 1/xi);
%     mu_xi = M1*fr(1,xi)
%     sig_xi = sqrt((M2-M1^2)*fr(2,xi) + 2*M1^2 - M2)
    p = 2/(1+beta);
    grnd = gamrnd(1/p,ones(sizeOut));                   % Step 1: Generate n random variates from gamma distribution with
                                                        %         shape parameter 1/p and scale parameter 1
    signrnd = sign(rand(sizeOut)-0.5);                  % Step 2: Generate n random signs (+1 or -1) with equal probability
    EP_rnd = signrnd.*(abs(grnd).^(1/p)).* ...          % Step 3: Compute n random variates from EP(0,1,beta)
        sqrt(gamma(1/p))./sqrt(gamma(3/p));
    w = xi/(xi+1/xi); signrndw = sign(rand(sizeOut)-w); % Step 4: Generate n random signs (+1 or -1) with probability 1-w and w
    SEP_rnd = -signrndw.*abs(EP_rnd)./(xi.^signrndw);   % Step 5: Compute n random variates from SEP(mu_xi,sig_xi,xi,beta)
    R = (SEP_rnd-mu_xi)./sig_xi;       
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
