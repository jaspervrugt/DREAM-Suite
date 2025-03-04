function R = SSTrnd(nu,xi,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draws samples, R, from the standardized skewed Student's t (SST)      %%
%% distribution, R ~ SST(nu,xi), where R is A-by-B-by...                 %%
%%                                                                       %%
%%  SYNOPSIS: R = SSTrnd(nu,xi,[A B ...])                                %%
%%   where                                                               %%
%%    nu     [input] degrees of freedom q > 2                            %%
%%    xi     [input] kurtosis xi > 0                                     %%
%%    R      [outpt] an A-by-B-by... array, R ~ SST(nu,xi,[A B])         %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2016                               %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

calc_method = 1; % Correct solution

if nargin < 2
    error('Requires at least two input arguments');
end

% Extract distribution parameters from input arguments
[err, sizeOut] = sizechk(2,nu,xi,varargin{:});
if err > 0
    error('Size information is inconsistent');
end

% Return NaN for invalid parameter values
if (nu <= 2 || xi <= 0 || isnan(nu) || isnan(xi))
    R = nan(sizeOut);
else
    % Compute first and second moment of SST distribution
    switch calc_method
        case 0 % OK for small nu; if nu > 100 or so, M1 and M2 not defined
            M = @(j,nu) gamma((j+1)/2)*gamma((nu-j)/2)*(nu-2)^(j/2) / ...
                (sqrt(pi)*gamma(nu/2));
            M1 = M(1,nu); M2 = M(2,nu);
        case 1 % OK for all nu --> log formulation avoids numerical overflow
            logM = @(j,nu) gammaln((j+1)/2) + gammaln((nu-j)/2) + ...
                (j/2)*log(nu-2) - 1/2*log(pi) - gammaln(nu/2);
            M1 = exp(logM(1,nu)); M2 = exp(logM(2,nu));
    end
    % Now compute mu_xi and sig_xi
    mu_xi = M1*(xi - 1/xi);                                         
    % Matches Scharnagl paper
    % mu_xi = gamma((nu-1)/2)*sqrt(nu-2)*(xi-1/xi)/(sqrt(pi)*gamma(nu/2))
    sig_xi = sqrt((M2-M1^2)*(xi^2 + xi^-2) + 2*M1^2 - M2);
    % = sqrt( -mu_xi^2 + xi^2 + 1/xi^2 - 1)
    % = sqrt( -mu_xi^2 + (xi^3 + 1/xi^3)/(xi + 1/xi) )

    switch calc_method
        case 0 % Brute force - we do not use default functions
            a = (-10:0.001:10)'; f_a = f_SST(a,nu,xi);
            cdf_sst = cumsum(f_a)/sum(f_a);
            ii = diff(cdf_sst) > 0; cdf_sst = cdf_sst(ii); a = a(ii);
            u_rnd = rand(sizeOut); R = interp1(cdf_sst,a,u_rnd);
        case 1 % More elegant solution - using statistics
            % We resort to non-standardized Student t-distribution 
            % WIKI: var(X) = sigma2*(nu/(nu-2)); for nu > 2
            % We want var(X) = 1, thus, sigma2 = (nu-2)/nu; 
            % --> sigma = sqrt((nu-2)/nu)
            ST_rnd = random('tlocationscale',...    % Step 1: Sample variates from standardized Student
                0,sqrt((nu-2)/nu),nu,sizeOut);      %         t distribution with nu degrees of freedom
                                                    %         --> matches f_ST(a|nu,1) in paper
            w = xi/(xi+1/xi);                       % Step 2: Create n random signs (+1 or -1) with prob. 1-w and w
            signrndw = sign(rand(sizeOut)-w);
            SST_rnd = -signrndw.* ...               % Step 3: Compute n random variates from SST(nu,xi)
                abs(ST_rnd)./(xi.^signrndw);
            R = (SST_rnd-mu_xi)./sig_xi;            % Step 4: Normalize to obtain n standardized variates from SST(nu,xi)
    end

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