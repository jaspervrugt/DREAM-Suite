function pdf = f_SGT(a,lambda,p,q)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine the density of the standardized SGT (zero-mean and unit     %%
%% standard deviation) for lambda, p and q at the values of a            %%
%%                                                                       %%
%%  SYNOPSIS: pdf = f_SGT(a,lambda,p,q)                                  %%
%%   where                                                               %%
%%    a      [input] array of size A-by-B-by...                          %%
%%    lambda [input] skewness of SGT distribution, -1 < lambda < 1       %%
%%    p      [input] controls kurtosis, p > 0                            %%
%%    q      [input] controls kurtosis, q > 2                            %%
%%    pdf    [outpt] an A-by-B-by... array of pdf = f_SGT(a,lambda,p,q)  %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2016                               %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 4
    error('Requires at least four input arguments');
end

% Determine mean and variance of SGT distribution
if p*q < 2 % mu_lpq and kappa_lpq do not admit zero-mean and unit variance
    error('Product of p and q must exceed two');
else % mu_lpq and kappa_lpq admit a zero-mean and unit variance
    % we expect the standardized partial residuals to be distributed 
    % according to the standardized skewed generalized t distribution, 
    % SGT(.|0,1,lambda,p,q)
    kappa_lpq = beta(1/p,q/p) / sqrt((1+3*lambda^2).*...
        beta(1/p,q/p) * beta(3./p,(q-2)/p) - ...
        4*lambda^2 * beta(2/p,(q-1)/p).^2);
    mu_lpq = 2*kappa_lpq * lambda * 1 * beta(2/p,(q-1)/p) / ...
        beta(1/p,q/p);
end

eps_sgt = a + mu_lpq;                                   % Shift SGT distribution
pdf = p ./ (2 * kappa_lpq * 1 * beta(1/p,q/p)) * ...    % SGT density eps_n
    (1 + ( abs(eps_sgt) ./ (kappa_lpq * 1 * ...
    (1 + lambda*sign(eps_sgt))) ).^p ).^(-(q+1)/p);

end