function pdf = f_SST(a,nu,xi,method)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Determine the density of the standardized SST (zero-mean and unit     %%
%% standard deviation) for nu and xi at the values of a                  %%
%%                                                                       %%
%%  SYNOPSIS: pdf = f_SST(a,nu,xi,method)                                %%
%%   where                                                               %%
%%    a      [input] array of size A-by-B-by...                          %%
%%    nu     [input] number of degrees of freedom: nu > 2                %%
%%    xi     [input] kurtosis xi > 0                                     %%
%%    method [input] OPT: calculation method                             %%
%%    pdf    [outpt] an A-by-B-by... array of pdf = f_SST(a,nu,xi)       %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2016                               %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 4, method = 1; end
% method = 1: universal solution that works for all nu
if nargin < 3
    error('Requires at least three input arguments');
end

% Compute first and second moment of SST distribution
switch method
    case 0 % OK for small nu; if nu > 100 or so, M1 and M2 not defined
        M = @(j,nu) gamma((j+1)/2)*gamma((nu-j)/2)*(nu-2)^(j/2) / ...
            (sqrt(pi)*gamma(nu/2));
        M1 = M(1,nu); M2 = M(2,nu);
    case 1 % OK for all nu --> log formulation to avoid numerical overflow
        logM = @(j,nu) gammaln((j+1)/2) + gammaln((nu-j)/2) + ...
            (j/2)*log(nu-2) - 1/2*log(pi) - gammaln(nu/2);
        M1 = exp(logM(1,nu)); M2 = exp(logM(2,nu));
end
% Now compute mu_xi and sig_xi
mu_xi = M1*(xi - 1/xi);     % Matches Scharnagl paper
        % mu_xi = gamma((nu-1)/2)*sqrt(nu-2)*(xi-1/xi)/(sqrt(pi)*gamma(nu/2))
        % log_mu_xi = gammaln((nu-1)/2) + 1/2*log(nu-2) ...
        % + log(xi-1/xi) - 1/2*log(pi) - gammaln(nu/2)
sig_xi = sqrt((M2-M1^2)*(xi^2 + xi^-2) + 2*M1^2 - M2);
        % = sqrt( -mu_xi^2+ xi^2 + 1/xi^2 - 1)

% Standardize a
a_sst = (mu_xi + sig_xi .* a) ./ ...
                (xi.^sign(mu_xi + sig_xi .* a));  

% We resort to non-standardized Student t-distribution (WIKI) -->
% var(X) = sigma2*(nu/(nu-2)); for nu > 2
% We want var(X) = 1, thus, sigma2 = (nu-2)/nu;

% Standardized Student's t density with skew
switch method % Compute log-likelihood and density
    case 0  % OK for small degrees of freedom, nu
        pdf = 2*sig_xi*gamma((nu+1)/2) ./ ...
            ((xi + 1/xi)*gamma(nu/2)*sqrt(pi*(nu-2))) * ...
            ( 1 + 1/(nu-2) * a_sst.^2).^(-(nu+1)/2);
    case 1  % OK for all degrees of freedom, nu > 2
        % gammaln(x) function avoids gamma(x) to inf for large x
        log_pdf = log(2*sig_xi) + gammaln((nu+1)/2) - log(xi + 1/xi) ...
            - gammaln(nu/2) - 1/2*log(pi*(nu-2)) ...
            - (nu+1)/2 * log ( 1 + 1/(nu-2) * a_sst.^2 );
        pdf = exp(log_pdf);
end

end

% trapz(a,pdf) % must integrate to one