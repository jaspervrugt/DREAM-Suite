function [pdf,logpdf] = f_SEP(a,beta,xi)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Determine the density of the standardized SEP (zero-mean and unit     %%
%% standard deviation) for beta and xi at the values of a                %%
%%                                                                       %%
%%  SYNOPSIS: [pdf,logpdf] = f_SEP(a,beta,xi)                            %%
%%   where                                                               %%
%%    a      [input] array of size A-by-B-by...                          %%
%%    beta   [input] skewness -1 < beta <= 1                             %%
%%    xi     [input] kurtosis xi > 0                                     %%
%%    pdf    [outpt] an A-by-B-by... array of pdf = f_SEP(a,beta,xi)     %%
%%    logpdf [outpt] an A-by-B-by... array of natural log of pdf         %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2014                               %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 3
    error('Requires at least three input arguments');
end

% Compute SEP parameters
A1 = gamma(3*(1+beta)/2);
A2 = gamma((1+beta)/2);
Cb = (A1/A2)^(1/(1+beta));
Wb = sqrt(A1)/((1+beta)*(A2^(1.5)));
M1 = gamma(1+beta)/sqrt(A1*A2);
M2 = 1; mu_xi = M1*(xi-1/xi);
sig_xi = sqrt((M2-M1^2)*(xi^2 + 1/xi^2) + 2*M1^2 - M2);

a_SEP = mu_xi + sig_xi * a;
a_skew = a_SEP ./ xi.^(sign(a_SEP));
% If beta --> -1, then function may return nan
% Note that if beta --> -1, xi does not matter anymore and
% the distribution is synonymous to uniform between -1/(2*pdf),1/(2*pdf)
% With beta --> -1, Cb goes to zero and pdf becomes equal to Wb
norm_cst = max(realmin , (2*sig_xi/(xi + 1/xi)) * Wb); % --> converges to Wb
switch Cb < 1e-100
    case 1 % SEP density is uniform (note: impl. for integral function)
        pdf = realmin * ones(size(a));          % initialize pdf to zero
        ii = abs(a_skew) <= (1/2 * 1/norm_cst); % in this range, pdf = Wb  
        pdf(ii) = norm_cst;                     % pdf equal to Wb
        logpdf = log(pdf);                      % logpdf
    otherwise % SEP density of eps_n (regular expression)
        logpdf = log(norm_cst) ...
            - Cb * abs(a_skew).^(2/(1+beta));
        pdf = exp(logpdf);
end

end

        % a_skew = (mu_xi + sig_xi.*eps_n) ./ ...   
        %     (xi.^sign(mu_xi + sig_xi.*eps_n));
        % loglik = n*log(2*sig_xi*Wb) - n*log(xi+1/xi) ...      
        %     - sum(log(std_e)) + (lambda-1)*sum(log(Y_meas+K),'omitnan') ...
        %     - sum(Cb * abs(a_skew).^(2/(1+beta)),'omitnan');

% f_SEP = (2*sig_xi/(xi + 1/xi)) * Wb * exp(-Cb * abs(a_skew).^(2/(1+beta)));     
