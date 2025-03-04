function X = SSTinv(P,nu,xi)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inverse of standardized Skewed Student t (SST) distribution.          %%
%% X = SSTinv(P,nu,xi) returns the inverse of the standardized SST's cdf %%
%% with nu and xi, at the values in P                                    %%
%%                                                                       %%
%%  SYNOPSIS: X = SSTinv(P,nu,xi)                                        %%
%%   where                                                               %%
%%    P      [input] array of size A-by-1-by...                          %%
%%    nu     [input] degrees of freedom q > 2                            %%
%%    xi     [input] kurtosis xi > 0                                     %%
%%    X      [outpt] an A-by-1-by... array with inverse of SST cdf       %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2016                               %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    error('SSTinv:Requires at least three input arguments');
end
if nu <= 2 || xi <= 0  
    warning('SSTinv: nu > 2 and xi > 0'); X = nan(size(P)); return
end
if isinf(nu)
    warning('SSTinv:Cannot specify inf for nu: Setting value of nu to 1e10');
    nu = 1e10;
end

calc_method = 1;    % [0]: brute force
                    % [1]: numerical solution
P = P(:); nP = numel(P); X = nan(nP,1);     % How many elements of P?

% Now return X values at given quantiles, P       
switch calc_method
    case 0 % brute force
        x = -50:0.001:50;                   % Cdf at many values of x
                                            % --> OK if mu = 0 and sigma = 1;
        pdf_x = f_SST(x,nu,xi);             % Compute PDF at x
        cdf_x = cumsum(pdf_x)/sum(pdf_x);   % cdf at x
        ii = [ 1 find(diff(cdf_x)>0) + 1]'; % Check duplicate points
        x = x(ii); cdf_x = cdf_x(ii);       % Remove - interp1 cannot handle
        X = interp1(cdf_x,x,P);             % Inverse CDF through interp1
    case 1 % numerical solution from pdf
        for i = 1:nP
            f = @(x) integral(@(x) ...      % Setup as root finding problem
                f_SST(x,nu,xi),-inf,x) ...
                - P(i);
            X(i) = fzero(f,0);              
        end
end