function X = SEPinv(P,beta,xi)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantile function of the standardized skewed exponential power (SEP)  %%
%% distribution. X = SEPinv(P,beta,xi) returns the inverse of the        %%
%% standardized SEP's CDF with beta and xi, at the values in P           %%
%%                                                                       %%
%%  SYNOPSIS: X = SEPinv(P,beta,xi)                                      %%
%%   where                                                               %%
%%    P      [input] array of size A-by-1-by...                          %%
%%    beta   [input] skewness -1 < beta <= 1                             %%
%%    xi     [input] kurtosis xi > 0                                     %%
%%    X      [outpt] an A-by-1-by... array with inverse of SEP CDF       %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2020                               %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    error('SEPinv:Requires at least three input arguments');
end
if beta <= -1 || beta > 1 || xi <= 0  
    warning('SEPinv:-1 < beta <= 1 and xi > 0'); X = nan(size(P)); return
end

calc_method = 1;    % [0]: brute force
                    % [1]: numerical solution using PDF integration
                    % [2]: numerical solution using CDF
P = P(:); nP = numel(P); X = nan(nP,1);     % How many elements of P?

% Now return X values at given quantiles, P       
switch calc_method
    case 0 % brute force
        x = -50:0.001:50;                   % Cdf at many values of x
                                            % --> OK if mu = 0 and sigma = 1;
        pdf_x = f_SEP(x,beta,xi);           % Compute PDF at x
        cdf_x = cumsum(pdf_x)/sum(pdf_x);   % cdf at x
        ii = [ 1 find(diff(cdf_x)>0) + 1]'; % Check duplicate points
        x = x(ii); cdf_x = cdf_x(ii);       % Remove - interp1 cannot handle
        X = interp1(cdf_x,x,P);             % Inverse CDF through interp1
    case 1 % more elegant
        for i = 1:nP
            f = @(x) integral(@(x) ...      % Setup as root finding problem
                f_SEP(x,beta,xi),-inf,x) ...
                - P(i);
            X(i) = fzero(f,0);              % Solve for zero point
        end
end

end