function X = LAPinv(P,mu,sigma)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inverse of Laplace distribution (LAP) cumulative density function     %%
%% (cdf). X = LAPinv(P,mu,b) returns the inverse of LAP's cdf with mu    %%
%% and b, at the values in P                                             %%
%%                                                                       %%
%%  SYNOPSIS: R = LAPinv(P,mu,sigma)                                     %%
%%   where                                                               %%
%%    P      [input] array of size A-by-1-by...                          %%
%%    mu     [input] mean of Laplace distribution         mu in R        %%
%%    sigma  [input] std. dev. of LAP distribution        sigma > 0      %%
%%    X      [outpt] an A-by-1-by... array with inverse of LAP cdf       %%
%%                                                                       %%
%% © Written by Jasper A. Vrugt, Dec. 2020                               %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    error('LAPinv:Requires at least three input arguments');
end
if sigma <= 0  
    warning('LAPinv: sigma > 0'); X = nan(size(P)); return
end

calc_method = 2;    % [0]: brute force
                    % [1]: numerical solution
                    % [2]: analytic solution
P = P(:); nP = numel(P); X = nan(nP,1);     % How many elements of P?
b = sigma / sqrt(2);
F_lap = @(x) .5 + .5 * sign(x-mu) ...       % Define CDF Laplace distribution
    .* (1 - exp(- abs(x-mu)/b));    

% Now return X values at given quantiles, P       
switch calc_method
    case 0 % brute force
        x = -10:0.001:10;                   % Compute cdf at many values of x
                                            % --> OK if mu = 0 and b = 1;
        cdf_x = F_lap(x);                   % Compute cdf at x
        ii = [ 1 find(diff(cdf_x)>0) + 1]'; % Check duplicate points
        x = x(ii); cdf_x = cdf_x(ii);       % Remove - interp1 cannot handle
        X = interp1(cdf_x,x,P);             % Inverse CDF through interp1
    case 1 % more elegant
        for i = 1:nP
            f = @(x) F_lap(x) - P(i);       % Setup as root finding problem
            X(i) = fzero(f,0);              % Solve for zero point
        end
    case 2 % most elegant: quantile function
        Fi_lap = @(x) mu - b*sign(x-.5) ...
            .* log(1-2*abs(x-.5));          % Quantile function LAP
        X = Fi_lap(P);                      % Return quantiles at P
end

end

% % F_lap = @(x) .5*(1 + sign(x-mu) ...        % Old formulation - why the square root of 2?
% %     .* (1 - exp(- abs(x-mu)/(b/sqrt(2)) )));
