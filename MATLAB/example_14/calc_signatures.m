function S = calc_signatures(y,method)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Calculate the summary metrics of measured/simulated discharge record               %%
%%                                                                                    %%
%% SYNOPSIS: S = calc_signatures(y,method)                                            %%
%% where                                                                              %%
%%  y           [input] nx1 vector of discharge values                                %%
%%  method      [input] Computational method for Flow Duration Curve (= FDC)          %%
%%   [a] 1          Weibull plotting position                                         %%
%%   [b] ≠ 1        Built-in plotting position                                        %%
%%  S           [outpt] 1x4 vector with summary metrics of streamflow record, y       %%
%%   S(1)           Annual runoff index                                               %%
%%   S(2)           Annual baseflow index                                             %%
%%   S(3)           Air-entry (= parameter 1) of 2-parameter VG-inspired FDC function %%
%%   S(4)           Slope (= parameter 2) of 2-parameter VG-inspired FDC function     %%
%%                  --> FDCFIT toolbox has many more functionalities and options      %%
%%                                                                                    %%
%% Check the following paper                                                          %%
%%   Vrugt, J.A. (2018), FDCFIT: Software toolbox for fitting a large class of flow   %%
%%       duration curves to streamflow data, manual                                   %%
%%   Sadegh, M., J.A. Vrugt, H.V. Gupta, and C. Xu (2016), The soil water             %%
%%       characteristic as new class of closed-form parametric expressions for the    %%
%%       flow duration curve, J. Hydrol., 535, pp. 438–456                            %%
%%                                                                                    %%
%% (c) Written by Jasper A. Vrugt, March 2014                                         %%
%% Los Alamos National Laboratory 			        	              %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 2, method = 1; end

% Store local variables in memory
persistent P ny fi

if isempty(P)           % Load data and define local variables
    daily_data = load('03451500.dly');  % Load the French Broad data
    P = daily_data(731:end,4);          % Measured precip (2y spin-up)
    ny = numel(y);                      % # streamflow data points
    fi = 0.925;                         % Low pass filter (baseflow index)
end
y = y(:);                           % Make sure y is a column vector
S = nan(1,4);                       % Initial S is empty (= row vector)
S(1) = sum(y)/sum(P);               % S1: (Annual) runoff ratio 
yb = nan(ny,1); yb(1) = y(1)/4;     % Initialize baseflow contribution
for j = 2 : ny          % Low pass filter
    yb(j) = min(fi * yb(j-1) + 1/2 * (1 - fi) * ...
        ( y(j) + y(j-1) ) , y(j) );
end
S(2) = sum(yb)/sum(y);              % S2: (Annual) baseflow index
[Y_s,E_p] = calc_FDC(y,ny,method);  % Sort flows & calc exceedance probs.
x0 = [1 1];                         % Initial parameter values Nelder-Mead
S(3:4) = fminsearch(@(x) ...        % Fit the FDC using van Genuchten WRF
    VG_Ret(x,Y_s,E_p,ny),x0);       % --> see Sadegh et al., JofH, 2015

end
% <><><><><><><><><><><><> End of primary function <><><><><><><><><><><><>

% <><><><><><><><><><><><><> Secondary functions <><><><><><><><><><><><><>
% SUBROUTINE 1
function [Y_s,E_p] = calc_FDC(Y,nY,method)
% Calculate exceedance probability of each flow level

switch method
    case 1
        Y_s = sort(Y);			        % Sort discharge in ascending order
        E_p = 1 - ((1:nY) - .5)'./nY;   % Weibull exceedance probabilities
    otherwise
        [F,Y_s] = ecdf(Y);              % Built-in empirical CDF
        E_p = 1 - F;                    % Exceedance probabilities
                                        % --> see Vrugt, WRR, 2024
end

end

% <><><><><><><><><><><><><> Secondary functions <><><><><><><><><><><><><>
% SUBROUTINE 2
function [RMSE,E_p] = VG_Ret(pars,Y_s,E_ptrue,nY)
% Now compute exceedance probabilities according to VG WRF and values of 
% pars. Next, compute the RMSE between the actual, E_ptrue, and the VG 
% predicted, E_ppred, exceedance probabilities. The FDCFIT toolbox gives
% access to a much larger class of FDC functions. 

alpha = pars(1); n = pars(2);               % Extract parameter values
E_p = (1 + (alpha * Y_s).^n).^(1/n - 1);	% E_p according to VG WRF
RMSE = sqrt(sum((E_p - E_ptrue).^2)/nY);	% Calculate RMSE

end
