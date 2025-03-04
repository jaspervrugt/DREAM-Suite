function [I,i,flag] = Haverkamp_I_patch(eta,plugin,rtol,kmax,itol)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Function solves for the I(t) relationship of Haverkamp using Newton's %%
%% method                                                                %%
%%  SYNOPSIS: I = Haverkamp_I_patch(eta,plugin)                          %%
%%            I = Haverkamp_I_patch(eta,plugin,t,rtol)                   %%
%%            I = Haverkamp_I_patch(eta,plugin,t,rtol,kmax)              %%
%%            I = Haverkamp_I_patch(eta,plugin,t,rtol,kmax,itol)         %%
%% where                                                                 %%
%%  eta     [input] 4x1 vector with Haverkamp parameter values           %%
%%                      S   [cm/h^0.5]                                   %%
%%                      Ks  [cm/h]                                       %%
%%                      ß   [-]                                          %%
%%                      Ki  [cm/h]                                       %%
%%  plugin  [input] structure with nx1 vector time in h                  %% 
%%  rtol    [input] OPT: tolerance on function value root   (def: 1e-12) %%
%%  kmax    [input] OPT: maximum # Newton iterations        (def: 20)    %%
%%  itol    [input] OPT: tolerance constant rate assumption (def: 1e-10) %%
%%  I       [outpt] cumulative infiltration, I(t) in cm                  %%
%%  i       [outpt] infiltration rate, i(t) in cm/h                      %%
%%  flag    [outpt] exit flag: [1] exact [2] approximate                 %%
%%                                                                       %%
%%  LITERATURE                                                           %%
%%  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and  %%
%%      H. Vereecken (2023), The time validity of Philip's two-term      %%
%%      infiltration equation: An elusive theoretical quantity? Vadose   %%
%%      Zone Journal, e20309, pp. 1-25,                                  %%
%%          https://doi.org/10.1002/vzj2.20309                           %%
%%  Vrugt, J.A. and y. Gao (2022), On the three-parameter infiltration   %%
%%      equation of Parlange et al. (1982): Numerical solution,          %%
%%      experimental design, and parameter estimation, Vadose Zone       %%
%%      Journal, 21:e20167, pp. 1-25,                                    %%
%%          https://doi.org/10.1002/vzj2.20167                           %%
%%                                                                       %%
%%  © Written by Jasper A. Vrugt, Jan. 2019                              %%
%%  University of California Irvine                                      %%
%%                                                                       %%
%%  © Written by Jasper A. Vrugt, Jan. 2019                              %%
%%  University of California Irvine                                      %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if nargin < 5, itol = 1e-10; end
if nargin < 4, kmax = 20; end
if nargin < 3, rtol = 1e-12; end

%% Unpack parameter values
S = eta(1);             % Units of cm/h^(1/2)
Ks = eta(2);            % Units of cm/h
switch plugin.model_setup
    case 1
        B = eta(3);    % Unitless
        Ki = 0;        % Units of cm/h
    case 2
        Ki = eta(3);   % Units of cm/h
        B = eta(4);    % Unitless
end
%% Initialization
dK = Ks - Ki; xi = dK/S^2; t = plugin.t(:);         % Compute dK, xi and t column vector
n = numel(t); [I,i] = deal(nan(n,1)); flag = 1;     % Initialize cum. inf., rate/flag
if t(1) == 0, I(1) = 0; j = 2; else, j = 1; end     % Set I(1) to zero if t(1) is zero
dt = [ t(1) ; diff(t) ];                            % Compute delta time
r = @(I,t) I-1/(2*xi)*log(exp(2*xi*B*(I-Ki*t))/B + ...
    (B-1)/B) - (dK*(1-B)+Ki)*t;                     % Residual function
dr = @(I,t) (1-(B*exp(2*xi*B*(I-Ki*t)))./ ...
    (exp(2*xi*B*(I-Ki*t))+B-1));                    % Derivative function
%% Dynamic part
I_up = S*sqrt(t) + t*Ks; I_low = max(I_up(1)/2,.1); % Upper/lower limits root bracket
gamma = 2*xi*B*(I_up-Ki*t);                         % Compute gamma
for m = j:n                                         % For loop
    if gamma(m) < 709.783                               % Check overflow condition
        y(1) = I_low; k = 1;                            % Root guess and iter.
        while (abs(r(y(k),t(m)))>rtol) && (k<kmax)          % While loop
            y(k+1) = y(k) - r(y(k),t(m))/dr(y(k),t(m));         % Next iterate
            k = k + 1;                                          % Increment iteration
        end                                                 % End while loop
        [I(m),I_low] = deal(y(k));                          % I(m) is equal to root
        A = exp(2*B*xi*(I(m)-Ki*t(m)));                     % Compute alpha
        i(m) = (((1-B)*Ks+B*Ki)*(A+B-1)-A*B*Ki)/(A+B-A*B-1);% Compute inf. rate
    elseif (m > 1)                                      % Patch for overflow
        if abs(i(m-1)-Ks) < itol                            % Constant inf rate?
            i(m) = Ks; I(m) = I(m-1) + dt(m)*i(m);              % Constant rate
        else                                                % Else
            flag = 2; break;                                    % flag = 2; break loop
        end                                                 % End if statement    
    else                                                % Else
        flag = 2; break;                                    % flag = 2; break loop
    end                                                 % End if statement
end                                                % End for loop

end
