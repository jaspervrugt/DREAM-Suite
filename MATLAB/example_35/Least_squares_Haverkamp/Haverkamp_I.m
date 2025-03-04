function [I,i,flag] = Haverkamp_I(eta,plugin,rtol,kmax)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function solves for the I(t) relationship of Haverkamp using Newton's method  %%
%%  SYNOPSIS: I = Haverkamp_I(eta,plugin)                                             %%
%%            I = Haverkamp_I(eta,plugin,rtol)                                        %%
%%            I = Haverkamp_I(eta,plugin,rtol,kmax)                                   %%
%% where                                                                              %%
%%  eta  [input]      4x1 vector with S [cm/h^.5], Ks [cm/h], ß [-], Ki [cm/h]        %%
%%  plugin [inpt]     structure with Haverkamp input variables                        %%
%%  rtol [input]      OPT: tolerance on function value at root      (default: 1e-12)  %%
%%  kmax [input]      OPT: maximum number of Newton iterations      (default: 20)     %%
%%  I    [outpt]      cumulative infiltration, I (cm), as function of time, t (h)     %%
%%  i    [outpt]      infiltration rate, i (cm/h), as function of time, t (h)         %%
%%  flag [outpt]      exit flag: [1] exact [2] approximate                            %%
%%                                                                                    %%
%%  LITERATURE                                                                        %%
%%  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and               %%
%%      H. Vereecken (2023), The time validity of Philip's two-term infiltration      %%
%%      equation: An elusive theoretical quantity? Vadose Zone Journal, e20309,       %%
%%      pp. 1-25, https://doi.org/10.1002/vzj2.20309                                  %%
%%  Vrugt, J.A. and y. Gao (2022), On the three-parameter infiltration equation of    %%
%%      Parlange et al. (1982): Numerical solution, experimental design, and          %%
%%      parameter estimation, Vadose Zone Journal, 21:e20167, pp. 1-25,               %%
%%      https://doi.org/10.1002/vzj2.20167                                            %%
%%                                                                                    %%
%%  © Written by Jasper A. Vrugt, Jan. 2019                                           %%
%%  University of California Irvine                                                   %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if nargin < 4, kmax = 20; end 
if nargin < 3, rtol = 1e-12; end
if numel(eta) < 3, eta(3) = plugin.B; end
if numel(eta) < 4, eta(4) = 0; end
%% Initialization
S = eta(1); Ks = eta(2); B = eta(3); Ki = eta(4);   % Unpack parameter values
dK = Ks-Ki; xi = dK/S^2; t = plugin.t(:);           % Compute dK and xi, # elements t
n = numel(t); I = nan(n,1); flag = 1;               % Initial cum. inf.
if dK <= 0, [I,i] = deal(inf(n,1)); flag = 0; return, end 
if t(1) == 0, I(1) = 0; j = 2; else, j = 1; end
r = @(I,t) I-1/(2*xi)*log(exp(2*xi*B*(I-Ki*t))/B + ...
    (B-1)/B) - (dK*(1-B)+Ki)*t;                     % Residual function
dr = @(I,t) (1-(B*exp(2*xi*B*(I-Ki*t)))./ ...
    (exp(2*xi*B*(I-Ki*t))+B-1));                    % Derivative function
%% Dynamic part
I_up = S*sqrt(t) + t*Ks; I_low = max(I_up(1)/2,.1); % Upper/lower limits root bracket
gamma = 2*xi*B*(I_up-Ki*t);                         % Compute gamma 
for m = j:n                                         % For loop
     if gamma(m) < 709.783                              % Check overflow condition   
        y(1) = I_low; k = 1;                                % Yes: root guess and iter.    
        while (abs(r(y(k),t(m)))>rtol) && (k<kmax)          % While loop
            y(k+1) = y(k) - r(y(k),t(m))/dr(y(k),t(m));         %#ok Next iterate
            k = k + 1;                                          % Increment iteration 
        end                                                 % End while loop    
        [I(m),I_low] = deal(y(k));                          % I(m) is equal to root
    else                                                % Otherwise 
        flag = 2; break                                     % No: flag=2, break loop 
    end                                                 % End if statement
end                                                 % End for loop
A = exp(2*B*xi*(I-Ki*t));                           % Compute alpha
i = (((1-B)*Ks+B*Ki)*(A+B-1)-A*B*Ki)./(A+B-A*B-1);  % Compute infiltration rate

end
