function [I,i,flag,K] = Haverkamp_IK(eta,plugin,rtol,kmax,method)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function solves for the I(t) relationship of Haverkamp using Newton's method  %%
%%  SYNOPSIS: I = Haverkamp_I(eta,t)                                                  %%
%%            I = Haverkamp_I(eta,t,rtol)                                             %%
%%            I = Haverkamp_I(eta,t,rtol,maxiter)                                     %%
%% where                                                                              %%
%%  eta  [input]       4x1 vector with S [cm/h^0.5], Ks [cm/h], beta [-], Ki [cm/h]   %%
%%  t    [input]       nx1 vector with time, t, in hours (h)                          %%
%%  rtol [opt. input]  tolerance on function value at root        (default: 1e-12)    %%
%%  kmax [opt. input]  maximum number of Newton iterations        (default: 20)       %%
%%  I    [output]      cumulative infiltration, I (cm), as function of time, t (h)    %%
%%  i    [output]      infiltration rate, i (cm/h), as function of time, t (h)        %%
%%  flag [output]      exit flag: [1] exact [2] approximate                           %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if nargin < 4, kmax = 20; end 
if nargin < 3, rtol = 1e-12; end
if numel(eta) == 3, eta(4) = 0; end
%% Initialization
S = eta(1); Ks = eta(2); B = eta(3); Ki = eta(4);   % Unpack parameter values
dK = Ks-Ki; xi = dK/S^2; t = plugin.t(:);           % Compute dK and xi, # elements t
n = numel(t); I = nan(n,1); I_low = 0.05; flag = 1;  % Initial cum. inf.
K = nan(n,1);
if t(1) == 0, I(1) = 0; K(1) = 0; j = 2; else, j = 1; end
r = @(I,t) I-1/(2*xi)*log(exp(2*xi*B*(I-Ki*t))/B + ...
    (B-1)/B) - (dK*(1-B)+Ki)*t;                     % Residual function
dr = @(I,t) (1-(B*exp(2*xi*B*(I-Ki*t)))./ ...
    (exp(2*xi*B*(I-Ki*t))+B-1));                    % Derivative function
%% Dynamic part
I_up = S*sqrt(t) + t*Ks;                            % Upper limit root bracket
gamma = 2*xi*B*(I_up-Ki*t);                         % Compute gamma 
for m = j:n                                         % For loop
     if gamma(m) < 709.783                              % Check overflow condition   
        Ii = [ I_low , 1/2*(I_low+I_up(m)) , I_up(m) ];
        y(1) = Ii(method); k = 1;                           % Yes: root guess and iter.    
        while (abs(r(y(k),t(m)))>rtol) && (k<kmax)          % While loop
            y(k+1) = y(k) - r(y(k),t(m))/dr(y(k),t(m));         % Next iterate
            k = k + 1;                                          % Increment iteration 
        end                                                 % End while loop    
        [I(m),I_low] = deal(y(k));                          % I(m) is equal to root
    else                                                % Otherwise 
        flag = 2; break                                     % No: flag=2, break loop 
    end                                                 % End if statement
    K(m) = k;
end                                                 % End for loop
A = exp(2*B*xi*(I-Ki*t));                           % Compute alpha
i = (((1-B)*Ks+B*Ki)*(A+B-1)-A*B*Ki)./(A+B-A*B-1);  % Compute infiltration rate