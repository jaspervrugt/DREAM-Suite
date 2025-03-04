function [t,flag] = Haverkamp_t_patch(eta,plugin,itol)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function solves for the t(I) relationship of Haverkamp                        %%
%%  SYNOPSIS: t = Haverkamp_t_patch(eta,I)                                            %%
%%  SYNOPSIS: t = Haverkamp_t_patch(eta,I,itol)                                       %%
%% where                                                                              %%
%%  eta  [input]       3x1 vector with S [cm/h^0.5], Ks [cm/h], and beta [-]          %%
%%  I    [input]       nx1 vector with cumulative infiltration, I, in cm              %%
%%  itol [opt. input]  tolerance on constant rate assumption      (default: 1e-10)    %%
%%  t    [output]      time, t, in hours corresponding to I (cm)                      %%
%%  flag [output]      exit flag: [1] exact [2] approximate                           %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
if nargin < 3, itol = 1e-10; end
%% Initialization part
S = eta(1); Ks = eta(2); B = eta(3); Ki = 0;        % Unpack parameter values
dK = Ks-Ki; xi = dK/S^2; I = plugin.I(:);           % Compute dK, xi, I is column vector
n = numel(I); t = nan(n,1); flag = 1;               % Initialize time
%% Dynamic part: Evaluate time expression
gamma = 2*B*xi*I;                                   % Compute gamma
if I(1) == 0; t(1) = 0; j = 2; else, j = 1; end     % t(1) equals zero if I(1) is zero
for m = j:n                                         % For loop
    if gamma(m) < 709.783                               % Numerical overflow?
        t(m) = 1./(dK*(B-1)*xi).*(1/2* ...                  % No: compute inf. time
            log(exp(gamma(m))/B + (B-1)/B) - xi*I(m));      
        A = exp(2*B*xi*(I(m)-Ki*t(m)));                     % Compute alpha
        i = (1-B)*Ks*(A+B-1)/(A+B-A*B-1);                   % Compute inf. rate        
    elseif (m > 1)                                      % Patch for overflow                    
        if abs(i-Ks) < itol                                 % Constant inf rate?
            t(m) = t(m-1) + (I(m)-I(m-1))/Ks; flag = 2;         % Constant rate
        else                                                % Else
            flag = 2; break                                     % flag = 2; break loop
        end                                                 % End if statement    
    else                                                % Else
        flag = 2; break                                     % flag = 2; break loop
    end                                                 % End if statement
end                                                 % End for loop