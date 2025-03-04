function [t,dtdI,flag] = Haverkamp_t(eta,plugin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% This function solves for the t(I) relationship of Haverkamp                        %%
%%  SYNOPSIS: t = Haverkamp_t(eta,plugin)                                             %%
%% where                                                                              %%
%%  eta  [input] 3x1 vector with S [cm/h^0.5], Ks [cm/h], and beta [-]                %%
%%  plugin [input] structure with Haverkamp input variables                           %%
%%  t    [outpt] time, t, in hours corresponding to I (cm)                            %%
%%  dtdI [outpt] derivative of time, t, with respect to I, in cm                      %%
%%  flag [outpt] exit flag: [1] finite [2] infinite time values                       %%
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

%% Initialization part
S = eta(1); Ks = eta(2); Ki = 0;                    % Unpack parameter values
if numel(eta)<3, B = plugin.B; else, B = eta(3); end
dK = Ks-Ki; xi = dK/S^2; I = plugin.I(:);           % Compute dK and xi
n = numel(I); t = nan(n,1); flag = 1;               % Initialize time
%% Dynamic part: Evaluate time expression
gamma = 2*B*xi*I;                                   % Compute gamma
if I(1) == 0; t(1) = 0; j = 2; else, j = 1; end     % t(1) equals zero if I(1) is zero
for m = j:n                                         % For loop
    if gamma(m) < 709.783                               % Provision numerical overflow
        t(m) = 1./(dK*(B-1)*xi).*(1/2* ...
            log(exp(gamma(m))/B + (B-1)/B) - xi*I(m));      % Compute inf. time
    else
        flag = 2; break                                 % Break for loop; flag = 2
    end                                                 % End if statement
end                                                 % End for loop
dtdI = (exp((2*B*I*Ks)/S^2)) ./ ((B - 1) .* ... 
    (exp((2*B*I*Ks)/S^2)/B + (B - 1)/B)) - 1/(Ks*(B - 1));  

end

