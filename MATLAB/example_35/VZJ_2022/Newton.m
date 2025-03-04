function [Ir,exitflag] = Newton(r,dr,I0,maxiter,tolfun)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Newton's method for root finding of infiltration equation of Haverkamp             %%
%%  SYNOPSIS: Ir = newton(fun,I0);                                                    %%
%%            [Ir,exitflag] = Newton(r,dr,I0);                                        %%
%%            [Ir,exitflag] = Newton(r,dr,I0,maxiter);                                %%
%%            [Ir,exitflag] = Newton(r,dr,I0,maxiter,tolfun);                         %%
%% where                                                                              %%
%%       r:        [input] anonymous function handle residual function                %%
%%       dr:       [input] anonymous function handle derivative of residual function  %%
%%       I0:       [input] initial value of cum. inf. in cm                           %%
%%       maxiter:  [optional input] maximum number of iterations (default: 20)        %%
%%       tolfun:   [optional input] tolerance on root function value (default: 1e-10) %%
%%       Ir:       [output] cumulative infiltration (in cm) at specified time         %%
%%       exitflag: [output] exit condition. Possible conditions are:                  %%
%%                          1: Newton found a zero point ( = root)                    %%
%%                         -2: Newton terminated with infinite root                   %%
%%                         -3: Newton terminated with infinite function value root    %%
%%                         -4: Newton terminated with imaginary root                  %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization part: Define settings of Secant method
if nargin < 5, tolfun = 1e-12; end
if nargin < 4, maxiter = 20; end
I(1) = I0;                                              % Set initial solution

%% Dynamic part: Iteratively refine estimate of root
for k = 2:maxiter
    dI = r(I(k-1))/dr(I(k-1));                              % Increment of I
    I(k) = I(k-1) - dI;                                     % Update root guess
    if abs(dI) < tolfun                                     % Check if converged
        break
    end
end
%% End of Dynamic Part

Ir = I(k); fIr = r(Ir);                                 % Final entry is root
if isfinite(Ir)                                         % Determine exitflag
    exitflag = 1;
elseif ~isfinite(Ir)
    exitflag = -2;
elseif ~isfinite(fIr)       
    exitflag = -3; 
elseif ~isreal(fIr)
    exitflag = -4;
end