function [Ir,exitflag] = secant(fun,Ii0,Ii1,maxiter,tolfun)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The Secant method for root finding of infiltration equation of Haverkamp           %%
%%  SYNOPSIS: Ir = secant(fun,Ii);                                                    %%
%%            [Ir,exitflag] = secant(fun,Ii);                                         %%
%%            [Ir,exitflag] = secant(fun,Ii,maxiter);                                 %%
%%            [Ir,exitflag] = secant(fun,Ii,maxiter,tolfun);                          %%
%% where                                                                              %%
%%       fun:      [input] anonymous function handle of right-hand-side Haverkamp     %%
%%       Ii:       [input] starting point (cumulative infiltration in cm)             %%
%%       maxiter:  [optional input] maximum number of iterations (default: 20)        %%
%%       tolfun:   [optional input] tolerance on root function value (default: 1e-10) %%
%%       Ir:       [output] cumulative infiltration (in cm) at specified time         %%
%%       exitflag: [output] exit condition. Possible conditions are:                  %%
%%                          1: secant found a zero point ( = root)                    %%
%%                          2: secant terminated with infinite root                   %%
%%                          3: secant terminated with infinite function value root    %%
%%                          4: secant terminated with imaginary function value root   %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization part: Define settings of Secant method
if nargin < 5, tolfun = 1e-12; end
if nargin < 4, maxiter = 20; end
exitflag = 1; n = 3;                                          % Exitflag and iteration
I(n-2) = Ii0; I(n-1) = Ii1; % + max(0.1,0.25*Ii);                  % First pair of triple
I(n) = (I(1)*fun(I(2))-I(2)*fun(I(1)))/(fun(I(2))-fun(I(1))); % Last point of triple

%% Dynamic part: Iteratively refine estimate of root 
while (abs(fun((I(n)))) > tolfun) && ((n-2) < maxiter)        % Secant method
    n = n + 1;                                                % Increment iteration
    I(n) = (I(n-2)*fun(I(n-1))-I(n-1)*fun(I(n-2))) ...        % Update guess of root
        / (fun(I(n-1))-fun(I(n-2)));
end 
%% End of Dynamic Part

Ir = I(n);                                                    % Final entry is root
fIr = fun(Ir);                                                % fun(root)   
if ~isfinite(Ir)                                              % determine exitflag
    exitflag = 2;
elseif ~isfinite(fIr)       
    exitflag = 3; 
elseif ~isreal(fIr)
    exitflag = 4;
end