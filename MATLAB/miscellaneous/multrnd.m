function [X,Y] = multrnd(n,p,m)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% MULTRND Multinomial random sequence of m simulations of k outcomes    %%
%%    with p probabiltites in n trials                                   %%
%%    Y = MULTRND(N,P,M) generates a random sequence of m simulations    %% 
%%    of k integers from multinomial distribution with n trials and k    %%
%%    outcomes, where the probability for each simulation (draw) equals  %%
%%                     n!                                                %%
%%           ---------------------  ×  p1^n1 × p2^n2 ×  . . .  × pk^nk   %%
%%           n1! × n2! × ... × nk!                                       %%
%%                                                                       %%
%%    Then, a single sample (n1, n2,..., nk) have a multinomial joint    %%
%%    distribution with parameters n and p1, p2,..., pk. Parameter n is  %%
%%    called the number of trials; parameters p1, p2,..., pk are called  %%
%%    the category probabilities; k is called the number of categories   %%
%%                                                                       %%
%% SYNOPSIS: [X,Y] = multrnd(n,p,m)                                      %%
%%  where                                                                %%
%%   n      [input] # of trials                                          %%
%%   p      [input] probabilities of categories (with sum(p) = 1)        %%
%%   m      [input] # of simulations            (default = 1)            %%
%%   X      [outpt] multinomial random deviates (= default)              %%
%%   Y      [outpt] multinomial probabilities sampled random deviates    %%
%%                                                                       %%
%% Example 1: We are interested to get a multinomial random values of 4  %%
%%            outcomes with probabilities p = [0.1 0.06 .7 .14] with     %%
%%            n = 2678 trials in m = 10 simulations                      %%
%%                                                                       %%
%% Calling on Matlab the function:                                       %%
%%            [X Y] = multrnd(n,p,m)                                     %%
%% We yield                                                              %%
%%   X = [ 271       152       1873      382                             %%
%%         249       154       1890      385                             %%
%%         266       172       1862      378                             %%
%%         290       147       1882      359                             %%
%%         247       153       1873      405                             %%
%%         291       155       1842      390                             %%
%%         268       141       1900      369                             %%
%%         248       158       1899      373                             %%
%%         267       181       1855      375                             %%
%%         259       175       1884      360 ]                           %%
%%                                                                       %%
%%   Y = [ 0.1012    0.0568    0.6994    0.1426                          %%
%%         0.0930    0.0575    0.7058    0.1438                          %%
%%         0.0993    0.0642    0.6953    0.1412                          %%
%%         0.1083    0.0549    0.7028    0.1341                          %%
%%         0.0922    0.0571    0.6994    0.1512                          %%
%%         0.1087    0.0579    0.6878    0.1456                          %%
%%         0.1001    0.0527    0.7095    0.1378                          %%
%%         0.0926    0.0590    0.7091    0.1393                          %%
%%         0.0997    0.0676    0.6927    0.1400                          %%
%%         0.0967    0.0653    0.7035    0.1344 ]                        %%
%%                                                                       %%
%% Written by A. Trujillo-Ortiz, R. Hernandez-Walls and A. Castro-Perez  %%
%%    Facultad de Ciencias Marinas                                       %%
%%    Universidad Autonoma de Baja California                            %%
%%    Apdo. Postal 453                                                   %%
%%    Ensenada, Baja California                                          %%
%%    Mexico                                                             %%
%%    Email: atrujo@uabc.mx                                              %%
%%    Copyright © January 12, 2005                                       %%
%%                                                                       %%
%% To cite this file, this would be an appropriate reference:            %%
%%    Trujillo-Ortiz, A., R. Hernandez-Walls and A. Castro-Perez (2005)  %%
%%    multrnd: Multinomial random sequence. A MATLAB file. WWW document. %%
%%    URL = http://www.mathworks.com/matlabcentral/fileexchange/...      %%
%%    loadFile.do?objectId=6788                                          %%
%%                                                                       %%
%% References:                                                           %%
%%    Abramowitz, M. and I.A. Stegun (1964), Handbook of Mathematical    %%
%%    Functions, Government Printing Office, 26.1.20. Available on       %%
%%    Internet at the address http://hcohl.shell42.com/as/frameindex.htm %%
%%                                                                       %%
%% Modified by Jasper A. Vrugt, Dec. 2016                                %%
%% University of California Irvine                                       %%
%% DREAM Package                                                         %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 3
    if nargin == 2
        m = 1;
        warning(['DREAM_PACKAGE:multrnd:Two input arguments only, ',...
            'and, thus, the value of m = 1 [# simulations] (default)']);
    else
        error(['DREAM_PACKAGE:multrnd:You must submit at least two ',...
            'input arguments, n and p.']);
    end
end

if (numel(n) ~= 1) || (fix(n) ~= n) || (n < 0)
   error(['DREAM_PACKAGE:multrnd:Scalar n [# trials] must be a ',...
       'positive integer.']);
end

np = numel(p); P = sum(p);
if P ~= 1
    error(['DREAM_PACKAGE:multrnd:Sum of category probabilities, ',...
       'p, must equal 1.']);
end

X = nan(m,np);      % Initialize return argument
s = cumsum(p/P);    % CDF of category probabilities    
for i = 1:m
    o = ones(1,n);
    r = rand(1,n);
    for j = 1:np
        o = o + (r > s(j)); 
    end
    for j = 1:np
        X(i,j) = sum(o == j); 
    end
end

Y = X./n;           % Compute category probabilities according to draws

end
