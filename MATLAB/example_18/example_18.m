%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                    %%
%% DDDDD   RRRR    EEEEE    AA    MM   MM      SSSSSS  UU   UU   II   TTTTTTTT  EEEEE %%
%% DDDDDD  RRRR    EEEEE   AAAA   MM   MM      SSSSS   UU   UU   II   TTTTTTTT  EEEEE %%
%% DD  DD  RR RR   EE     AA  AA  MMM MMM      SS      UU   UU   II      TT     EE    %%
%% DD  DD  RR RR   EEE    AA  AA  MMMMMMM ---- SS      UU   UU   II      TT     EEE   %%
%% DD  DD  RRRRR   EEE    AAAAAA  MMM MMM ---- SSSSSS  UU   UU   II      TT     EEE   %%
%% DD  DD  RR RR   EE     AAAAAA  MM   MM          SS  UU   UU   II      TT     EE    %%
%% DDDDDD  RR  RR  EEEEE  AA  AA  MM   MM       SSSSS  UUUUUUU   II      TT     EEEEE %%
%% DDDDD   RR  RR  EEEEE  AA  AA  MM   MM      SSSSSS  UUUUUUU   II      TT     EEEEE %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                    %%
%% Example 18: Lotka-Volterra model: Informal likelihood (GLUE)                       %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. and K.J. Beven (2018), Embracing equifinality with efficiency:       %%
%%       Limits of Acceptability sampling using the DREAM_{(LOA)} algorithm,  Journal %%
%%       of Hydrology,  559 , pp. 954-971, doi:10.1016/j.jhydrol.2018.02.026          %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.d = 4;         % Dimension of the problem
DREAMPar.lik = 32;      % Informal likelihood function
DREAMPar.GLUE = 10;     % Value (default) of likelihood shape parameter

%% Initial sampling and parameter ranges
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling: reflection
Par_info.names = {'\alpha','\beta','\gamma','\delta'};
Par_info.min = [   0     0    0     0    ]; % Lower bound parameters
Par_info.max = [   1     10   1     10   ]; % Upper bound parameters
%% Load calibration data vector
Meas_info.Y = load('abundances.txt');   % Load food web dataset

%% Define name of function (.m file) for posterior exploration
Func_name = 'lotka_volterra';

%% Optional settings
options.parallel = 'no';   % Run chains in parallel
options.modout = 'yes';     % Store results

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_kzs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 2500;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 8000;                          % # generations
    case {'dream_kzs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 1000;                          % # generations
        DREAMPar.a_1 = 0;
        Meas_info.Sigma = 1;
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar, ...
    Par_info,Meas_info,options);