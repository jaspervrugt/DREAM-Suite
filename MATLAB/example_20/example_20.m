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
%% Example 20: Analytic solution of heat flow equation: Limits of Acceptability       %%
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
DREAMPar.lik = 23;      % Likelihood 23 (limits of acceptability)

%% Initial sampling and parameter ranges
Par_info.initial = 'latin';         % Latin hypercube sampling
Par_info.boundhandling = 'reflect'; % Explicit boundary handling: reflection
Par_info.names = {'T_{\rm a}','A_{0}','\phi','d'};
Par_info.min = [ 10  0   -12    0  ]; % Lower bound parameters
Par_info.max = [ 30  20   12   500 ]; % Upper bound parameters

% T_a = Annual average temperature    [°C]
% A_0 = Temperature amplitude         [°C]
% fi  = Phase constant                [rad]
% d   = Characteristic damping depth  [cm]

%% Define name of function (.m file) for posterior exploration
Func_name = 'heatflow';

%% Load calibration data vector (summary metrics with which simulations are compared)
Meas_info.S = load('temp_data.txt');  % Temperature data are summary statistics

%% Optional settings (limits of acceptability)
options.epsilon = 2;        % degrees Celsius
options.modout = 'yes';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
% method = 'mtdream_zs';
method = 'dream_dzs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 2500;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 5000;                          % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 200*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options);
