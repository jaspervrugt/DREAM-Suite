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
%% Example 36: Inverse modeling of the parameters of a two class flocculation model   %%
%%             using temperature data and limits of acceptability                     %%
%%                                                                                    %%
%% ---------------------------------------------------------------------------------- %%

clc; clear; close all hidden;           % clear workspace and figures

%% Problem settings defined by user
DREAMPar.d = 2;                         % Dimension of the problem
DREAMPar.lik = 23;                      % Likelihood 23 (limits accptblty)

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'fold';        % Explicit boundary handling
Par_info.names = {'M','\tau_{\rm c}'};
Par_info.min = [ 0  1.8 ];              % Minimum parameter values
Par_info.max = [ 1  2.5 ];              % Minimum parameter values
% M:        Empirical erosion parameter [kg/m2/s]
% τ_c:      Critical shear stress for erosion [Pa]

%% Define name of function (.m file) for posterior exploration
Func_name = 'flocsedtran_1DV';

%% Load experimental data
data = readmatrix(['' ...               % Temperature data
    'MeasData_BLK2018.txt']);
[n,K] = size(data);                     % # rows and # columns of data
ts_obs = 3600 * data(1:n,1);            % observation times
data = data(1:n,3:K);                   % isolate right columns

%% Define training data
Meas_info.S = data(:);                  % Temperature data, summary metrics

%% Setup numerical simulation model
plugin = model_properties;              % plugin structure of model prpties
plugin.ts_obs = ts_obs;                 % Observation time

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                % # Markov chains
        DREAMPar.T = 2500;              % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                 % # Markov chains
        DREAMPar.T = 3000;              % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 200*ones(1,DREAMPar.d);    % # discrete steps
end

%% Optional settings
options.modout = 'yes';                 % Return model simulations?
options.epsilon = [ 2 ];                % Limits of acceptability
options.parallel = 'yes';               % Run each chain on different core
options.save = 'yes';                   % Save workspace DREAM during run
options.print = 'yes';                  % Figures printed to screen

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar, ...
    Par_info,Meas_info,options,[],plugin);

%% Save results
save results.mat chain output FX Z logL method Func_name DREAMPar ...
    Par_info Meas_info options plugin
