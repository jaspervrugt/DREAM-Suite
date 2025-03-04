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
%% Example 30: Predator prey interactions                                             %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% DOWNLOAD AT:
%% http://izt.ciens.ucv.ve/ecologia/Archivos/ECO_POB%202007/ECOPO2_2007/Cariboni%20et%20al%202007.pdf

%% Problem settings defined by user
DREAMPar.d = 4;                             % Dimension of the problem
DREAMPar.lik = 11;                          % Gaussian likelihood function

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling
Par_info.names = {'r','\alpha','m','\theta'};
Par_info.min = [ 0.8  0.2  0.6  0.05 ]; % minimum parameter values
Par_info.max = [ 1.8  1.0  1.0  0.15 ]; % maximum parameter values

%% Define name of function (.m file) for posterior exploration
Func_name = 'pred_prey';

%% Define model setup
plugin.u0 = [ 8 2 ];    % initial count of each of the two species
plugin.ode_options = odeset('abstol',1e-5,...   % settings ode solver
    'reltol',1e-5,'maxstep',0.1,...             % state trajectories subject to chaos)
    'initialstep',0.1);
plugin.t = 0:1:60;      % specify the print time step and the end time

%% Load measured species abundances of predator and prey species (= 1 column)
%% Measured data --> x = [ 1 0.5 0.8 0.1 ]; data = pred_prey(x,plugin);
%% save data.txt -ascii data
load data.txt;
Y(1:61,1) = data(1:61,1) + normrnd(0,0.1*std(data(1:61,1)),61,1);
Y(62:122,1) = data(62:122,1) + normrnd(0,0.1*std(data(62:122,1)),61,1);
%% Now save data
Meas_info.Y = Y;

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                        % # Markov chains
        DREAMPar.T = 2500;                      % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                         % # Markov chains
        DREAMPar.T = 5000;                      % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 250*ones(1,DREAMPar.d);    % # discrete steps
end

%% Optional settings
options.modout = 'yes';         % Return model simulations of samples?
options.parallel = 'no';        % Run each chain on a different core
options.save = 'yes';           % Save workspace DREAM during run

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);

% Postprocess: Create one matrix of Markov chains
pars = genparset(chain);
% Now apply burn-in to parameter values
M = size(pars,1); pars = pars(floor(.8*M):M,1:DREAMPar.d);
% Apply burn-in to simulated abundances
Y = FX(floor(.8*M):M,:);
