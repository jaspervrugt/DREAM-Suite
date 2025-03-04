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
%% Example 25: Bottom-up control hypothesis: Bedrock depth from high-resolution       %%
%%             topographic data and a geomorphologic model of landscape evolution     %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Gomes, G.J.C., J.A. Vrugt, and E.A. Vargas Jr., Towards improved prediction of   %%
%%       the bedrock depth underneath hillslopes: Bayesian Inference of the bottom-up %%
%%       control hypothesis using high-resolution topographic data, Water Resources   %%
%%       Research, 52, 3085–3112, doi:10.1002/2015WR018147, 2016                      %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

DREAMPar.lik = 13;                          % Gaussian log-likelihood function

% Provide information parameter space and initial sampling
Par_info.initial = 'latin';                 % Latin hypercube sampling
Par_info.boundhandling = 'reflect';         % Explicit boundary handling
par_names = {'\Phi','\lambda_{1}',...
    'S_{\rm c}','\lambda_{2}'};
fpar_mod   = [  0.001   1.0  1.0  10 ];     % Default values parameters
parmin_mod = [  0.0001  0.1  0.8  1  ];     % Lower bound parameters
parmax_mod = [  0.1     3.0  1.5  20 ];     % Upper bound parameters
% Define name of function (.m file) for posterior exploration
Func_name = 'DTB_model';

% Load depth to bedrock calibration data
load DTB_data.txt

global LV       % Global for likelihood functions 13,14,16,17,44,45
% Be careful with its use
switch DREAMPar.lik
    case 13 %% Normal distribution
        % index:           1   2    3   4
        % parname:        s0  s1  phi1 phi2
        fpar_nuis   =  [  0.1  0    0   0 ];
        parmin_nuis =  [   0   0    0   0 ];
        parmax_nuis =  [   5   1    1   1 ];
        LV.filename = 'Normal'; id_nuis = [ 1 2 ];
    case 16 %% Laplace distribution
        % index:           1   2    3
        % parname:        s0  s1  phi1
        fpar_nuis   =  [  0.1  0    0 ];
        parmin_nuis =  [   0   0    0 ];
        parmax_nuis =  [   5   1    1 ];
        LV.filename = 'Laplace'; id_nuis = [2 3];
    case 17 %% SST
        % index:          1   2    3   4    5    6
        % parname:       s0  s1  nu   xi phi1  phi2
        fpar_nuis   =  [ 0.1  0  1e10  1    0    0 ];
        parmin_nuis =  [  0   0    2   0.1  0    0 ];
        parmax_nuis =  [  5   1   100  10   1    1 ];
        LV.filename = 'SL'; id_nuis = [2 3 4];
    case 44 %% GL+
        % index:          1   2    3   4    5    6
        % parname:       s0  s1  beta xi phi1  phi2
        fpar_nuis   =  [ 0.1  0    0   1    0    0 ];
        parmin_nuis =  [  0   0   -1  0.1   0    0 ];
        parmax_nuis =  [  5   1    1   10   1    1 ];
        LV.filename = 'GL_plus'; id_nuis = [2 3 4];
    case 45 %% SGT
        % index:          1   2    3   4   5    6    7
        % parname:       s0  s1  labda p   q  phi1 phi2
        fpar_nuis   =  [ 0.1  0    0   2  1e10  0    0 ];
        parmin_nuis =  [  0   0   -1  0.5   2   0    0 ];
        parmax_nuis =  [  5   1    1  100  100  1    1 ];
        LV.filename = 'UL'; id_nuis = [2 3 4 5];
end

Meas_info.Y = DTB_data(:,1);    % Measured discharge data
Meas_info.sigma2 = 'constant';  % Constant measurement error variance
%Meas_info.Sigma = 1.9;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
nuis_names = nuis_var_names(DREAMPar);              % Extract names nuis. variables
names = [par_names,nuis_names];                     % Combine with parameter names
n_modpar = numel(fpar_mod);                         % # model parameters
LV.id_vpar = [1:n_modpar , id_nuis + n_modpar ];    % Model parameter + nuisance variable selection
LV.fpar = [fpar_mod fpar_nuis];                     % Merge default values parameters and nuisance variables
parmin = [parmin_mod parmin_nuis];
parmax = [parmax_mod parmax_nuis];                  % Merge the parameter names and nuisance variables
Par_info.min = parmin(LV.id_vpar);
Par_info.max = parmax(LV.id_vpar);                  % Min/max values of parameter selection
Par_info.names = names(LV.id_vpar);                 % Names of parameters and nuisance variables
DREAMPar.d = size(Par_info.min,2);                  % # parameters
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Optional settings
options.modout = 'yes';                 % Store the model output
options.save = 'yes';                   % Save DREAM output during the run

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_d';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                % # Markov chains
        DREAMPar.T = 10000;             % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                 % # Markov chains
        DREAMPar.T = 15000;             % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 200*ones(1,DREAMPar.d);    % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options);
