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
%% Example 24: Soil water retention functions as parametric expressions for the flow  %%
%%             duration curve. A separate toolbox, called FDCFIT, has been developed  %%
%%             for FDC-Fitting. This includes 15 different FDC functions and many     %%
%%             other functionalities                                                  %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Sadegh, M., J.A. Vrugt, X. Cu, and H.V. Gupta (2016), The soil water             %%
%%       characteristic as new class of parametric expressions of the flow duration   %%
%%       curve, Journal of Hydrology, 535, 438-456, doi:10.1016/j.jhydrol.2016.01.027 %%
%%   Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model calibration and       %%
%%       evaluation: Approximate Bayesian computation, Water Resources Research, 49,  %%
%%       4335–4345, doi:10.1002/wrcr.20354                                            %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

DREAMPar.lik = 13;        % Log-likelihood with AR(1) of residuals

%% Initial sampling and parameter ranges
Par_info.initial = 'latin';         % Latin hypercube sampling
Par_info.boundhandling = 'reflect'; % Explicit boundary handling : reflection

% Par. names :    a    b    c
par_names = {'a','b','c'};
fpar_mod   = [ 10   5.0  1 ];       % Default values of a, b and c
parmin_mod = [ 1e-6 1.0  1e-6 ];    % Lower bound coefficients of VG-3
parmax_mod = [ 100  10   10   ];    % Upper bound coefficients of VG-3

%% Define name of function (.m file ) for posterior exploration
Func_name = 'FDC_vg';

global LV       % Global for likelihood functions 13,14,16,17,44,45
% Be careful with its use
switch DREAMPar.lik
    case 13 %% Normal distribution
        % index:           1   2    3   4
        % parname:        s0  s1  phi1 phi2
        fpar_nuis   =  [  0.1  0    0   0 ];
        parmin_nuis =  [   0   0    0   0 ];
        parmax_nuis =  [   1   1    1   1 ];
        LV.filename = 'Normal'; id_nuis = 3;
    case 16 %% Laplace distribution
        % index:           1   2    3
        % parname:        s0  s1  phi1
        fpar_nuis   =  [  0.1  0    0 ];
        parmin_nuis =  [   0   0    0 ];
        parmax_nuis =  [   1   1    1 ];
        LV.filename = 'Laplace'; id_nuis = 3;
    case 17 %% SST
        % index:          1   2    3   4    5    6
        % parname:       s0  s1  nu   xi phi1  phi2
        fpar_nuis   =  [ 0.1  0  1e10  1    0    0 ];
        parmin_nuis =  [  0   0    2   0.1  0    0 ];
        parmax_nuis =  [  1   1   100  10   1    1 ];
        LV.filename = 'SL'; id_nuis = [3 4 5];
    case 44 %% GL+
        % index:          1   2    3   4    5    6
        % parname:       s0  s1  beta xi phi1  phi2
        fpar_nuis   =  [ 0.1  0    0   1    0    0 ];
        parmin_nuis =  [  0   0   -1  0.1   0    0 ];
        parmax_nuis =  [  1   1    1   10   1    1 ];
        LV.filename = 'GL_plus'; id_nuis = [3 4 5];
    case 45 %% SGT
        % index:          1   2    3   4   5    6    7
        % parname:       s0  s1  labda p   q  phi1 phi2
        fpar_nuis   =  [ 0.1  0    0   2  1e10  0    0 ];
        parmin_nuis =  [  0   0   -1  0.5   2   0    0 ];
        parmax_nuis =  [  1   1    1  100  100  1    1 ];
        LV.filename = 'UL'; id_nuis = [1 3 4 5 7];
end

% Load record of daily streamflow data ( French Broad )
data = load('03451500.dly');
% Compute the daily flow duration curve, FDC - discharge in column 6
[e,y,p_0] = efdc(data,'daily',6);
Meas_info.Y = y;                % Measured discharge data
Meas_info.sigma2 = 'constant';  % Constant measurement error variance
plugin.E = e;                   % Exceedance probabilities of measured discharge record

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

%% Optional settings
options.parallel = 'no';    % Model is so fast -> parallel slows down
options.IO = 'no';          % No file writing used with FDC_VG model
options.save = 'yes';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 10000;                         % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 15000;                         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);
