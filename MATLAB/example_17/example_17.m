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
%% Example 17: Conceptual watershed model with normal likelihood with AR(1) and       %%
%%             hetoroscedasticity                                                     %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Bikowski, J., J.A. Huisman, J.A. Vrugt, H. Vereecken, and J. van der             %%
%%       Kruk (2012), Inversion and sensitivity analysis of ground penetrating radar  %%
%%       data with waveguide dispersion using deterministic and Markov chain Monte    %%
%%       Carlo methods, Near Surface Geophysics, Special issue "Physics-based         %%
%%       integrated characterization", 10(6), 641-652,                                %%
%%       doi:10.3997/1873-0604.2012041, 2012                                          %%
%%   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson (2009),             %%
%%       Equifinality of formal (DREAM) and informal (GLUE) Bayesian approaches in    %%
%%       hydrologic modeling?, Stochastic Environmental Research and Risk Assessment, %%
%%       23(7), 1011-1026, doi:10.1007/s00477-008-0274-y                              %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.lik = 13; % Normal log-likelihood - AR(1)

%% Parameter names of hymod
par_names = {'S_{\rm u,max}','\beta','\alpha','K_{\rm s}','K_{\rm f}'};
% index:          1     2     3     4     5
% parname:      Sumax  beta  alfa   Ks    Kf
fpar_mod   =  [  200   1.00  0.5   0.01  0.50 ];
parmin_mod =  [   50   0.10  0.0   1e-5  0.25 ];
parmax_mod =  [ 1000   10.0  1.0   0.25  5.00 ];

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';                 % Latin hypercube sampling
Par_info.boundhandling = 'reflect';         % Explicit boundary handling
%Par_info.norm = 1;                          % Work in normalized space

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

%% Define name of function (.m file) for posterior exploration
Func_name = 'hymod';

%% Load hydrologic data
daily_data = load('03451500.dly');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Define warm-up period, first and last index of measured discharge record
plugin.data.wmp = 364; idx_t1 = 365; idx_t2 = size(daily_data,1);
idx_p = idx_t1 - plugin.data.wmp : idx_t2;
plugin.data.P = daily_data(idx_p,4); plugin.data.Ep = daily_data(idx_p,5);
plugin.tout = 0:1:size(plugin.data.Ep,1);
plugin.idx = 1 + ( plugin.data.wmp : numel(idx_p) )';
%% Initial conditions
plugin.y0 = 1e-5 * ones(6,1);
%% Integration options
plugin.options.InitialStep = 1;     % initial time-step (d)
plugin.options.MaxStep     = 1;     % maximum time-step (d)
plugin.options.MinStep     = 1e-3;  % minimum time-step (d)
plugin.options.RelTol      = 1e-3;  % relative tolerance
plugin.options.AbsTol      = 1e-3;  % absolute tolerances (mm)
plugin.options.Order       = 2;     % 2nd order accurate method (Heun)
%% Measured streamflow data in mm/d
Meas_info.Y = daily_data(idx_t1:idx_t2,6);
Meas_info.sigma2 = 'nonconstant';   % Non-constant measurement error variance
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

%% Optional setting
%options.parallel = 'yes';  % Run in parallel
options.modout = 'yes';     % Return model (function) simulations of samples (yes/no)?
options.save = 'yes';       % Save memory during run (for restart trials)

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;            % # Markov chains
        DREAMPar.T = 5000;          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;             % # Markov chains
        DREAMPar.T = 25000;         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);
