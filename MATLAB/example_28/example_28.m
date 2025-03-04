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
%% Example 28: Conceptual model parameter and rainfall estimation                     %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten (2005),     %%
%%       Improved treatment of uncertainty in hydrologic modeling: Combining the      %%
%%       strengths of global optimization and data assimilation, Water Resources      %%
%%       Research, 41, W01017, doi:10.1029/2004WR003059                               %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.nCR = 10;      % Increase subspace sampling (default is 3)
DREAMPar.delta = 2;     % Delta is 2
DREAMPar.lik = 16;      % Model output is simulation: Laplacian likelihood

% DREAMPar.lik = 13;    % Gaussian likelihood using discharge data only
% DREAMPar.lik = 16;    % Laplace likelihood using discharge data only
% DREAMPar.lik = XX;    % discharge and rainfall data in likelihood (not defined herein)

%% Parameter names of hmodel
par_names = {'I_{\rm max}','S_{\rm u,max}','Q_{\rm s,max}',...
    '\alpha_{\rm e}','\alpha_{\rm f}','K_{\rm f}','K_{\rm s}'};
% index:         1     2     3      4    5     6     7
% parname:      Imax  Smax  Qsmax  alE  alF  Kfast Kslow
fpar_mod   =  [  1    100   10     100   0     2    70  ];
parmin_mod =  [ 0.5   10     0     1e-6  -10   0     0  ];    % If 'latin', min values
parmax_mod =  [ 10    1000  100    100   10   10    150 ];    % If 'latin', max values

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';         % Latin hypercube sampling
Par_info.boundhandling = 'reflect'; % Explicit boundary handling
%Par_info.norm = 1;                 % Work in normalized space

global LV                           %#ok Global for likelihood functions 13,14,16,17,44,45
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
        LV.filename = 'Laplace'; id_nuis = [ 1 2 ];
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
Func_name = 'rainfall_runoff';

%% Define the experimental data to be used for inference
data = load_data_dly('09497500');       % salt river data
id = [1 2 3 6 5 4]; data = data(:,id);  % Reorder data, 1: year
%                                                       2: month
%                                                       3: day
%                                                       4: discharge [mm/d]
%                                                       5: potential ET [mm/d],
%                                                       6: rainfall [mm/d]
plugin.Tmax = 2 * 365 + 65;             % max time --> 2 years with spin-up of 65 days (= illustration only)
Meas_info.Y = data(65:plugin.Tmax,4);   % Define measured discharge

%% Detect rainfall storms and add multipliers to hmodel parameter vector
plugin.nmod = numel(parmin_mod);        % # hmodel parameters
plugin.P = data(1:plugin.Tmax,6);       % Daily rainfall data
plugin.Y = data(1:plugin.Tmax,4);       % Daily discharge data
plugin.Ep = data(1:plugin.Tmax,5);      % Daily PET data
[plugin,mult,mult_names,fpar_mult] = ...
    check_rainfall(plugin);             % Isolate storm events
parmin_mod = [parmin_mod mult.min];     % Add multiplier minimum
parmax_mod = [parmax_mod mult.max];     % Add multiplier maximum
par_names = [par_names,mult_names];     % Add multiplier names
fpar_mod = [fpar_mod , fpar_mult];      % Add default multipliers

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuis_names = nuis_var_names(DREAMPar);          % Extract names nuis. variables
names = [par_names,nuis_names];                 % Combine with parameter names
n_modpar = numel(fpar_mod);                     % # model parameters
LV.id_vpar = [1:n_modpar , id_nuis + n_modpar]; % Model parameter + nuisance variable selection
LV.fpar = [fpar_mod fpar_nuis];                 % Merge default values parameters and nuisance variables
parmin = [parmin_mod parmin_nuis];              % Merge the parameter names and nuisance variables
parmax = [parmax_mod parmax_nuis];
Par_info.min = parmin(LV.id_vpar);              % Min/max values of parameter selection
Par_info.max = parmax(LV.id_vpar);
Par_info.names = names(LV.id_vpar);             % Names of parameters and nuisance variables
DREAMPar.d = size(Par_info.min,2);              % # parameters
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define hmodel settings for numerical simulation
plugin.y0 = 1e-5 * ones(5,1);                   % initial states (five reservoirs almost empty)
plugin.tout = 0:plugin.Tmax;                    % output simulation times of hmodel
plugin.data = struct('P',plugin.P,'Ep',plugin.Ep);
%% Integration settings for hmodel
plugin.options.InitialStep = 1;                 % initial time-step (d)
plugin.options.MaxStep     = 1;                 % maximum time-step (d)
plugin.options.MinStep     = 1e-6;              % minimum time-step (d)
plugin.options.RelTol      = 1e-6;              % relative tolerance
plugin.options.AbsTol      = 1e-6*ones(5,1);    % absolute tolerances (mm)
plugin.options.Order       = 2;                 % 2nd order accurate method (Heun)

%% We specify the discharge measurement error
Meas_info.Sigma = 0.1 * Meas_info.Y + 0.01;

%% Define optional settings
options.modout = 'yes';     % Return model (function) simulations of samples?
options.parallel = 'no';   % Evolve each chain on a different node?
options.save = 'yes';       % Save workspace of DREAM en route to posterior
options.restart = 'no';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 50;                            % # Markov chains
        DREAMPar.T = 20000;                         % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 4;                             % # Markov chains
        DREAMPar.T = 200000;                        % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);
