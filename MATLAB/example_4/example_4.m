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
%% Example 4: Distribution-adaptive likelihood functions applied to conceptual        %%
%%            watershed models using measured discharge data                          %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A., D.Y. de Oliveira, G. Schoups, and C.G.H. Diks (2022), On the use of %%
%%       distribution-adaptive likelihood functions: Generalized and universal        %%
%%       likelihood functions, scoring rules and multi-criteria ranking, Journal of   %%
%%       Hydrology, 615, Part B, 2022, doi:10.1016/j.jhydrol.2022.128542.             %%
%%       https://www.sciencedirect.com/science/article/pii/S002216942201112X          %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Schoups, G., J.A. Vrugt, F. Fenicia, and N.C. van de Giesen (2010), Corruption   %%
%%       of accuracy and efficiency of Markov Chain Monte Carlo simulation by         %%
%%       inaccurate numerical implementation of conceptual hydrologic models, Water   %%
%%       Resources Research, 46, W10530, doi:10.1029/2009WR008648                     %%
%%   Schoups, G., and J.A. Vrugt (2010), A formal likelihood function for parameter   %%
%%       and predictive inference of hydrologic models with correlated,               %%
%%       heteroscedastic and non-Gaussian errors, Water Resources Research, 46,       %%
%%       W10531, doi:10.1029/2009WR008933                                             %%
%%   Vrugt, J.A., C.J.F. ter Braak, H.V. Gupta, and B.A. Robinson (2009),             %%
%%       Equifinality of formal (DREAM) and informal (GLUE) Bayesian approaches in    %%
%%       hydrologic modeling?, Stochastic Environmental Research and Risk Assessment, %%
%%       23(7), 1011-1026, doi:10.1007/s00477-008-0274-y                              %%
%%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and        %%
%%       J.M. Hyman (2009), Accelerating Markov chain Monte Carlo simulation by       %%
%%       differential evolution with self-adaptive randomized subspace sampling,      %%
%%       International Journal of Nonlinear Sciences and Numerical Simulation, 10(3), %%
%%       271-288                                                                      %%
%%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), %%
%%       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       %%
%%       backward with Markov chain Monte Carlo simulation, Water Resources Research, %%
%%       44, W00B09, doi:10.1029/2007WR006720                                         %%
%%   Vrugt, J.A., H.V. Gupta, W. Bouten and S. Sorooshian (2003), A Shuffled Complex  %%
%%       Evolution Metropolis algorithm for optimization and uncertainty assessment   %%
%%       of hydrologic model parameters, Water Resour. Res., 39 (8), 1201,            %%
%%       doi:10.1029/2002WR001642.                                                    %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.lik = 44;      % 13: Normal likelihood function with AR(2)-model of residuals
                        % 14: Generalized likelihood function (Obsolete)
                        % 16: Laplacian likelihood function with AR(2)-model of residuals
                        % 17: Skewed Student t likelihood function
                        % 44: Generalized likelihood function PLUS
                        % 45: Universal likelihood function
%% Please refer to Vrugt et al. (2022) for theory on likelihood functions

%% Define name of function (.m file) for posterior exploration
Func_name = 'hymod';

%% Parameter names of hmodel
par_names = {'S_{\rm u,max}','\beta','\alpha','K_{\rm s}','K_{\rm f}'};

% index:          1     2     3     4     5
% parname:      Sumax  beta  alfa   Ks    Kf
fpar_mod   =  [  200   1.00  0.5   0.01  0.50 ];
parmin_mod =  [   50   0.10  0.0   1e-5  0.25 ];
parmax_mod =  [ 1000   10.0  1.0   0.25  5.00 ];

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';         % Latin hypercube sampling
Par_info.boundhandling = 'reflect'; % Explicit boundary handling
%Par_info.norm = 1;                 % Work in normalized space

global LV       %#ok Global for likelihood functions 13,14,16,17,44,45
% Be careful with its use
switch DREAMPar.lik
    case 13 %% Normal distribution
        % index:           1   2    3   4
        % parname:        s0  s1  phi1 phi2
        fpar_nuis   =  [  0.1  0    0   0 ];
        parmin_nuis =  [   0   0    0   0 ];
        parmax_nuis =  [   1   1    1   1 ];
        LV.filename = 'Normal'; id_nuis = [3];
    case 14 %% GL (OBSOLETE, DO NOT USE)
        % index:           1    2    3   4   5   6    7    8    9  10  11
        % parname:       std0 std1 beta xi  mu1 phi1 phi2 phi3 phi4 K lambda
        fpar_nuis   =  [  0.1   0    0   1   0   0    0    0    0   0   1  ];
        parmin_nuis =  [   0    0   -1  0.1  0   0    0    0    0   0  0.1 ];
        parmax_nuis =  [   1    1    1  10  100  1    1    1    1   1   1  ];
        LV.filename = 'GL'; id_nuis = [1 2 3 4 6];
    case 16 %% Laplace distribution
        % index:           1   2    3
        % parname:        s0  s1  phi1
        fpar_nuis   =  [  0.1  0    0 ];
        parmin_nuis =  [   0   0    0 ];
        parmax_nuis =  [   1   1    1 ];
        LV.filename = 'Laplace'; id_nuis = [ 3 ];
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
        LV.filename = 'UL'; id_nuis = [1 3 4 5 6];
end

%% Load data of Leaf River
load bound.txt;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Define warm-up period, first and last index of measured discharge record
plugin.data.wmp = 64; idx_t1 = 65; idx_t2 = 3717;
%% FOR TESTING ONLY
%idx_t2 = 1000;
idx_p = idx_t1 - plugin.data.wmp : idx_t2;
plugin.data.P = sum(bound(idx_p,6:9),2); plugin.data.Ep = bound(idx_p,5);
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
%% Measured streamflow data in mm/d (from m3/s with area of 1944)
F = 1944 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);
Meas_info.Y = 1/F * bound(idx_t1:idx_t2,4); %%
Meas_info.sigma2 = 'nonconstant';   % Non-constant measurement error variance
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

nuis_names = nuis_var_names(DREAMPar);          % Extract names nuis. variables
names = [par_names,nuis_names];                 % Combine with parameter names
n_modpar = numel(fpar_mod);                     % # model parameters
LV.id_vpar = [1:n_modpar , id_nuis + n_modpar]; % Model parameter + nuisance variable selection
LV.fpar = [fpar_mod fpar_nuis];                 % Merge default values parameters and nuisance variables
parmin = [parmin_mod parmin_nuis];
parmax = [parmax_mod parmax_nuis];              % Merge the parameter names and nuisance variables
Par_info.min = parmin(LV.id_vpar);
Par_info.max = parmax(LV.id_vpar);              % Min/max values of parameter selection
Par_info.names = names(LV.id_vpar);             % Names of parameters and nuisance variables
DREAMPar.d = size(Par_info.min,2);              % # parameters

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_kzs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                        % # Markov chains
        DREAMPar.T = 2500;                      % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                         % # Markov chains
        DREAMPar.T = 10000;                     % # generations
    case {'dream_kzs'}
        DREAMPar.N = 3;                         % # Markov chains
        DREAMPar.T = 5000;                      % # generations
        DREAMPar.M = 24;                        % # archive samples Kalman jump
        DREAMPar.a_1 = 0;                       % Activate Kalman jump immediately
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 250*ones(1,DREAMPar.d);    % # discrete steps
end

%% Optional settings
options.modout = 'no';               % Return model simulations samples?
%options.parallel = 'yes';            % Run each chain on a different core
options.parallel = 'no';            % Run each chain on a different core
options.save = 'yes';                % Save workspace DREAM during run
options.modout = 'yes';

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);
