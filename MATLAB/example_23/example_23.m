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
%% Example 23: SAC-SMA model: Limits of Acceptability                                 %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. and K.J. Beven (2018), Embracing equifinality with efficiency:       %%
%%       Limits of Acceptability sampling using the DREAM_{(LOA)} algorithm,  Journal %%
%%       of Hydrology,  559 , pp. 954-971, doi:10.1016/j.jhydrol.2018.02.026          %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Vrugt, J.A., Gupta H.V., Dekker, S.C., Sorooshian, S., Wagener, T. and W.        %%
%%       Bouten (2006), Application of stochastic parameter optimization to the       %%
%%       Sacramento Soil Moisture Accounting model, Journal of Hydrology, 325, pp.    %%
%%       288-307                                                                      %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.d = 13;                        % Dimension of the problem
DREAMPar.thinning = 1;                  % Only store each 5th sample
DREAMPar.lik = 23;                      % Limits of acceptability likelihood

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Initial sample from Latin hypercube sampling

%% Give the parameter ranges (minimum and maximum values)
Par_info.names = {'S1_{\rm Fmax}','S1_{\rm Tmax}','S2_{\rm FPmax}','S2_{\rm FSmax}', ...
    'S2_{\rm Tmax}','\alpha','\psi','k_{\rm i}','\kappa','\nu_{\rm p}','\nu_{\rm s}', ...
    'a_{\rm cmax}','k_{\rm f]'};
% Define parameter ranges
% index:            1       2        3      4        5       6      7     8      9   10   11     12   13
%% parname:      S1_Fmax S1_Tmax S2_FPmax S2_FSmax S2_Tmax alfa    psi   ki   kappa  nu_p nu_s a_cmax kf
%% parname:       uzfwm   uztwm    lzfpm   lzfsm    lztwm  zperc  rexp  uzk  pfree  lzpk lzsk   acm   kf
fpar_mod     =  [ 200      100      200    100      100    100     2    10.0  0.02  0.10 0.02   0.1  0.5 ];
Par_info.min =  [  10      10       10      10       10     1      1    0.01  0.05  1e-3 1e-3   0.05 0.0 ];
Par_info.max =  [ 500      500     1000   1000      500    250      5   1000  0.95  0.25 0.25   0.95 1.0 ];
Par_info.boundhandling = 'reflect';     % Explicit boundary handling

%% Define modelName
Func_name = 'sacsma';

%% Define plugin structure for hmodel
load bound.txt;  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Define warm-up period, first and last index of measured discharge record
plugin.data.wmp = 64; idx_t1 = 65; idx_t2 = 795; %3717;
%% FOR TESTING ONLY
idx_p = idx_t1 - plugin.data.wmp : idx_t2;
plugin.data.P = sum(bound(idx_p,6:9),2); plugin.data.Ep = bound(idx_p,5);
plugin.tout = 0:1:size(plugin.data.Ep,1);
plugin.idx = 1 + ( plugin.data.wmp : numel(idx_p) )';
%% Initial conditions
plugin.y0 = 1e-5 * ones(9,1);
%% Integration options
plugin.options.InitialStep = 0.2;   % initial time-step (d)
plugin.options.MaxStep     = 1;     % maximum time-step (d)
plugin.options.MinStep     = 1e-5;  % minimum time-step (d)
plugin.options.RelTol      = 1e-3;  % relative tolerance
plugin.options.AbsTol      = 1e-3;  % absolute tolerances (mm)
plugin.options.Order       = 2;     % 2nd order accurate method (Heun)
%% Measured streamflow data in mm/d (from m3/s with area of 1944)
F = 1944 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);
Meas_info.S = 1/F * bound(idx_t1:idx_t2,4); %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Number of parameters
d = numel(Par_info.min);

%% Optional settings
options.modout = 'yes';             % Return model simulations (yes/no)?
options.epsilon = 1.2*Meas_info.S;  % Limits of acceptability

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 5000;                          % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 25000;                         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 500*ones(1,DREAMPar.d);        % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options);
