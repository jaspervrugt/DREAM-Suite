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
%% Example 14: Approx. Bayesian Comp. with watershed signatures & box-car likelihood  %%
%%             Note: FDCFIT toolbox has a much larger built-in class of FDC functions %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Sadegh, M., J.A. Vrugt, H.V. Gupta, and C. Xu (2016), The soil water             %%
%%       characteristic as new class of closed-form parametric expressions for the    %%
%%       flow duration curve, J. Hydrol., 535, pp. 438–456                            %%
%%   Lochbuhler, T., J.A. Vrugt, M. Sadegh, and N. Linde (2014), Summary statistics   %%
%%       from training images as prior information in probabilistic inversion,        %%
%%       Geophysical Journal International, 201, 157-171, doi:10.1093/gji/ggv008      %%
%%   Sadegh, M., and J.A. Vrugt (2014), Approximate Bayesian computation using Markov %%
%%       chain Monte Carlo simulation: DREAM_(ABC), Water Resources Research,         %%
%%       doi:10.1002/2014WR015386.                                                    %%
%%   Vrugt, J.A., and M. Sadegh (2013), Toward diagnostic model calibration and       %%
%%       evaluation: Approximate Bayesian computation, Water Resources Research, 49,  %%
%%       4335–4345, doi:10.1002/wrcr.20354.                                           %%
%%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson (2008), %%
%%       Treatment of input uncertainty in hydrologic modeling: Doing hydrology       %%
%%       backward with Markov chain Monte Carlo simulation, Water Resources Research, %%
%%       44, W00B09, doi:10.1029/2007WR006720                                         %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Problem settings defined by user
DREAMPar.d = 7;                         % Dimensionality of target
DREAMPar.lik = 22;                      % ABC with summary metrics

%% Provide information parameter space and initial sampling
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling

%% Parameter names of hmodel
Par_info.names = {'I_{\rm max}','S_{\rm u,max}','Q_{\rm s,max}',...
    '\alpha_{\rm e}','\alpha_{\rm f}','K_{\rm f}','K_{\rm s}'};
% parname:       Imax  Smax  Qsmax   alE   alF   Kfast  Kslow
Par_info.min = [ 0.5   10     0    1e-6   -10     0      0    ];    % If 'latin', min values
Par_info.max = [ 10   1000   100   100     10     10    150   ];    % If 'latin', max values

%% Define name of function (.m file) for posterior exploration
Func_name = 'rainfall_runoff';

%% Define the observed streamflow data
daily_data = load('03451500.dly'); %Meas_info.Y = daily_data(731:end,6);

% Define warm-up period and end of training data period
plugin.data.wmp = 730; idx_t1 = 731; idx_t2 = size(daily_data,1);
%% FOR TESTING ONLY
idx_p = idx_t1 - plugin.data.wmp : idx_t2;
%% Define the structure plugin -- as second argument to model function
plugin.daily_data = daily_data;
plugin.data.P = daily_data(idx_p,4); plugin.data.Ep = daily_data(idx_p,5);
plugin.tout = 0:1:size(plugin.data.Ep,1);
plugin.idx = 1 + ( plugin.data.wmp : numel(idx_p) )';
%% Initial conditions
plugin.y0 = 1e-5 * ones(5,1);
%% Integration options
plugin.options.InitialStep = 1;     % initial time-step (d)
plugin.options.MaxStep     = 1;     % maximum time-step (d)
plugin.options.MinStep     = 1e-3;  % minimum time-step (d)
plugin.options.RelTol      = 1e-3;  % relative tolerance
plugin.options.AbsTol      = 1e-3;  % absolute tolerances (mm)
plugin.options.Order       = 2;     % 2nd order accurate method (Heun)
%% Measured streamflow data in mm/d
Meas_info.Y = daily_data(idx_t1:idx_t2,6); %%
Meas_info.sigma2 = 'nonconstant';   % Non-constant measurement error variance

%% Now calculate summary metrics from the discharge data and define as prior distribution
Meas_info.S = calc_signatures( daily_data(731:end,6) )';

%% Optional settings
options.DB = 'no';                      % Diagnostic Bayes: ABC with summary metrics as prior
options.epsilon = [ 0.01 ];             % Epsilon of the noisy ABC implementation (illustration values)
options.parallel = 'yes';               % Run each chain on different core
options.modout = 'yes';                 % Store model simulations
options.save = 'yes';

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                            % # Markov chains
        DREAMPar.T = 10000;                         % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                             % # Markov chains
        DREAMPar.T = 10000;                         % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);       % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);
