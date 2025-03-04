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
%% Example 27: Simultaneous Optimization and Data Assimilation (SODA) of the single   %%
%%             layer interception model of Bouten (1992) using the Ensemble Kalman    %%
%%             Filter coupled with MCMC-based parameter estimation                    %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten (2005),     %%
%%       Improved treatment of uncertainty in hydrologic modeling: Combining the      %%
%%       strengths of global optimization and data assimilation, Water Resources      %%
%%       Research, 41, W01017, doi:10.1029/2004WR003059                               %%
%%   Vrugt, J. A., S.C. Dekker, and W. Bouten, Identification of rainfall             %%
%%       interception model parameters from measurements of throughfall and forest    %%
%%       canopy storage (2003), Water Resources Research, 39(9), 1251,                %%
%%       doi:10.1029/2003WR002013                                                     %%
%%                                                                                    %%
%% Check video lectures on my website: faculty.sites.uci.edu/jasper/teaching          %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% The SODA method uses an inner state estimation loop with the Ensemble Kalman
% Filter (EnKF) inside an outer parameter estimation method using MCMC simulation with
% the DREAM algorithm. Thus, MCMC simulation proposes a candidate point. This proposal
% enters the EnKF and produces a vector of mean state forecasts at each time t. This
% results in a likelihood of the parameter values, which is used by the MCMC method to
% generate a new candidate vector.
%
% This example does not appear in Vrugt et al. (2005) but is used in my lectures on
% data assimilation and the SODA methodology. The measured forcing data is shifted
% in comparison to the forcing data used in the interception model. This makes it
% impossible to describe accurately the observed storage data. Data assimilation will
% be essential in removing the impact of the forcing data errors on the simulated
% storage (= state variable) by the single-layer interception model. This guarantees
% a better compliance between the inferred parameter estimates and their true values.

% I originally developed SODA to better treat structural and input data errors during
% Bayesian parameter estimation. If such errors are not treated properly then model
% parameters will compensate for such systematic errors. The hypothesis in the SODA
% paper was that the time series of state updates (= innovations) should contain
% valuable information about epistemic (= model structural) errors.

% I recommend that you check the following lecture to understand in more detail
% the concepts, theory, and implementation of data assimilation:
% https://www.youtube.com/watch?v=l5cb4ONIB-c&feature=youtu.be&disable_polymer=true
% including methods such as SODA.

%% Problem settings defined by user
DREAMPar.d = 5;                         % Dimension of the problem
DREAMPar.lik = 11;                      % Execute Gaussian likelihood!

%% Determine the parameter ranges (depends on data, u and copula choice
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling
Par_info.names = {'a','b','c','d','\phi_{\rm mod}'};
Par_info.min = [ 0.5 1    0.5 0.1  0 ]; % minimum parameter values
Par_info.max = [ 1.0 1000 5.0 2.0  1 ]; % maximum parameter values

%% Define name of function (.m file) for posterior exploration
Func_name = 'ENKF_Storage';

%% Load the data ( reference run without measurement error )
all_data = readmatrix('Data_Canopy_Storage.xlsx','Sheet','Sheet1', ...
    'DataRange','B4:E84');

%% Unpack the observed storage
obsStorage = all_data(:,4);

%% Now define the measurement error matrix (= scalar)
R = diag(0.01);

%% Single column of "measured states" (= storage) for likelihood function
Meas_info.Y = obsStorage;

%<><><><><><><><><><><><><><><> ENKF SETUP <><><><><><><><><><><><><><><><>
plugin.obsTime = all_data(:,1); % Unpack observation times
plugin.Precip = all_data(:,2);  % Unpack precipitation time series
plugin.Potevap = all_data(:,3); % Unpack potential evaporation
plugin.R = R;                   % Measurement error covariance matrix
plugin.C = 0.2;                 % Model error covariance matrix
plugin.m = 1;                   % # state variables
plugin.N = 50;                  % # ensemble members in EnKF? (quick result)
plugin.Y = obsStorage';         % Observed storage measurements
plugin.NT = size(plugin.Y,2);   % # measurement times
plugin.options = odeset;        % options for ODE solver
plugin.options.InitialStep = ...        % Initial time step
    1e-2 * max(diff(plugin.obsTime));
plugin.options.MaxStep = ...            % Maximum time step
    max(diff(plugin.obsTime));
%<><><><><><><><><><><><><><> END ENKF SETUP <><><><><><><><><><><><><><><>

% Additional DREAM Package options
options.modout = 'yes';         % Return model simulations (yes/no)?
options.parallel = 'yes';       % Run each chain on a different core
options.save = 'yes';           % Save workspace DclREAM during run

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_dzs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                        % # Markov chains
        DREAMPar.T = 1000;                      % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                         % # Markov chains
        DREAMPar.T = 3000;                      % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 900*ones(1,DREAMPar.d);    % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);
