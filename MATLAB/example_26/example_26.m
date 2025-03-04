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
%% Example 26: Simultaneous Optimization and Data Assimilation (SODA) of the Lorenz   %%
%%             1963 model using the Ensemble Kalman Filter coupled with MCMC-based    %%
%%             parameter estimation                                                   %%
%%                                                                                    %%
%% Check the following papers                                                         %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, 75, pp. 273-316, doi:10.1016/j.envsoft.2015.08.013             %%
%%   Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten, Improved   %%
%%       treatment of uncertainty in hydrologic modeling: Combining the strengths of  %%
%%       global optimization and data assimilation, Water Resources Research, 41,     %%
%%       W01017, doi:10.1029/2004WR003059, 2005                                       %%
%%   Lorenz, E.N., Deterministic non-periodic flow, Journal of Atmospheric Sciences,  %%
%%       20, 130 141, 1963                                                            %%
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
% The present example is equivalent to case study 1 of Vrugt et al. (2005) and was
% used to benchmark the SODA methodology. This is a toy problem.

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
DREAMPar.d = 4;                         % Dimension of the problem
DREAMPar.lik = 11;                      % Normal likelihood!

%% Determine the parameter ranges
Par_info.initial = 'latin';             % Latin hypercube sampling
Par_info.boundhandling = 'reflect';     % Explicit boundary handling
Par_info.names = {'\sigma','\rho','\beta','\phi_{\rm mod}'};
Par_info.min = [ 0  0  0   0 ];         % minimum of parameter values
Par_info.max = [ 30 50 10  1 ];         % maximum of parameter values

%% Define name of function (.m file) for posterior exploration
Func_name = 'EnKF_Lorenz63';

%% Load the data ( reference run without measurement error )
Y_ref = readmatrix('Data_Lorenz_System.xlsx','Sheet','Sheet1', ...
    'DataRange','C4:E164');

%% How many measurement times?
n = size(Y_ref,1);

%% Now define the measurement error matrix (See Vrugt et al., 2005; Page 9)
R = diag([ 2 2 2 ]);

%% Draw n x 3 matrix of measurement errors for each observation time
meas_error = mvnrnd([0 0 0],R,n);

%% Generate n x 3 matrix with perturbed reference state
Y_meas = Y_ref + meas_error;

%% Single column of "measured states" for likelihood computation
Meas_info.Y = Y_meas(:);
Meas_info.sigma2 = 'constant';

%<><><><><><><><><><><><><><><> ENKF SETUP <><><><><><><><><><><><><><><><>
plugin.m = 3;                   % Number state variables Lorenz 1963 model
plugin.R = R;                   % m x m measurement error covariance matrix
plugin.C = [2.000 0.000 0.000;  % m x m model error covariance matrix
    0.000 12.13 0.000;  % (from Evensen, 1994)
    0.000 0.000 12.31]; % (see Vrugt et al., 2005; Page 9)
plugin.N = 10;                  % # ensemble members EnKF (quick result!)
plugin.Y = Y_meas';             % Measured states for EnKF updating step
plugin.NT = size(plugin.Y,2);   % # measurement times
plugin.dt = 0.25;               % Time step
plugin.options = odeset;            % options for ODE solver
plugin.options.InitialStep = 1e-3;  % Define initial time step
plugin.options.MaxStep = 1e-1;      % Define maximum time step
%<><><><><><><><><><><><><><> END ENKF SETUP <><><><><><><><><><><><><><><>

%% Additional options for DREAM Package
options.modout = 'yes';             % Return model simulations (yes/no)?
options.parallel = 'yes' ;          % Run each chain on a different core
options.save = 'yes';               % Save workspace DREAM during run

%% Define method to use {'dream','dream_zs','dream_d','dream_dzs','mtdream_zs'}
method = 'dream_zs';

switch method
    case {'dream','dream_d'}
        DREAMPar.N = 10;                        % # Markov chains
        DREAMPar.T = 2000;                      % # generations
    case {'dream_zs','dream_dzs','mtdream_zs'}
        DREAMPar.N = 3;                         % # Markov chains
        DREAMPar.T = 5000;                      % # generations
end

if strcmp(method,'dream_d') || strcmp(method,'dream_dzs')
    Par_info.steps = 1000*ones(1,DREAMPar.d);   % # discrete steps
end

%% Run DREAM-Suite
[chain,output,FX,Z,logL] = DREAM_Suite(method,Func_name,DREAMPar,...
    Par_info,Meas_info,options,[],plugin);
