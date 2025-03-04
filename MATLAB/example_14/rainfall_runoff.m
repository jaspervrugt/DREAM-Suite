function S_mod = rainfall_runoff(x,plugin)
% Rainfall runoff model (= hmodel) using C++ executable with time variable
% integration time step. This C-code was written by Gerrit Schoups and
% modified by J.A. Vrugt. Similar implementations are available for Hymod,
% GR4J, NAM, HyMOD and the SAC-SMA model. 

%% Store local variables in memory
% persistent data y0 options tout
% 
%% Load the data and define local variables - only once
% if isempty(tout)
%     % Load the French Broad data
%     plugin.daily_data = load('03451500.dly');
%     % First two years are warm-up
%     plugin.data.idx = [731:size(daily_data,1)]';
%     % Define the PET and precipitation.
%     plugin.data.P = daily_data(:,4); plugin.data.Ep = daily_data(:,5);
%     %% Initial conditions
%     plugin.y0 = 1e-5 * ones(5,1);
%     %% Integration options
%     plugin.options.InitialStep = 1;                 % initial time-step (d)
%     plugin.options.MaxStep     = 1;                 % maximum time-step (d)
%     plugin.options.MinStep     = 1e-5;              % minimum time-step (d)
%     plugin.options.RelTol      = 1e-3;              % relative tolerance
%     plugin.options.AbsTol      = 1e-3*ones(5,1);    % absolute tolerances (mm)
%     plugin.options.Order       = 2;                 % 2nd order accurate method (Heun)
%     %% Simulation time
%     plugin.tout = 0 : size(plugin.data.P,1);
% end


%% Run model
Y = hmodel(x,plugin);
Y(1:3)
Y(end-2:end)
% Compute watershed signatures from simulated discharge record
S_mod = calc_signatures(Y)';