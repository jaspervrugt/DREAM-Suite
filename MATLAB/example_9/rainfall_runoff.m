function Y = rainfall_runoff(x)
% Rainfall runoff model 

% Store local variables in memory
persistent plugin

% Load the data and define local variables - only once
if isempty(plugin)
    % Load the French Broad data
    daily_data = load('03451500.dly');
    % First two years are warm-up
    idx = (731:size(daily_data,1) + 1)';
    % Define the PET and precipitation.
    data.P = daily_data(:,4); data.Ep = daily_data(:,5);
    %% Initial conditions
    y0 = 1e-5 * ones(5,1);
    %% Integration options
    options.InitialStep = 1;                 % initial time-step (d)
    options.MaxStep     = 1;                 % maximum time-step (d)
    options.MinStep     = 1e-5;              % minimum time-step (d)
    options.RelTol      = 1e-3;              % relative tolerance
    options.AbsTol      = 1e-3*ones(5,1);    % absolute tolerances (mm)
    options.Order       = 2;                 % 2nd order accurate method (Heun)
    %% Running time
    tout = 0 : size(data.P,1);
    %% Dictionary plugin
    plugin = struct('idx',idx,'data',data,'y0',y0,'options',options,'tout',tout);
    
end
plugin.idx'
%% Run model
Y = hmodel(x,plugin);