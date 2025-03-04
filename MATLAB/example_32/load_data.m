function [Meas_info,plugin] = load_data
% This function loads the discharge data and defines the structures 
% Meas_info and plugin for use in DREAM package

load daily_data_maurer_1980_2008.mat;
t_daily = datetime(daily_data(:,3),daily_data(:,2),daily_data(:,1));
P_daily = daily_data(:,4);    % Precipitation in mm/d
Ep_daily = daily_data(:,5);   % Potential evaporation in mm/d
Q_daily = daily_data(:,6);    % Discharge in mm/d
%% Calibration period
t1 = datetime(1999,10,01); t2 = datetime(2004,09,30);
idx_t1 = find(t_daily==t1); idx_t2 = find(t_daily==t2);
Meas_info.Y = Q_daily(idx_t1:idx_t2);
Meas_info.sigma2 = 'nonconstant';   % Non-constant measurement error variance
%% Define plugin structure for hmodel
plugin.data.wmp = 1 * 365;
%% Do not change after this; if idx_t1, idx_t2 defined then all OK
idx_p = idx_t1 - plugin.data.wmp : idx_t2;
plugin.data.P = P_daily(idx_p); plugin.data.Ep = Ep_daily(idx_p);
%plugin.daily_data = daily_data(idx_p,:);
plugin.tout = 0:1:size(plugin.data.Ep,1);
plugin.idx = 1 + ( plugin.data.wmp : numel(idx_p) )';
plugin.n = size(Meas_info.Y,1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% hymod initial conditions
plugin.y0 = 1e-5 * ones(6,1);
%% hymod integration settings
plugin.options.InitialStep = 1;     % initial time-step (d)
plugin.options.MaxStep     = 1;     % maximum time-step (d)
plugin.options.MinStep     = 1e-3;  % minimum time-step (d)
plugin.options.RelTol      = 1e-3;  % relative tolerance
plugin.options.AbsTol      = 1e-3;  % absolute tolerances (mm)
plugin.options.Order       = 2;     % 2nd order accurate method (Heun)
