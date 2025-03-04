function [data_hydrus] = Load_data
% Load observational data_hydrus and define initial and boundary conditions

% Load observational for hydrus
load meas.mat
data_hydrus.water = meas.water(meas.waterind);
data_hydrus.hoy = meas.hoy(meas.waterind);
data_hydrus.ObsNode = 1;

% Provide data_hydrus needed to modify initial condition
load ProfileDat.mat
data_hydrus.initial = initial_conditions;

% Provide data_hydrus needed to modify the lower boundary condition
load AtmosphIn.mat
data_hydrus.boundcon = boundcon;