function [SimRR] = hymodFORTRAN(x) 
% Runs the FORTRAN HYMOD model and returns the simulated discharge

% Define local variables
persistent F MaxT

% Calculate area factor - only once
if isempty(F)
    % Define MaxT
    MaxT = 795;
    % Area factor to translate HYMOD output in mm/d to m3/s (calibration data); (area Leaf River is 1944 km2)
    F = 1944 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);
end

% Write the parameter values to a file Param.in
dlmwrite('Param.in',x,'delimiter',' ');

% Execute the model -- this model reads the current parameter values from Param.in
dos('HYMODsilent.exe');

try 
    % Load the output of the model 
    SimRR = F * load('Q.out'); SimRR = SimRR(65:MaxT);
catch
    % Assign bad values - model might have crashed
    SimRR = 0 * [65:MaxT]';
end