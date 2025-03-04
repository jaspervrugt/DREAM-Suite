function Y = rainfall_runoff(x,plugin)
% Run hmodel using C++ source code. The vector x stores the model
% parameters, and possibly also rainfall multipliers

mult = ones(plugin.Tmax,1);     % Unpack the multipliers
for j = 1 : size(plugin.id,1)   % Extract values of multipliers
    mult(plugin.id(j,1):plugin.id(j,2)) = x(plugin.nmod + j);
end
plugin.data.P = mult.* plugin.P;    % New hyetograph

mult'

%% Assign parameters
func.Imax  = x(1);      % interception storage capacity (mm)
func.Sumax = x(2);      % unsaturated zone storage capacity (mm)
func.Qsmax = x(3);      % maximum percolation rate (mm/d)
func.aE    = x(4);      % evaporation coefficient
func.aF    = x(5);      % runoff coefficient
func.aS    = 1e-6;      % percolation coefficient
func.Kf    = x(6);      % fast-flow response time (d)
func.Ks    = x(7);      % slow-flow response time (d)

%% Run hmodel in MATLAB or C
Y = hmodel(x,plugin); 
