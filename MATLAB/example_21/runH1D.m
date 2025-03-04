function sim_hydrus = runH1D(x,data_hydrus)
% Run the HYDRUS-1D model with new parameters and initial and boundary conditions and get simulated time series of soil water contents

% Update to values of "x" the parameter values in HYDRUS input file SELECTOR.IN
ModifySelectorIn(x(1:5))

% Modify file PROFILE.DAT using latest value of bottom boundary condition
data_hydrus.initial(:,3) = x(6);
ModifyProfileDat(data_hydrus);

% Modify ATMOSPH.IN
data_hydrus.boundcon(:,7) = x(6);
ModifyAtmosphIn(data_hydrus);

% Run HYDRUS-1D in dos command
[status,output] = dos('H1D_CALC.EXE');

% Read OBS_NODE.OUT
sim_hydrus = ReadObsNodeOut;