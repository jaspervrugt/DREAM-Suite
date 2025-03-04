function sim_hydrus = ReadObsNodeOut
% Read text file containing the time series of simulated soil water contents

% Open OBS_NODE.OUT
fid = fopen('OBS_NODE.OUT');

% Go to data section
flag = [];
while isempty(flag)
	str = fgetl(fid);
	flag = findstr(str,'time');
end
% Read simulated soil water contents
cols = 4; rows = Inf;
data = fscanf(fid,'%f',[cols rows]);
% close file
fclose(fid);

% Transpose data
data = data';

% Store simulated soil water contents in "sim_hydrus"
sim_hydrus.hoy = data(:,1); sim_hydrus.water = data(:,3);