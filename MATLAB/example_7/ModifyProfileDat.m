function ModifyProfileDat(data_hydrus)
% Modify text file containing initial condition

%global EXAMPLE_dir

% Open PROFILE.DAT
fid = fopen('PROFILE.DAT','r+');

% Go to profile section
flag = [];
while isempty(flag)
	str = fgetl(fid);
	flag = findstr(str,'Mat');
end

% Insert new initial condition
fseek(fid,0,'cof');
fprintf(fid,'%5.0f %15.6e %15.6e %3.0f %3.0f %15.6e %15.6e %15.6e %15.6e\n',data_hydrus.initial');
fseek(fid,0,'eof');

% Close file
fclose(fid);