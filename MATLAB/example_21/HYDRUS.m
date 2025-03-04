function sim_Y = HYDRUS(x,data_hydrus)
% Runs the HYDRUS model and computes the log-likelihood

% Store local variable (is initally empty after declaration)
persistent fid

% Only write the directory file (level_01.dir) one time in given directory
if isempty(fid)
    
    % Write level_01.dir (needed to execute HYDRUS-1D)
    fid = fopen('level_01.dir','w+'); fprintf(fid,'%s',pwd); fclose(fid);
    
end

% No transformation of parameters needed --> see example 7

try

    % Run HYDRUS-1D
	sim_hydrus = runH1D(x,data_hydrus);
	
	% Filter simulated water contents to return values at measurement times
    ind = zeros(size(data_hydrus.hoy));
	for i = 1:size(data_hydrus.hoy,1)
		ind(i) = find(sim_hydrus.hoy == data_hydrus.hoy(i));
	end

	% Return simulated soil moisture contents
    sim_Y = sim_hydrus.water(ind);
    
catch	
	
	% Return "bad" simulated value if HYDRUS-1D did not converge properly
    % --> this will result in a very low log-likelihood value
    sim_Y = zeros(size(data_hydrus.hoy));

end