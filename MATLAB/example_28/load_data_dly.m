function data = load_data_dly(ID_watershed);
% This function loads the data of a watershed if in dly format

% Open the file
evalstr = strcat('fid = fopen(''',ID_watershed,'.dly''',');');

% Now open the file
eval(evalstr);

% Check if file is existent
if fid == -1,
    % Return error
    error('Wrong ID_watershed');
else
    % Initialize data
    data = NaN(1e5,8);
    
    % Now read the data
    counter = 1;
    
    while 1 == 1,
        
        try
            
            % Now get the next line
            tline = fgetl(fid);
            
            % Now extract year, month, and day
            year = str2num(tline(1:4)); month = str2num(tline(5:6)); day = str2num(tline(7:8));
            
            % Now get all data
            dat = str2num(tline(9:end)); 
                
            % Now store (streamflow is third column of dat)
            data(counter,1:8) = [year month day dat];
            
            % Update counter
            counter = counter + 1;
            
        catch
            
            % Break the while loop
            fclose(fid); break;
            
        end
        
    end
    
    % Now remove remaining NaN values
    idx = find(isnan(data(:,1))); data = data ( 1 : idx(1)-1 , 1 : 8 );
    
end