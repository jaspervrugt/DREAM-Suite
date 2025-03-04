function [data_new,true_pars] = preprocess_data(data,Parameters,n_soil,...
    n_data,interpolation_method,interpolation_scheme,I_max)
% This function preprocesses the data and returns cell array with final
% data set and corresponding true values of Ks [cm/h] and S [cm/h^1/2]

if nargin == 6
    I_max = [];
end

% Initialize data_new and true_pars
data_new = cell(12,1); true_pars = nan(12,2); 

% Lets create the final data set from original cell array data of .mat file
for soil_type = 1:n_soil
    % Extract time (hours) versus cumulative infiltration (cm)
    dat = data{soil_type};
    % remove nan values (appear at end)
    dat = dat(~isnan(dat(:,2)),1:2);
    % remove duplicate values of time
    dat = dat([ 1 ; find(diff(dat(:,1)) > 0) + 1],1:2);
    % remove duplicate values of infiltration
    dat = dat([ 1 ; find(diff(dat(:,2)) > 0) + 1],1:2);
    % now switch among interpolation methods
    switch interpolation_method
        case 1 
            % What is t_meas, t_end, and I_meas?
            t_meas = dat(1:end,1); t_end = dat(end,1); I_meas = dat(1:end,2);
            % Divide time interval in n_data + 1 values
            switch lower(interpolation_scheme)
                case 'linear'
                    t_int = linspace(0,t_end,n_data+1);
                case 'logarithmic'
                    t_int = logspace(log10(1),log10(t_end+1),n_data+1) - 1;
                case 'square_root'
                    % Compute dT
                    dT = sqrt(t_end - 0)/n_data; 
                    % Now we need to take n_data steps with dT 
                    t_sqrt_t_int = 0:dT:sqrt(t_end);
                    % Now back transform to linear space
                    t_int = t_sqrt_t_int.^2; t_int(n_data+1) = t_end;
                otherwise
                      error('do not know this interpolation scheme')
            end
            % Interpolate I_meas to t_int to yield I_int
            I_int = interp1(t_meas,I_meas,t_int);
        case 2
            % What is t_meas and I_end?
            t_meas = dat(1:end,1); I_meas = dat(1:end,2);
            % Normalize measured infiltration with I_max
            I_norm = I_meas/I_max;
            % Get corresponding measurement times
            switch lower(interpolation_scheme)
                case 'linear'
                    I_norm_int = linspace(0,1,n_data+1);
                    % Note: now diff(I_int) should produce one and same value
                case 'logarithmic'
                    I_norm_int = logspace(log10(1),log10(2),n_data+1) - 1;
                    % Note: now diff(I_int) should increase with time
                otherwise
                    error('do not know this interpolation scheme')
            end
            % Now determine corresponding measurement times
            t_int = interp1(I_norm,t_meas,I_norm_int);
            % And denormalize cumulative infiltration values
            I_int = I_norm_int * I_max; 
    end
    % Combine (t_int,I_int) and remove time zero with I_int = 0;
    data_new{soil_type} = [t_int(2:n_data+1)' I_int(2:n_data+1)'];
    % Extract corresponding values of Ks [cm/h] and S [cm/h^(1/2)]
    true_pars(soil_type,1:2) = [Parameters{soil_type+2,9:10}];
end