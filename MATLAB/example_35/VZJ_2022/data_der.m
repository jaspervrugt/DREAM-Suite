function [t_meas,I_meas,dIdt,dtdI,n_data] = data_der(dat)
% Compute derivative from data

% Define plugin structure
n_data = size(dat,1);
% Store time and infiltration
t_meas = dat(1:n_data,1); I_meas = dat(1:n_data,2);
% define method
method = 'diff'; % diff
switch method
    case 'der'
        % Compute infiltration rate
        dIdt = diff(I_meas)./diff(t_meas);
        % Compute 1/infiltration rate
        dtdI = diff(t_meas)./diff(I_meas);
        % Return mean of data and cum. inf.
        t_meas = 1/2*(t_meas(1:n_data-1) + t_meas(2:n_data));
        I_meas = 1/2*(I_meas(1:n_data-1) + I_meas(2:n_data));
        % One soil produces 2 similar end times - remove this data point!
        % Soil_type = 228 (1st approach - problem) - maybe more - did not check later
        % SOil_type = 303 (2nd approach - problem)
        ii = find(isnan(dIdt));
        if isempty(ii)
            n_data = n_data - 1;
            alarm = 0;
        else
            dIdt = dIdt(1:ii-1);
            dtdI = dtdI(1:ii-1);
            t_meas = t_meas(1:ii-1);
            I_meas = I_meas(1:ii-1);
            n_data = n_data - 2;
            alarm = 1
        end
    case 'diff'
        % Must adjust the 
        dIdt = diff(I_meas); n_data = n_data - 1;
        dtdI = diff(t_meas);
end