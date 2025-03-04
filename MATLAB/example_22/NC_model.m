function [ Q ] = NC_model(x)
% Nash-Cascade unit hydrograph: Serial configuration 3 linear reservoirs

% Store local variables in memory
persistent T T_max P A

if isempty(T)
    T_max = 25;                         % Maximum time          [d]
    T = 1:T_max;                        % Define time           [d]
    P = [ 10 25 8 2 zeros(1,21) ]';     % Precipitation record  [mm/d]
                    % Define flow matrix    [mm/d]
end

k = x(1); n = x(2);  % Unpack recession constant and number of reservoirs

if k < 1
    % Write to screen
    disp('NC_model: Recession constant < 1 day --> numerical errors possible')
end

IUH = 1 / ( k * gamma(n) )  * ...       % Instantaneous unit hydrograph
    ( T / k ).^(n-1) .* exp( - T / k );
for t = 1 : T_max   % Loop over time
    id = 1 : T_max - (t - 1);           % Index
    A(t,t:T_max) = P(t) * IUH(id);      % Compute flow
end
Q = sum(A)';                            % Discharge: column vector
