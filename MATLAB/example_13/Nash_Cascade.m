function [SimRR] = Nash_Cascade(k);
% Nash-Cascade unit hydrograph -- series of three linear reservoirs

% Store local variables in memory
persistent n Time maxT Precip F

if isempty(n)
    
    % Load the French Broad data
    daily_data = load('03451500.dly');
    % Define maximum time
    maxT = 365;
    % Define the precipitation
    Precip = daily_data(1:maxT,4);
    % Area factor to translate Nash-Cascade output from mm/d to m3/s
    F = 767 * (1000 * 1000 ) / (1000 * 60 * 60 * 24);
    % Define number of linear reservoirs
    n = 3;
    % Define Time
    Time = [ 1 : maxT ];
   
end

% -------------------------------------------------------------------------
%                           Model script
% -------------------------------------------------------------------------

if k < 1
    % Write to screen
    disp('Recession constant < 1 day --> numerical errors possible')
end

% Define help matrix
A = zeros(maxT,maxT);

% Calculate unit hydrograph
IUH = 1 / ( k * gamma(n) )  * ( Time / k ).^(n-1) .* exp( -Time / k );

% Now loop over time
for t = 1 : maxT
    
    % Define idx
    idx = [ 1 : maxT - (t - 1) ];
    
    % Calculate flow
    A(t,t:maxT) = Precip(t) * IUH(idx);
    
end

% Now determine total flow
SimRR = F * sum(A)';