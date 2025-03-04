function [dSdt] = InterceptionModelEquations(t,S,flag,P_int,E0_int,a,b,c,d)

% Note: P and E0 are vectors with rainfall and potential evapotranspiration
% They are thus time varying. Each iteration we thus need to derive their
% values from interpolation. At current time, rainfall is equal to: 
% P_int = max ( interp1 ( input(:,1) , input(:,2) , t ) , 0 ); % --> first column of input is time, second is rainfall, third is E0
% Same for E0
% E0_int = max ( interp1 ( input(:,1) , input(:,3) , t ) , 0 ); 

% NOTE: I assume average forcing values for time period between time t and t+1
% Advantage: interp1 is slow --> saves time - not affecting much numerics in this example

% Calculate interception (unknown parameter x rainfall intensity)
I = a * P_int;              % --> in mm/day

% Calculate drainage (only if storage larger than storage capacity)
if S > c,
    D = b*( S - c);         % --> in mm/day
else
    D = 0;
end

% Calculate evaporation
E = d * E0_int * S / c;      % --> in mm/day

% Now calculate the change in storage 
dSdt = I - D - E;           % --> in mm/day