function [ Y ] = lotka_volterra ( x )
% Lotka-Volterra model of population dynamics

persistent dydt y0 tout;    % Retain variables in memory after first call

if isempty(dydt)
    %% Initial conditions (initial population of prey and predator population)
    y0 = [30 4];
    %% Print time in years (start at t = 0 )
    tout = [ 0 : 1/12 : 20 ];
    % Define dydt --> the Lotka-Volterra equations differential equations
    % dydt = inline('[ alpha * y(1) - beta * y(1) * y(2) ; -gamma * y(2) + delta * y(1) * y(2) ]', ...
    %   't', 'y','alpha','beta','gamma','delta');
    % We apply a trick so we can use built-in ODE solver with inline function and additional parameters
    % We pass parameters as state variables (three to six) and specify
    % their partial derivatives as zero so their values do not change
%    dydt = inline('[ y(3)*y(1)-y(4)*y(1)*y(2) ; -y(5)*y(2)+y(6)*y(1)*y(2) ; zeros(4,1) ]', 't', 'y');
    dydt = @(t,y) [ y(3)*y(1)-y(4)*y(1)*y(2) ; -y(5)*y(2)+y(6)*y(1)*y(2) ; zeros(4,1) ];
end

%% Unpack the parameters of x in order of alpha, beta, gamma and delta and add to initial states
y = [ y0 x ]';
%% Now use ode45 time-variable integration method to solve system of two equations (monthly print step)
[ t , Y ] = ode45 ( dydt , tout , y ); 
% Remove the first row of Y (initial state at t = 0) and return output as one vector
Y = Y(2:end,1:2); Y = Y(:);