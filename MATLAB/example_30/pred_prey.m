function Y = pred_prey(x,plugin)
% This function runs the 1-predator-1-prey model

% use built-in ode45 function to solve 2 differential equations from t = 0 to t = 335
[ ~ , u ] = ode45(@(t,u) PP(t,u,x),plugin.t,plugin.u0,plugin.ode_options);
% now return the simulated abundances of both species as one vector
Y = u(:);
end

% Secundairy function
function u_out = PP(t,u,pars)
% Solves the coupled ordinary differential equations of the 1-predator-1-prey population model

% Assign the states
x = u(1); y = u(2);

% intrinsic rate of prey natural increase
r = pars(1);
% proportionality constant linking prey mortality to the number of prey and predators
alfa = pars(2);
% mortality rate of predators
m = pars(3);
%  proportionality constant linking the increase in predators to the number of predators and prey
theta = pars(4);
% maximum number of preys that the environment can support (carrying capacity?)
K = 50;

% Now compute each of the two ordinary differential equations
dxdt = r * x * ( 1 - x/K ) - alfa * x * y;
dydt = -m * y + theta * x * y;

% Now assemble the dx/dt and dy/dt in one vector u
u_out = [ dxdt ; dydt ];
end