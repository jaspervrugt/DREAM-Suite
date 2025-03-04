function [ rho ] = ABC_func(x)
% Simple one-dimensional example to illustrate ABC

% Example function
%if rand < 1/2
%    rho = normrnd(x(1),1/10); 
%else
%    rho = normrnd(x(1),1);
%end;

% Draw 100 numbers from standard normal distribution
Y = normrnd(x(1),1,1,100);

% Define distance function
if rand < 1/2
    % Take the mean of the absolute values of x
    rho = abs(mean(Y));
else
    rho = abs(Y(1));
end