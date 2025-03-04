function [mx,my] = getmxmy(x,y)
% Determines line for quantile-quantile plot

% Interquartile range for x
q1x = prctile(x,25); q3x = prctile(x,75); dx = q3x - q1x; 
% Interquartile range for y
q1y = prctile(y,25); q3y = prctile(y,75); dy = q3y - q1y; 
% Determine slope
slope = dy./dx;
% Determine midpoint of IQR - for x and y
x_center = (q1x + q3x)/2; y_center = (q1y + q3y)/2;
% min and max of x
x_max = max(x); x_min = min(x);
% compute corresponding y-values: max
y_max = y_center + slope.*(x_max - x_center);
% compute corresponding y-values: min
y_min = y_center - slope.*(x_center - x_min);
% return arguments, mx and my
mx = [x_min; x_max]; my = [y_min; y_max];