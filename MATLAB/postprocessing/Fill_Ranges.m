function [h] = Fill_Ranges(x,y_low,y_up,color,ax)
% Fill the ranges with a given color
% 
% SYNOPSIS: bounds(x,ylow,yupper);
%           bounds(x,ylow,yupper,color);
%
% Input:    x = vector with "x" values
%           ylow = vector with lower range values
%           yupper = vector with upper range values
%           color = filling color
%
if nargin < 5
    ax = gca;
end
    
% Now create a vector x1
X = [x(:); flipud(x(:)); x(1)];

% And corresponding matrix y1
Y = [y_up(:); flipud(y_low(:)); y_up(1)];

% Now fill area with "fill" function 
h = fill(ax,X,Y,color);

% Set color
set(h,'edgecolor',color);