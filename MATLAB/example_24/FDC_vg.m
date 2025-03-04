function q = FDC_vg(x,plugin)
% This function returns the streamflow according to van Genuchten soil 
% water characteristic. Exceedance probabilities are used as input

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Written by Jasper A. Vrugt                          %
%                   University of California Irvine                       %
%                           December 2014                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameter values
a = x(1); b = x(2); c = x(3);
% Exceedance probability in --> discharge out
q = (1/a) * (plugin.E.^(-1/c)-1).^(1/b);