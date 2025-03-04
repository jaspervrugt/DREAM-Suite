function s = logmeanexp(x,dim)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%LOGMEANEXP Returns log(mean(exp(a),dim)) while avoiding numerical        %
% underflow                                                               %
%                                                                         %
% dim = 1 for taking the mean along each column of a (default)            %
% dim = 2 for taking the mean along each row of a                         %
%                                                                         %
% © Modified by Jasper A. Vrugt, Aug. 2015                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 2, dim = 1; end
s = logsumexp(x, dim) - log(size(x, dim));

end