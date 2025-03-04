function nbins = calcnbins(x, method, minb, maxb)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the "ideal" number of bins using different methods
%
% nbins = CALCNBINS(x, method) calculates the ideal number of bins to use
%   in a histogram using a choice of methods
% SYNOPSIS
%       nbins = calcnbins(X)
%       nbins = calcnbins(X,method);
%       nbins = calcnbins(X,minb);
%       nbins = calcnbins(X,minb,maxb);
%
% REQUIRED INPUT ARGUMENTS
%   x       vector of data points
% OPTIONAL INPUT ARGUMENTS
%   method  string with choice of method        [default: method = 'middle']
%       'fd'        A single integer is returned, and calcnbins uses the
%                   Freedman-Diaconis method, based upon the inter-quartile
%                   range and number of data points
%       'scott'     A single integer is returned, and calcnbins uses
%                   Scott's method, based upon the sample standard
%                   deviation and number of data
%       'sturges'   A single integer is returned, and calcnbins uses
%                   Sturges' method, based upon the number of data
%       'middle'    A single integer is returned. calcnbins uses all three
%                   methods, then picks the middle (median) value
%       'all'       A structure is returned with fields 'fd', 'scott' and
%                   'sturges', each containing the calculation from the
%                   respective method
%   minb    smallest acceptable number of bins  [default: minb = 1]
%   maxb    largest acceptable number of bins   [default: minb = inf]
%
% Literature
%   D. Freedman, and P. Diaconis (1981), On the histogram as a
%       density estimator: L2 theory, Zeitschrift für
%       Wahrscheinlichkeitstheorie und verwandte Gebiete, 57 (4), pp. 453-476.
%   D.W. Scott (1979), On optimal and data-based histograms, Biometrika,
%       66 (3), pp. 605-610.
%   H.A. Sturges (1926). "The choice of a class interval". Journal of the
%       American Statistical Association, pp. 65-66.
%
% Warning/error is produced if
%   method is not written as full name
%   x has complex numbers (imaginary component will be ignored)
%   x is a matrix or multidimensional array
%
% Note, the number of bins (nbins) is just an estimate - not necessarily a
% definite choice!
%
% EXAMPLE
% y = randn(1000,1);
% nb = calcnbins(y, 'all')
%    nb =
%             fd: 23
%          scott: 18
%        sturges: 11
% calcnbins(y)
%    ans =
%        18
% subplot(3, 1, 1); hist(y, nb.fd);
% subplot(3, 1, 2); hist(y, nb.scott);
% subplot(3, 1, 3); hist(y, nb.sturges);
% y2 = rand(100,1);
% nb2 = calcnbins(y2, 'all')
%    nb2 =
%             fd: 5
%          scott: 5
%        sturges: 8
% hist(y2, calcnbins(y2))
%
% Author: Richard Cotton    Date: 2008/10/24    Version 1.5
%         Jasper Vrugt      Date: 2021/11/18    Version 2.0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input checking
narginchk(1, 4);
if ~isnumeric(x) && ~islogical(x)
    error('calcnbins:invalidX', 'The x argument must be numeric or logical.')
end
if ~isreal(x)
    x = real(x);
    warning('calcnbins:complexX', 'Imaginary parts of x will be ignored.');
end
% Ignore dimensions of x.
if ~isvector(x)
    x = x(:);
    warning('calcnbins:nonvectorX', 'x will be coerced to a vector.');
end
nanx = isnan(x);
if any(nanx)
    x = x(~nanx);
    warning('calcnbins:nanX', 'Values of x equal to NaN will be ignored.');
end
if nargin < 2 || isempty(method)
    method = 'middle';
end
if ~ischar(method)
    error('calcnbins:invalidMethod', 'method must be a char array.');
end
validmethods = {'fd','scott','sturges','all','middle'};
methodmatches = strncmpi(validmethods,method,10);
nmatches = sum(methodmatches);
if nmatches < 1
    error('calnbins:unknownMethod', 'The method specified is unknown');
end
method = validmethods{methodmatches};
if nargin < 3 || isempty(minb)
    minb = 1;
end
if nargin < 4 || isempty(maxb)
    maxb = Inf;
end

% Switch among methods
switch method
    case 'fd'
        nbins = calc_fd(x);
    case 'scott'
        nbins = calc_scott(x);
    case 'sturges'
        nbins = calc_sturges(x);
    otherwise
        nbins(1) = calc_fd(x);
        nbins(2) = calc_scott(x);
        nbins(3) = calc_sturges(x);
end
nbins = confine2range(nbins, minb, maxb);

if strcmp(method,'middle')
    nbins = median(nbins);
elseif strcmp(method,'all')
    nbins = struct('fd',nbins(1),'scott',nbins(2),'sturges',nbins(3));
end

end

%% Secondary functions
function nbins = calc_fd(x)
h = iqr(x); % diff(prctile0(x, [25; 75])); % inter-quartile range
if h == 0
    h = 2*median(abs(x-median(x))); % twice median absolute deviation
end
if h > 0
    nbins = ceil((max(x)-min(x))/(2*h*length(x)^(-1/3)));
else
    nbins = 1;
end
%    nbins = confine2range(nbins, minb, maxb);
end

function nbins = calc_scott(x)
h = 3.5*std(x)*length(x)^(-1/3);
if h > 0
    nbins = ceil((max(x)-min(x))/h);
else
    nbins = 1;
end
%    nbins = confine2range(nbins, minb, maxb);
end

function nbins = calc_sturges(x)
nbins = ceil(log2(length(x)) + 1);
%    nbins = confine2range(nbins, minb, maxb);
end

function y = confine2range(x, lower, upper)
y = ceil(max(x, lower)); y = floor(min(y, upper));
end

% % function y = prctile0(x, prc)
% % % Simple version of prctile that only operates on vectors, and skips
% % % the input checking (In particluar, NaN values are now assumed to
% % % have been removed.)
% % lenx = length(x);
% % if lenx == 0
% %     y = [];
% %     return
% % end
% % if lenx == 1
% %     y = x;
% %     return
% % end
% % 
% %     function foo = makecolumnvector(foo)
% %         if size(foo, 2) > 1
% %             foo = foo';
% %         end
% %     end
% % 
% % sortx = makecolumnvector(sort(x));
% % posn = prc.*lenx/100 + 0.5;
% % posn = makecolumnvector(posn);
% % posn = confine2range(posn, 1, lenx);
% % y = interp1q((1:lenx)', sortx, posn);
% % end