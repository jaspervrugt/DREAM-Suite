function r = stblrnd(alpha,beta,gamma,delta,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% STBLRND alpha-stable random number generator                          %%
%%    R = STBLRND(alpha,beta,gamma,delta) draws a sample from the Levy   %%
%%    alpha-stable distribution with characteristic exponent alpha,      %%
%%    skewness beta, scale parameter gamma and location parameter delta. %%
%%    alpha, beta, gamma and delta must be scalars in following ranges:  %%
%%        0 < alpha <= 2                                                 %%
%%       -1 <= beta <= 1                                                 %%
%%        0 < gamma < inf                                                %%
%%     -inf < delta < inf                                                %%
%%                                                                       %%
%% R = STBLRND(alpha,beta,gamma,delta,M,N,...) or                        %%
%% R = STBLRND(alpha,beta,gamma,delta,[M,N,...]) returns M-by-N-by-array %%
%%                                                                       %%
%% References:                                                           %%
%% [1] Chambers, J.M., C.L. Mallows, and B.W. Stuck (1976), A Method for %%
%%     Simulating Stable Random Variables, JASA, Vol. 71, No. 354. pp.   %%
%%     340-344                                                           %%
%% [2] Weron, A. and R. Weron (1995), Computer Simulation of Levy        %%
%%     alpha-Stable Variables and Processes, Lecture Notes in Physics,   %%
%%     Vol. 457, pp. 379-392                                             %%
%%                                                                       %%
%% Modified by Jasper A. Vrugt, April 2017                               %%
%% University of California Irvine                                       %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

warning on; warning(['DREAM_PACKAGE:stblrnd:Code must be benchmarked. ',...
    'PDF from stblpdf does not match stblrnd for given values of ',...
    'alpha, beta, gamma and delta']);
% Check number of input arguments
if nargin < 4
    error(['DREAM_PACKAGE:stblrnd:TooFewInputs','Requires at least ', ...
        'four input arguments.']);
end
% Check parameter values
if alpha <= 0 || alpha > 2 || ~isscalar(alpha)
    error(['DREAM_PACKAGE:stblrnd:BadInputs',' "alpha" must be a ', ...
        'scalar which lies in the interval (0,2]']);
end
if abs(beta) > 1 || ~isscalar(beta)
    error(['DREAM_PACKAGE:stblrnd:BadInputs',' "beta" must be a ', ...
        'scalar which lies in the interval [-1,1]']);
end
if gamma < 0 || ~isscalar(gamma)
    error(['DREAM_PACKAGE:stblrnd:BadInputs',' "gamma" must be a ', ...
        'non-negative scalar']);
end
if ~isscalar(delta)
    error('DREAM_PACKAGE:stblrnd:BadInputs',' "delta" must be a scalar');
end
% Determine size of output
[err,sizeOut] = genOutsize(4,alpha,beta,gamma,delta,varargin{:});
if err > 0
    error(['DREAM_PACKAGE:stblrnd:InputSizeMismatch',...
        'Size information is inconsistent.']);
end

% <><><><><><><><><><><><><><> Draw sample <><><><><><><><><><><><><><><><>

% Check if simple case, if so be quick, if not do the general algorithm
if alpha == 2                           % Gaussian distribution
    r = sqrt(2) * randn(sizeOut);
elseif alpha==1 && beta == 0            % Cauchy distribution
    r = tan( pi/2 * (2*rand(sizeOut) - 1) );
elseif alpha == .5 && abs(beta) == 1    % Levy distribution (= Pearson V)
    r = beta ./ randn(sizeOut).^2;
elseif beta == 0                        % Symmetric alpha-stable
    V = pi/2 * (2*rand(sizeOut) - 1);
    W = -log(rand(sizeOut));
    r = sin(alpha * V) ./ ( cos(V).^(1/alpha) ) .* ...
        ( cos( V.*(1-alpha) ) ./ W ).^( (1-alpha)/alpha );
elseif alpha ~= 1                       % General case, alpha not 1
    V = pi/2 * (2*rand(sizeOut) - 1);
    W = - log( rand(sizeOut) );
    const = beta * tan(pi*alpha/2);
    B = atan( const );
    S = (1 + const * const).^(1/(2*alpha));
    r = S * sin( alpha*V + B ) ./ ( cos(V) ).^(1/alpha) .* ...
        ( cos( (1-alpha) * V - B ) ./ W ).^((1-alpha)/alpha);
else                                    % General case, alpha = 1
    V = pi/2 * (2*rand(sizeOut) - 1);
    W = - log( rand(sizeOut) );
    piover2 = pi/2;
    sclshftV =  piover2 + beta * V ;
    r = 1/piover2 * ( sclshftV .* tan(V) - ...
        beta * log( (piover2 * W .* cos(V) ) ./ sclshftV ) );
end

% Apply scale and shift parameters
if alpha ~= 1
    r = gamma * r + delta;
else
    r = gamma * r + (2/pi) * beta * gamma * log(gamma) + delta;
end

end
% <><><><><><><><><><><><> End of primary function <><><><><><><><><><><><>

% <><><><><><><><><><><><><> Secondary functions <><><><><><><><><><><><><>
% SUBROUTINE 1
function [err,commonSize,numElements] = genOutsize(nparams,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Function to find output size                                          %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

try
    tmp = 0;
    for argnum = 1:nparams
        tmp = tmp + varargin{argnum};
    end
    if nargin > nparams+1
        tmp = tmp + zeros(varargin{nparams+1:end});
    end
    err = 0;
    commonSize = size(tmp);
    numElements = numel(tmp);

catch
    err = 1;
    commonSize = [];
    numElements = 0;
end

end
% <><><><><><><><><><><> End of secondary function <><><><><><><><><><><><>

