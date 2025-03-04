function p = stblpdf(x,alpha,beta,gam,delta,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% STBLPDF(x,alpha,beta,gamma,delta) returns pdf of stable distribution  %%
%%    with characteristic exponent, alpha, skewness, beta, scale         %%
%%    parameter, gam, and location parameter, delta, at the values in X. %%
%%    We use the parameterization of stable distribtuions used in [2].   %%
%%    The characteristic function phi(t) of a S(alpha,beta,gam,delta)    %%
%%    random variable has the form                                       %%
%%                                                                       %%
%% phi(t) = exp(-gam^alpha |t|^alpha [1-i beta (tan(pi alpha/2) sign(t)] %%
%%                 + i delta t )  if alpha ~= 1                          %%
%%                                                                       %%
%% phi(t) = exp(-gam |t| [1+i beta (2/pi) (sign(t)) log|t|] + i delta t  %%
%%                                if alpha = 1                           %%
%%                                                                       %%
%% The size of P is the size of X. alpha, beta, gam and delta must be    %%
%% scalars                                                               %%
%%                                                                       %%
%% P = STBLPDF(X,alpha,beta,gam,delta,tol) computes the pdf to within an %%
%%     absolute error of tol                                             %%
%%                                                                       %%
%% The algorithm works by computing the numerical integrals in Theorem   %%
%% 1 in [1] using MATLAB's integral (not quadv) function. The integrands %%
%% are smooth non-negative functions, but for certain parameter values   %%
%% can have sharp peaks which might be missed.  To avoid this, STBLPDF   %%
%% locates the maximum of this integrand and breaks the integral into    %%
%% two pieces centered around this maximum (this is exactly the idea     %%
%% suggested in [1] )                                                    %%
%%                                                                       %%
%% If abs(alpha - 1) < 1e-5,  alpha is rounded to 1                      %%
%%                                                                       %%
%% P = STBLPDF(...,'quick') skips the step of locating the peak in the   %%
%% integrand, and thus is faster, but is less accurate deep into tails   %%
%% of the pdf. This option is useful for plotting. In place of 'quick',  %%
%% stblpdf also excepts a logical true or false (for quick or not quick) %%
%%                                                                       %%
%% See also: STBLRND, STBLCDF, STBLINV, STBLFIT                          %%
%%                                                                       %%
%% References:                                                           %%
%% [1] Nolan, J.P. (1997), Numerical Calculation of Stable Densities and %%
%%     Distribution Functions, Commun. Statist. - Stochastic Modles,     %%
%%     13(4), 759-774                                                    %%
%% [2] Samorodnitsky, G., and M.S. Taqqu (1994), Stable non-Gaussian     %%
%%     random processes: stochastic models with infinite variance,       %%
%%     CRC Press                                                         %%
%%                                                                       %%
%% Modified by Jasper A. Vrugt, April 2017                               %%
%% University of California Irvine                                       %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 5
    error(['DREAM_PACKAGE:stblpdf:TooFewInputs','Requires at least ',...
        'five input arguments.']);
end

% Check parameters
if alpha <= 0 || alpha > 2 || ~isscalar(alpha)
    error(['DREAM_PACKAGE:stblpdf:BadInputs',' "alpha" must be a ',...
        'scalar which lies in the interval (0,2]']);
end
if abs(beta) > 1 || ~isscalar(beta)
    error(['DREAM_PACKAGE:stblpdf:BadInputs',' "beta" must be a scalar',...
        'which lies in the interval [-1,1]']);
end
if gam < 0 || ~isscalar(gam)
    error(['DREAM_PACKAGE:stblpdf:BadInputs',' "gam" must be a ',...
        'non-negative scalar']);
end
if ~isscalar(delta)
    error('DREAM_PACKAGE:stblpdf:BadInputs',' "delta" must be a scalar');
end

% Warn if alpha is very close to 1 or 0
if ( 1e-5 < abs(1 - alpha) && abs(1 - alpha) < .02) || alpha < .02
    warning('DREAM_PACKAGE:stblpdf:ScaryAlpha',...
        'Difficult to approximate pdf for alpha close to 0 or 1')
end

calc_method = 2;
% warnings will happen during call to QUADV [= integral], and it's okay
warning('off');

% Check and initialize additional inputs
quick = false;
tol = [];
for i=1:length(varargin)
    if strcmp(varargin{i},'quick')
        quick = true;
    elseif islogical(varargin{i})
        quick = varargin{end};
    elseif isscalar(varargin{i})
        tol = varargin{i};
    end
end

if isempty(tol)
    if quick
        tol = 1e-8;
    else
        tol = 1e-12;
    end
end


% <><><><><><><><><><><><><><> Compute pdf <><><><><><><><><><><><><><><><>

% Check if simple case, if so be quick, if not do the general algorithm
if alpha == 2                           % Gaussian distribution
    x = (x - delta)/gam;                % Standardize
% %     p = 1/sqrt(4*pi) * exp(-.25*x.^2); % ~ N(0,2)
    p = 1/sqrt(4*pi) * exp(-1/4*x.^2);  % ~ N(0,2)
    p = p/gam;                          % rescale
elseif alpha==1 && beta == 0            % Cauchy distribution
    x = (x - delta)/gam;                % Standardize
    p = (1/pi) * 1./(1 + x.^2);
    p = p/gam;                          % rescale
elseif alpha == .5 && abs(beta) == 1    % Levy distribution
    x = (x - delta)/gam;                % Standardize
    p = zeros(size(x));
    if  beta ==1
        p(x <= 0) = 0;
        p(x > 0) = sqrt(1/(2*pi)) * exp(-.5./x(x>0)) ./ x(x>0).^1.5;
    else
        p(x >= 0) = 0;
        p(x < 0) = sqrt(1/(2*pi)) * exp(.5./x(x<0)) ./ (-x(x<0)).^1.5;
    end
    p = p/gam;                      % rescale
elseif abs(alpha - 1) > 1e-5        % Gen. Case, alpha ~= 1
    xold = x; % Save for later
    % Standardize in (M) parameterization ( See equation (2) in [1] )
    x = (x - delta)/gam - beta * tan(alpha*pi/2);
    % Compute pdf
    p = zeros(size(x));
    zeta = - beta * tan(pi*alpha/2);
    theta0 = (1/alpha) * atan(beta*tan(pi*alpha/2));
    A1 = alpha*theta0;
    A2 = cos(A1)^(1/(alpha-1));
    exp1 = alpha/(alpha-1);
    alpham1 = alpha - 1;
    c2 = alpha ./ (pi * abs(alpha - 1) * ( x(x>zeta) - zeta) );
    V = @(theta) A2*(cos(theta)./sin(alpha*(theta + theta0))).^exp1 .* ...
        cos( A1 + alpham1*theta ) ./ cos(theta);
    % x > zeta, calculate integral using QUADV [= integral]
    if any(x(:) > zeta)
        xshift = (x(x>zeta) - zeta) .^ exp1;
        if beta == -1 && alpha < 1
            p(x > zeta) = 0;
        elseif ~quick % Locate peak in integrand and split up integral
            g = @(theta) xshift(:) .* V(theta) - 1;
            R = repmat([-theta0, pi/2 ],numel(xshift),1);
            if abs(beta) < 1
                theta2 = bisectionSolver(g,R,alpha);
            else
                theta2 = bisectionSolver(g,R,alpha,beta,xshift);
            end
            theta2 = reshape(theta2,size(xshift));
            % change variables so the two integrals go from
            % 0 to 1/2 and 1/2 to 1.
            theta2shift1 = 2*(theta2 + theta0);
            theta2shift2 = 2*(pi/2 - theta2);
            g1 = @(theta) xshift .* ...
                V(theta2shift1 * theta - theta0);
            g2 = @(theta) xshift .* ...
                V(theta2shift2 * (theta - .5) + theta2);
            zexpz = @(z) max(0,z .* exp(-z)); % use max incase of NaN
            switch calc_method
                case 1      % old implementation with quadv function
                    p(x > zeta) = c2 .* ...
                        (theta2shift1 .* quadv(@(theta) ...
                        zexpz(g1(theta)),0,0.5,tol) + theta2shift2 .* ...
                        quadv(@(theta) zexpz(g2(theta)),0.5,1,tol));
                otherwise   % new approach using built-in integral function
                    p(x > zeta) = c2 .* ...
                        (theta2shift1 .* integral(@(theta) ...
                        zexpz(g1(theta)),0,0.5,'AbsTol',tol, ...
                        'ArrayValued',1) + theta2shift2 .* ...
                        integral(@(theta) zexpz(g2(theta)), ...
                        0.5,1,'AbsTol',tol,'ArrayValued',1));
            end
        else  % be quick - calculate integral without locating peak
            % Use a default tolerance of 1e-6
            g = @(theta) xshift * V(theta);
            zexpz = @(z) max(0,z .* exp(-z)); % use max incase of NaN
            switch calc_method
                case 1      % old implementation with quadv function
                    p(x > zeta) = c2 .* quadv(@(theta) zexpz( g(theta) ), ...
                        -theta0 , pi/2, tol );
                otherwise   % new approach using built-in integral function
                    p(x > zeta) = c2 .* integral(@(theta) zexpz(g(theta)), ...
                        -theta0,pi/2,'AbsTol',tol,'ArrayValued',1);
            end
        end
        p(x > zeta) = p(x>zeta)/gam; % rescale
    end
    % x = zeta, this is easy
    if any( abs(x(:) - zeta) < 1e-8 )
        p(abs(x - zeta) < 1e-8) = max(0,gamma(1 + 1/alpha)*...
            cos(theta0)/(pi*(1 + zeta^2)^(1/(2*alpha))));
        p(abs(x - zeta) < 1e-8) = p( abs(x - zeta) < 1e-8 )/gam; % rescale

    end
    % x < zeta, recall function with -xold, -beta, -delta
    % This doesn't need to be rescaled.
    if any(x(:) < zeta)
        p(x < zeta) = stblpdf(-xold(x<zeta),alpha,-beta,...
            gam , -delta , tol , quick);
    end
else                    % Gen case, alpha = 1
    x = (x - (2/pi) * beta * gam * log(gam) - delta)/gam; % Standardize
    % Compute pdf
    piover2 = pi/2;
    twooverpi = 2/pi;
    oneoverb = 1/beta;
    theta0 = piover2;
    % Use logs to avoid overflow/underflow
    logV = @(theta) log(twooverpi * ((piover2 + beta *theta)./ ...
        cos(theta))) + (oneoverb * (piover2 + beta *theta) .* tan(theta));
    c2 = 1/(2*abs(beta));
    xterm = ( -pi*x/(2*beta));
    if ~quick   % Locate peak in integrand and split up integral
        % Use a default tolerance of 1e-12
        logg = @(theta) xterm(:) + logV(theta) ;
        R = repmat([-theta0, pi/2 ],numel(xterm),1);
        theta2 = bisectionSolver(logg,R,1-beta);
        theta2 = reshape(theta2,size(xterm));
        % change variables so the two integrals go from
        % 0 to 1/2 and 1/2 to 1.
        theta2shift1 = 2*(theta2 + theta0);
        theta2shift2 = 2*(pi/2 - theta2);
        logg1 = @(theta)  xterm + ...
            logV(theta2shift1 * theta - theta0);
        logg2 = @(theta)  xterm + ...
            logV(theta2shift2 * (theta - .5) + theta2);
        zexpz = @(z) max(0,exp(z) .* exp(-exp(z))); % use max incase of NaN
        switch calc_method
            case 1      % old implementation with quadv function
                p = c2 .* (theta2shift1 .* quadv(@(theta) ...
                    zexpz(logg1(theta)),0,0.5,tol) + theta2shift2 .* ...
                    quadv(@(theta) zexpz(logg2(theta)),0.5,1,tol));
            otherwise   % new approach using built-in integral function
                p = c2 .* ...
                    (theta2shift1 .* integral(@(theta) ...
                    zexpz(logg1(theta)),0,0.5,'AbsTol',tol, ...
                    'ArrayValued',1) + theta2shift2 .* ...
                    integral(@(theta) zexpz(logg2(theta)), ...
                    0.5,1,'AbsTol',tol,'ArrayValued',1));
        end
    else    % be quick - calculate integral without locating peak
        % Use a default tolerance of 1e-6
        logg = @(theta) xterm + logV(theta);
        zexpz = @(z) max(0,exp(z) .* exp(-exp(z))); % use max incase of NaN
        switch calc_method
            case 1      % old implementation with quadv function
               p = c2 .* quadv(@(theta) zexpz(logg(theta)), ...
                   -theta0,pi/2,tol);
            otherwise   % new approach using built-in integral function
               p = c2 .* integral(@(theta) zexpz(logg(theta)), ...
                   -theta0,pi/2,'AbsTol',tol,'ArrayValued',1);
        end
    end
    p = p/gam; % rescale
end
p = real(p); % just in case a small imaginary piece crept in
% This might happen when (x - zeta) is really small

end
% <><><><><><><><><><><><> End of primary function <><><><><><><><><><><><>

% <><><><><><><><><><><><><> Secondary functions <><><><><><><><><><><><><>
% SUBROUTINE 1
function X = bisectionSolver(f,R,alpha,varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Solves equation g(theta)-1 = 0 in stblpdf using vectorized bisection  %%
%% and a tolerance of 1e-5. The solution to this equation is used to     %%
%% increase accuracy in the calculation of a numerical integral          %%
%%                                                                       %%
%% If alpha ~= 1 and 0 <= beta < 1, the equation always has a solution   %%
%%                                                                       %%
%% If alpha > 1 and beta <= 1, then g is monotone decreasing             %%
%%                                                                       %%
%% If alpha < 1 and beta < 1, then g is monotone increasing              %%
%%                                                                       %%
%% If alpha = 1,  g is monotone increasing if beta > 0 and monotone      %%
%% decreasing is beta < 0. Input alpha = 1 - beta to get desired results %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 2
    error(['DREAM_PACKAGE:bisectionSolver:TooFewInputs',...
        'Requires at least two input arguments.']);
end
noSolution = false(size(R,1));
% if ~isempty(varargin)
%     beta = varargin{1};
%     xshift = varargin{2};
%     if abs(beta) == 1
%         V0=(1/alpha)^(alpha/(alpha-1))*(1-alpha)*cos(alpha*pi/2)*xshift;
%         if alpha > 1
%             noSolution = V0 - 1 %>= 0;
%         elseif alpha < 1
%             noSolution = V0 - 1 %<= 0;
%         end
%     end
% end

tol = 1e-6; maxiter = 30;
[N,M] = size(R);
if M ~= 2
    error('DREAM_PACKAGE:bisectionSolver:BadInput',...
        '"R" must have 2 columns');
end
a = R(:,1);
b = R(:,2);
X = (a+b)/2;
try
    val = f(X);
catch ME
    error(['DREAM_PACKAGE:bisectionSolver:BadInput',...
        'Input function inconsistint with rectangle dimension'])
end
if size(val,1) ~= N
    error(['DREAM_PACKAGE:bisectionSolver:BadInput',...
        'Output of function must be a column vector with ',...
        'dimension of input']);
end
% Main loop
val = inf; iter = 0;
while( max(abs(val)) > tol && iter < maxiter )
    X = (a + b)/2;
    val = f(X);
    l = (val > 0);
    if alpha > 1
        l = 1-l;
    end
    a = a.*l + X.*(1-l);
    b = X.*l + b.*(1-l);
    iter = iter + 1;
end
if any(noSolution(:))
    X(noSolution) = (R(1,1) + R(1,2))/2;
end

end
% <><><><><><><><><><><> End of secondary function <><><><><><><><><><><><>
