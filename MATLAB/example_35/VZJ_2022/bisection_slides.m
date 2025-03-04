function xr = bisection_slides(fun,xmin,xmax,nx)
% Root finding using the bracketing method 
%
% SYNOPSIS: xr = bisection_slides(fun,xmin,xmax)
%           xr = bisection_slides(fun,xmin,xmax,nx)
% INPUTS:   fun:    inline function of f(x) 
%           xmin:   minimum value of x
%           xmax:   maximum value of x
%           nx:     number of iterations        (default: 15)
% OUTPUTS:  xr:     estimate of root 
% 
if nargin < 4, nx = 15; end
a = xmin; b = xmax;             % lower and upper value of x
fa = fun(a); fb = fun(b);       % evaluate the function value at x=a and x=b
for i = 1:nx
    c = (b+a)/2; fc = fun(c);               % compute mid point of x=a and x=b, and fc = f(c)
    if ( sign(fa) ~= sign(fc) )             % compare sign function at x=a and x=c
        b = c; fb = fc;                     % not equal --> root in [a,c], thus b = c;
    else
        a = c; fa = fc;                     % equal --> root in [c,b]; thus; a = c;
    end
end
xr = (a+b)/2;                   % return the midpoint of last a and b as our root
[a xr b]; a-b