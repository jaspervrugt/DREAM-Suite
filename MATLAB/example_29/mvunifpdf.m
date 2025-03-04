function pdf = mvunifpdf(x,a,b)
%MVUNIFPDF Computes the density of the multivariate uniform distribution 
% on the interval for a <= x <= b
% Written by Jasper A. Vrugt 
% To illustrate use of multivariate prior in example 29 of DREAM Package

if (any(x<a) || any(x>b))
    pdf = realmin; % Outside domain, pdf is zero but return realmin instead
else
    pdf = prod(unifpdf(x,a,b));
end

end