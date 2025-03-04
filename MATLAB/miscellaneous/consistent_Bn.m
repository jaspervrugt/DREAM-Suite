function [Beta_n] = consistent_Bn(Jn,q)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns a consistent estimator of the variability matrix
% according to Newey and West (1987)
%
% SYNOPSIS: [BBn,Bn,Jn] = consistent_Bn(Jn,q)
%   INPUT
%       Jn      nxp Jacobian matrix of the p parameters
%       q       scalar with maximum lag of correlation/dependence of
%               entries of the score function
%   OUTPUT
%       Beta_n  pxp variability matrix according to Newey and West (1987)
%
% Written by Jasper A. Vrugt
% April 19, 2024
% Laguna Niguel, United States
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How many samples, n, and parameters, d
[n,d] = size(Jn);
% Initialize covariance terms at lags 1,2,...,q
Beta_j = nan(d,d,q);
% Compute mean variability matrix (= Bn in sandwich form)
Bn = 1/n * Jn(1:n,1:d)' * Jn(1:n,1:d);
% Now compute covariance terms at lags 1,...,q (Newey and West, 1987)
for j = 1:q
    t = j + 1;
    % Compute mean variability matrix (= Beta_j/n in our paper)
    Beta_j(:,:,j) = 1/(n-j) * Jn(t:n,1:d)' * Jn(1:n-j,1:d);
end
% For first lag, j = 1 we need to compute the first-lag correction between
% the score function at t and t-1 (= t-j).
% We can solve this in a loop as follows
% a = zeros(d,d); for t = 2:n, a = a + Jn(t,1:d)'*Jn(t-1,1:d); end; a = a/n;
% But you can also do this at once
% a2 = Jn(2:n,1:d)'*Jn(1:n-1,1:d)/n
% a is equal to a2

% Now compute correction
dB = zeros(d,d);
for j = 1:q
    dB = dB + (1 - j/(q + 1)) * (Beta_j(:,:,j) + Beta_j(:,:,j)');
end
% Compute average variability matrix (= Beta_n in our paper)
Beta_n = Bn + dB;

end