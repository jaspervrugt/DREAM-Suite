function loglik = banana_lik(x)
% d-variate banana shaped target distribution (twisted normal)
%   Haario, H., E. Saksman, and J. Tamminen (1999), Adaptive proposal 
%       distribution for random walk Metropolis algorithm. Computational 
%       Statistics 14, 375–395, https://doi.org/10.1007/s001800050022

% Store local variables in memory
persistent b d invC log_Z

if isempty(invC)    % Local memory
    d = size(x,2);                                  % Target dimensionality
    b = 0.1;                                        % Target nonlinearity
    C = eye(d,d); C(1,1) = 100; invC = inv(C);      % Target covariance
    log_Z = log(((2*pi)^(-d/2)) * det(C)^(-1/2));   % Integration constant
    if ( d > 150 ), log_Z = 0; end                  % Log_Z to zero if -inf 
                                                    % --> inconsequential                
end

x(:,2) = x(:,2) + b * x(:,1).^2 - 100 * b;  % Banana-shaped nonlinearity 
                                            % between x(1) and x(2)

switch d
    case 1      % univariate normal
        loglik = log_Z - 1/2 * invC * (x.^2)';
    otherwise   % multivariate normal
        loglik = log_Z - 1/2 * sum(x'.*(invC*x'));
end

end