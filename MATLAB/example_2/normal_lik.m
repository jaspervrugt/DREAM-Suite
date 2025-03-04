function loglik = normal_lik(x)
% d-variate normal log-likelihood, N(µ,Σ) with µ = 0 and Σ = dxd covariance 
% matrix with 1:d on main diagonal and 0.5 correlation between dimensions 

comp_method = 2;    % 1: Exact using built-in functions [slow >10 year ago]
                    % 2: Some protection against underflow  
                    % 3: Some protection against underflow [faster, Note 1]
% What is underflow? L(x) → zero when x is far removed from target mean
%                         → log(L(x)) goes to -∞.
% Note 1: comp_method = 2: N = 1e5; d = 100; X = 20*rand(N,d) - 10;  
% MATLAB crashes, matrices become too large
% --> use comp_method 3 (= faster for large N and d)

persistent d mu C invC log_Z        % Local variables in memory
                                    % Careful: change of d causes problems

if isempty(invC)                    % Store invrnt target info local memory
    d = size(x,2);                  % # parameters
    mu = zeros(1,d);                % µ of zero: 1xd vector
    R = .5 * eye(d,d) ...           % dxd target correlation matrix R
            + .5 * ones(d,d); 
    C = nan(d,d);                   % Initialize dxd covariance matrix Σ
    for i = 1:d                     % Populate target covariance matrix
        for j = 1:d
            C(i,j) = R(i,j) * sqrt(i*j);
        end
    end
    invC = inv(C);                  % Inverse of covariance matrix
    log_Z = log(((2*pi)^(-d/2)) ... % Normalization constant of target
        * det(C)^(-1/2));
%   if (d > 150), log_Z = 0; end    % log_Z becomes -inf for large d 
                                    % --> set log_Z to zero for large d
end
x = x - mu;                         % Scale candidate point(s)

switch comp_method
    case 1  % Built-in function: no protection against underflow
        loglik = log(mvnpdf(x,zeros(1,d),C))';
    case 2  % Protection & faster  [calc_method = 1: slow > 10 year ago]
        loglik = log_Z ...          % (Nxd)*(dxd)*(dxN) = NxN
            -.5 *diag(x *invC*x')'; % Take diagonal elements                           
    case 3  % Protection & fastest & works for large N and large d
        if d == 1
            loglik = log_Z ...      % univariate Gaussian 
                - .5*invC*(x'.^2);  % 1x1 - 1x1 * 1xN = 1x1 or 1xN 
        else 
            loglik = log_Z ...      % multivariate Gaussian
                - .5 * sum( ...     % (Nxd)*((dxd)*(dxN))=(dxN).*(dxN)=dxN
                x'.*(invC*x'));     % sum(dxN) = 1xN = loglik N cand. pnts.           
        end
end

end
