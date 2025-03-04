function [eta_opt,SSR_opt] = simplex_hav(I_meas,t_meas,approach,N)
% Implements simplex method
plugin.I = I_meas; plugin.t = t_meas;
switch approach
    case 1
        F = @(eta) sum ((I_meas - Haverkamp_I(eta,plugin)).^2);
    case 2
        F = @(eta) sum ((t_meas - Haverkamp_t(eta,plugin)).^2);
end

% % % Define bounds on parameters
% % eta_min = [ 0 0 0 ];
% % eta_max = [ 50 50 2 ];
% % % how many parameters? 21.9477   17.1336    0.0000
% % d = 3;
% % % Initial value?
% % Eta = [20*rand(N,1) 50*rand(N,1) 2*rand(N,1)]';
% % eta_opt = nan(N,d); SSR_opt = nan(N,1);
% % % Loop over different eta values
% % for trial = 1:N
% %     eta_opt(trial,1:d) = fminsearch(@(eta) F(eta),Eta(1:d,trial));
% %     SSR_opt(trial,1) = F(eta_opt(trial,1:d));
% % end


eta_min = [ 0 0 0 ];
eta_max = [ 50 50 2 ];
% Simplex_method
p = 3; % number of parameters of model ( = Rosenbrock function )
% Number of vertices simplex ( = p + 1 );
m = p + 1;
% set algorithmic paraeter values
alpha = 1; gamma = 2; rho = 1/2; sigma = 1/2; iter = 1;
% define initial simplex
eta = repmat(eta_min,m,1) + rand(m,p) .* (repmat(eta_max,m,1) - repmat(eta_min,m,1) ); 
% Now transpose 
eta = eta';
% % eta(1,:) = -1.5 + 3 * rand(1,m); 
% % eta(2,:) = -1 + 3 * rand(1,m);
% Write frame

% calculate F of each of m = p + 1 points
F_x = nan(1,m);
for j = 1 : m
    F_x(j) = F(eta(:,j));
end % Note -> can vectorize this as well - but for clarity
[na,idx] = sort(F_x);
% Initialize termination criterion
notconverged = 1; maxiter = 200000;

% while statement
while ( notconverged == 1 ) && ( iter < maxiter )
    % sort the function values from low to high
    [na,idx] = sort(F_x); F_x = F_x(idx); eta = eta(1:p,idx);
    % Calculate centroid of simplex (without worst point)
    x0 = mean(eta(1:p,1:m-1),2);
    % Now determine maximum value of chi
    lb = eta_min'; ub = eta_max';
    % determine jump
    jump = ( x0 - eta(1:p,m) );
    for j = 1 : p
        if jump(j) < 0
            chi_min(j,1) = (ub(j) - x0(j))/jump(j);
            chi_max(j,1) = (lb(j) - x0(j))/jump(j);
        elseif jump(j) == 0
            chi_min(j,1) = -inf; chi_max(j,1) = inf;
        else
            chi_min(j,1) = (lb(j) - x0(j))/jump(j);
            chi_max(j,1) = (ub(j) - x0(j))/jump(j);
        end
    end
    
    % REFLECTION STEP AND CALCULATE F(R)
    R = x0 + alpha * ( x0 - eta(1:p,m) ); F_R = F(R);
    % NOW CHECK
    if ( F_x(1) < F_R ) && ( F_R < F_x(m-1) )
        % Replace worst point with R and store function value
        eta(1:p,m) = R; F_x(m) = F_R; step = 'reflection';
    elseif F_R < F_x(1)
        % EXPANSION STEP AND CALCULATE F(E)
        E = x0 + gamma * ( x0 - eta(1:p,m) ); F_E = F(E);
        % Now check
        if F_E < F_R
            eta(1:p,m) = E; F_x(m) = F_E; step = 'expansion';
        else
            eta(1:p,m) = R; F_x(m) = F_R; step = 'reflection';
        end
    elseif F_R >= F_x(m-1) % (now clear that F_R > F(m-1) )
        % CONTRACTION STEP
        C = x0 + rho * ( x0 - eta(1:p,m) ); F_C = F(C);
        if F_C < F_x(m)
            eta(1:p,m) = C; F_x(m) = F_C; step = 'contraction';
        else
            % REDUCTION STEP
            for j = 2:m
                eta(1:p,j) = eta(1:p,1) + sigma * (eta(1:p,j) - eta(1:p,1)); F_x(j) = F(eta(1:p,j));
            end
            step = 'reduction';
        end
    end
    % update iter
    iter = iter + 1;

    % calculate difference between theta_old and theta
    notconverged = max(max(abs(bsxfun(@minus,eta,mean(eta,2))))) > 1e-6;
    
end

[~,idx] = sort(F_x); F_x = F_x(idx); eta = eta(1:p,idx);
SSR_opt = F_x;
eta_opt = eta;
