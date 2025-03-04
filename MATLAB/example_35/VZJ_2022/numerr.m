% Draw from prior ranges
close all;
N = 1000; d = 3;
eta_min = ones(1,3) * 1e-5; eta_max = [ 50 100 2 ];
[flag_Iexp,flag_t,flag_I,flag_Iexp,err_I,err_t,err_J1,err_J2] = deal(nan(N,1)); ss = zeros(N,1);
% Create initial parameter vectors
Eta = repmat(eta_min,N,1) + lhsdesign(N,d,'criterion','maximin','iteration',50) .* ...
    ( repmat(eta_max,N,1) - repmat(eta_min,N,1) );

case_to_run = 1;    % [1] numerr_results.mat: first case -> looks at numerical error t and I solutions
                    % [2] newton_iter.mat: second case --> looks at number of iterations as function of dt
switch case_to_run
    case 1
        % Specify the respective infiltration times
        t = [0:0.1:24]'; plugin.t = t; n = numel(plugin.t);
        % Run haverkamp_t for exact solution
        %plugin.I = [0:0.1:10]';
        %n = numel(plugin.I);
        % Define time
        
        % % tic
        % % for i = 1:N
        % %     [I,inf_rate_Iexp,flag_Iexp] = Haverkamp_I_exp(Eta(i,1:3),plugin);        % Solve infiltration - with constant rate assumption
        % %     plugin.I = I;
        % %     % Now lets use the (I,t) curve
        % %     [t,flag_t(i)] = Haverkamp_t(Eta(i,1:3),plugin);          % Solve t-form without constant rate
        % %     err_t(i,1) = sum(abs(t(1:end) - plugin.t(1:end)));
        % %     [I,inf_rate_I,flag_I(i)] = Haverkamp_I(Eta(i,1:3),plugin);        % Solve I-form without constant rate
        % %     err_I(i,1) = sum(abs(I(1:end) - plugin.I(1:end)));      % Compute infiltration error
        % %     % Compare Jacobian matrices - numerical and analytic
        % %     J_num = jac_Haverkamp_num(Eta(i,1:3),plugin,1);
        % %     J_anal = jac_Haverkamp_anal(Eta(i,1:3),plugin,1);
        % %     err_J1(i,1) = sum(sum(abs(J_num - J_anal)));
        % %     J_num = jac_Haverkamp_num(Eta(i,1:3),plugin,2);
        % %     J_anal = jac_Haverkamp_anal(Eta(i,1:3),plugin,2);
        % %     err_J2(i,1) = sum(sum(abs(J_num - J_anal)));
        % % end
        
        % Second loop - now with fixed infiltration to 20 cm
        [flag_t,flag_I,err_I,err_t] = deal(nan(N,1)); ss = zeros(N,1);
        plugin.I = [0:0.01:10]'; n = numel(plugin.I); t_time = 0; I_time = 0;
        for i = 1:N
            %     [I,inf_rate_Iexp,flag_Iexp] = Haverkamp_I_exp(Eta(i,1:3),plugin);        % Solve infiltration - with constant rate assumption
            %     plugin.I = I;
            % Now lets use the (I,t) curve
            ts = cputime;
            [t,flag_t(i)] = Haverkamp_t(Eta(i,1:3),plugin);          % Solve t-form without constant rate
            t_end = cputime; t_time = t_time + (t_end-ts);
            %    err_t(i,1) = sum(abs(t(1:end) - plugin.t(1:end)));
            if flag_t(i) == 2
                plugin.t = fix_time(t,n);
            else
                plugin.t = t;
            end
            ts = cputime;
            [I,inf_rate_I,flag_I(i)] = Haverkamp_I(Eta(i,1:3),plugin);        % Solve I-form without constant rate
            t_end = cputime; I_time = I_time + (t_end-ts);
            if flag_I(i) == 2
                [idx_ii,ss(i)] = fix_I(I,n,inf_rate_I,Eta(i,1:3));
            else
                idx_ii = n;
            end
            err_I(i,1) = sum(abs(I(1:idx_ii) - plugin.I(1:idx_ii)));      % Compute infiltration error
            % Jacobian matrices: infiltration form
            % %     J_num = jac_Haverkamp_num(Eta(i,1:3),plugin,1);
            % %     J_anal = jac_Haverkamp_anal(Eta(i,1:3),plugin,1);
            % %     err_J1(i,1) = sum(sum(abs(J_num(1:idx_ii,:) - J_anal(1:idx_ii,:))));
            % %     % Jacobian matrices: time form
            % %     J_num = jac_Haverkamp_num(Eta(i,1:3),plugin,2);
            % %     J_anal = jac_Haverkamp_anal(Eta(i,1:3),plugin,2);
            % %     err_J2(i,1) = sum(sum(abs(J_num(1:idx_ii,:) - J_anal(1:idx_ii,:))));
        end
        [t_time I_time]
        
        idx_t1 = find(flag_t == 1); idx_t2 = find(flag_t == 2); not_converged_t = numel(idx_t2);
        % How many solutions did not converge according to I
        idx_I1 = find(flag_I == 1); idx_I2 = find(flag_I == 2); not_converged_I = numel(idx_I2);
        [not_converged_t not_converged_I]
        % How many solutions were at inf_rate == Ks?
        sum(ss)
        subplot(1,3,1),plot(Eta(idx_I2,1),Eta(idx_I2,2),'ro'); xlabel('S'); ylabel('Ks'); axis square
        subplot(1,3,2),plot(Eta(idx_I2,1),Eta(idx_I2,3),'ro'); xlabel('S'); ylabel('\beta'); axis square
        subplot(1,3,3),plot(Eta(idx_I2,2),Eta(idx_I2,3),'ro'); xlabel('Ks'); ylabel('\beta'); axis square
        % Compute ratio of eta1 and eta2
        xi_err = Eta(idx_I2,2)./Eta(idx_I2,1).^2;
        gamma = 2*Eta(:,3).*Eta(:,2)./Eta(:,1).^2;
        gamma(idx_I2) > 709.78/10
        figure(2),hist(2*Eta(:,3).*Eta(:,2)./Eta(:,1).^2,20);
        figure(3)
        for z = 1:numel(idx_t2)
            t = Haverkamp_t(Eta(idx_t2(z),1:3),plugin);
            plot(t,plugin.I,'r'); if z == 1; hold on; end
        end
        % SAVE RESULTS
        save numerr_results.mat Eta err_I idx_I1 idx_I2 idx_t1 idx_t2 flag_t flag_I t_time I_time not_converged_I not_converged_t
        
    case 2 % Check effect of dt on number of iterations Newton's method
        %dI = [0.01 0.02 0.05 0.1 0.2 0.5 1.0 2.0]; M = numel(dI);
        dt = [0.01 0.02 0.04 0.05 0.10 0.2 0.4 0.5 1.0 2.0 4 5]; M = numel(dt);
        % Now run the same samples but check the number of iterations
        sumK = zeros(M,3); err_Imethod = nan(N,3,M);
        for u = 1:M
            u
            plugin.t = [0:dt(u):20]'; n = numel(plugin.t); T = N*ones(1,3);
            % Iterate over MC samples
            for i = 1:N
                for method = 1:3
                    [I,inf_rate,flag,K] = Haverkamp_IK(Eta(i,1:3),plugin,1e-12,20,method);
                    if flag == 2
                        %sumK(u,method,u) = nan;
                        err_Imethod(i,method,u) = nan; T(method) = T(method)-1; 
                    else
                        if any(isnan(K))
                            inan = find(isnan(K)); inan = inan(1)-1;
                            [i inan method]
                        else
                            inan = n;
                        end
                        sumK(u,method) = sumK(u,method) + sum(K(1:inan));
                        % cannot compute error - forward modeling only
                        err_Imethod(i,method,u) = I(inan);
                        %                        err_Imethod(i,method,u) = sum(abs(I(1:inan) - plugin.I(1:inan)));
                    end
                end
            end
            % Normalize sumK with number of samples
            sumK(u,1:3) = sumK(u,1:3)./(T*n)            
        end
        % Print sumK
        sumK
        % How many prematurely converged runs for each dt?
        for u = 1:M
            Not_A_Num(u,1:3) = sum(isnan(err_Imethod(1:N,1:3,u)));
        end
        
        % Save results for plotting in fig_paper2
        save iter_Newton.mat dt sumK err_Imethod
        % NOTE: WITH FIXED TIMES: MANY MORE RUNS OF METHOD 1 (I_0 = I_low)
        % do not converge than those of method 2 (midpoint) and 3 (upper)
        % What do we do with this?
end

% Now test the solution of Haverkamp_I2