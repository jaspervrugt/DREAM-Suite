function [eta_opt,SSR_opt,K] = least_squares_haverkamp_der(t_meas, ...
    I_meas,dIdt,dtdI,approach,eta_min,eta_max,nu,N,eta_fun,k_max)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Least squares estimation of parameters of Haverkamp infiltration equation. This    %%
%% function differs from least_squares_haverkamp in that 
%%  SYNOPSIS: I = least_squares_haverkamp_der(t_meas,I_meas,dIdt,dtdI,approach, ...   %%
%%                    eta_min,eta_max)                                                %%
%%            I = least_squares_haverkamp_der(t_meas,I_meas,dIdt,dtdI,approach, ...   %%
%%                    eta_min,eta_max,nu)                                             %%
%%            I = least_squares_haverkamp_der(t_meas,I_meas,dIdt,dtdI,approach, ...   %%
%%                    eta_min,eta_max,nu,N)                                           %%
%%            I = least_squares_haverkamp_der(t_meas,I_meas,dIdt,dtdI,approach, ...   %%
%%                    eta_min,eta_max,nu,N,eta_fun)                                   %%
%%            I = least_squares_haverkamp_der(t_meas,I_meas,dIdt,dtdI,approach, ...   %%
%%                    eta_min,eta_max,nu,N,eta_fun,k_max)                             %%
%% where                                                                              %%
%%  eta  [input]       4x1 vector with S [cm/h^0.5], Ks [cm/h], beta [-], Ki [cm/h]   %%
%%  t    [input]       nx1 vector with time, t, in hours (h)                          %%


%%  eta  [input]       4x1 vector with S [cm/h^0.5], Ks [cm/h], beta [-], Ki [cm/h]   %%
%%  t    [input]       nx1 vector with time, t, in hours (h)                          %%
%%  rtol [opt. input]  tolerance on function value at root        (default: 1e-12)    %%
%%  kmax [opt. input]  maximum number of Newton iterations        (default: 20)       %%
%%  I    [output]      cumulative infiltration, I (cm), as function of time, t (h)    %%
%%  i    [output]      infiltration rate, i (cm/h), as function of time, t (h)        %%
%%  flag [output]      exit flag: [1] exact [2] approximate                           %%
%%                                                                                    %%
%%  LITERATURE                                                                        %%
%%  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and               %%
%%      H. Vereecken (2023), The time validity of Philip%s two-term infiltration      %%
%%      equation: An elusive theoretical quantity? Vadose Zone Journal, e20309,       %%
%%      pp. 1-25, https://doi.org/10.1002/vzj2.20309                                  %%
%%  Vrugt, J.A. and y. Gao (2022), On the three-parameter infiltration equation of    %%
%%      Parlange et al. (1982): Numerical solution, experimental design, and          %%
%%      parameter estimation, Vadose Zone Journal, 21:e20167, pp. 1-25,               %%
%%      https://doi.org/10.1002/vzj2.20167                                            %%
%%                                                                                    %%
%%  © Written by Jasper A. Vrugt, Jan. 2019                                           %%
%%  University of California Irvine                                                   %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 11, k_max = 250; end
if nargin < 10, eta_fun = 1e-6; end
if nargin < 9, N = 50; end
if nargin < 8, nu = 2; end

% Approach one or two
switch approach
    case 1 % Infiltration form
        plugin.t = t_meas; f = @(eta) Haverkamp_I_der(eta,plugin); Y_meas = dIdt; 
    case 2 % Time form
        plugin.I = I_meas; f = @(eta) Haverkamp_t_der(eta,plugin); Y_meas = dtdI;
end

% how many parameters?
d = 3;
% Determine method
%method = 'gauss-newton';
method = 'levenberg-marquardt';
% warning off % for badly scaled d_eta matrix
%Par_info.mu = [ 10 10 1 ];
%Par_info.cov = [ 5 0 0; 0 5 0; 0 0 0.3 ];
%Eta = mvnrnd(Par_info.mu,Par_info.cov,N)';
eta_mini = [ 1e-3 1e-3 1e-3 ]'; eta_maxi = [ 25 25 2 ]';
Eta = repmat(eta_mini,1,N) + lhsdesign(d,N) .* ( repmat(eta_maxi,1,N) - repmat(eta_mini,1,N) );
% Initialize eta_opt, K and SSR_opt
eta_opt = nan(N,d); [K,SSR_opt] = deal(nan(N,1));
% Loop over different eta values
for trial = 1:N
    eta = Eta(1:d,trial); 
    % Now compute nxp Jacobian matrix at current iterate
    Jf = jac_Haverkamp_der_anal(eta,plugin,approach); Jf = Jf(2:end,:);
    % Compute initial lambda
    lambda = max(diag(Jf'*Jf)); 
    % Initialize shift vector
    d_eta = inf * ones(d,1);
    % Now initialize loop counter
    k = 2;
    % Now while loop
    while (max(abs(d_eta)) > eta_fun) && (k <= k_max)
        % Compute nx1 residual vector
        e = Y_meas - f(eta(1:d,k-1));
        % Compute initial SSR both methods (only used by LM but for check of GN)
        if k == 2
            OF_eta = sum(e.^2); OF_start = OF_eta;
        end       
        % Now compute nxp Jacobian matrix at current iterate
        Jf = jac_Haverkamp_der_anal(eta(1:d,k-1),plugin,approach); Jf = Jf(2:end,:);
        % Compute update vector
        switch method
            case 'gauss-newton'
                % Add small perturbation to avoid badly scaled matrix
                d_eta = (Jf'*Jf + eps*eye(d,d))\(Jf'*e);
            case 'levenberg-marquardt'
                % Compute dxd dampling matrix
                D = lambda * diag(diag(Jf'*Jf));
                % Compute shift vector
                d_eta = (Jf'*Jf + D)\(Jf'*e);
        end
        % Now check whether we stay in bound
        [var_sigma_min,var_sigma_max] = deal(nan(d,1));
        for j = 1:d
            switch sign(d_eta(j))
                case -1
                    var_sigma_min(j) = (eta_max(j) - eta(j,k-1)) ./ d_eta(j);
                case 0
                    var_sigma_min(j) = -inf;
                case 1
                    var_sigma_min(j) = (eta_min(j) - eta(j,k-1)) ./ d_eta(j);
            end
        end
        for j = 1:d
            switch sign(d_eta(j))
                case -1
                    var_sigma_max(j) = (eta_min(j) - eta(j,k-1)) ./ d_eta(j);
                case 0
                    var_sigma_max(j) = inf;
                case 1
                    var_sigma_max(j) = (eta_max(j) - eta(j,k-1)) ./ d_eta(j);
            end
        end
        var_sigma = min(var_sigma_max); %; var_sigma, pause
        if var_sigma > 1
            var_sigma = 1; %fprintf('step length set to 1\n');
            % do nothing
        else
            % no boundary handling here
       %     var_sigma = 0.9*var_sigma;
        end
% %       var_sigma = 1;
        % Determine update
        switch method
            case 'gauss-newton'
                % Update parameter values
                eta(1:d,k) = eta(1:d,k-1) + var_sigma * d_eta;
            case 'levenberg-marquardt'
                % Update parameter values
                eta_try = eta(1:d,k-1) + d_eta; 
                % if out of bound? 
                % NEW IMPLEMENTATION - FROM LATEX PAPER
                for u = 1:d
                    if eta_try(u) < eta_min(u)
%                        eta_try(u) = eta(u,k-1) + 0.95 * var_sigma * d_eta(u);
                        eta_try(u) = eta_min(u) + 0.05*(eta(u,k-1) - eta_min(u));
                    end
                    if eta_try(u) > eta_max(u)
%                        eta_try(u) = eta(u,k-1) + 0.95 * var_sigma * d_eta(u);
                        eta_try(u) = eta_max(u) - 0.05*(eta_max(u) - eta(u,k-1));
                    end
                end
                % Compute nx1 residual vector of eta_try
                e_try = Y_meas - f(eta_try);
                % Compute objective function of eta_try
                OF_try = sum(e_try.^2);
                % Now check
                if OF_try < OF_eta
                    % Accept proposed eta and decrease lambda
                    eta(1:d,k) = eta_try; lambda = lambda/nu; OF_eta = OF_try;
                else
                    % Keep current value and increase lambda
                    eta(1:d,k) = eta(1:d,k-1); lambda = lambda * nu;
                end
        end
        % update counter
        k = k + 1;
    end
    %[OF_start OF_eta]
    % Store
    K(trial,1) = k-1;
    eta_opt(trial,1:d) = eta(1:d,k-1);
    switch method
        case 'gauss-newton'
            e = I_meas - f(eta_opt(trial,1:d)); OF_eta = sum(e.^2);
    end
    SSR_opt(trial,1) = OF_eta;
end