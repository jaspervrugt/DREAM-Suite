function [X_opt,F_opt,K_opt] = least_squares_Haverkamp(d,t_meas, ...
    I_meas,app,nu,N,xtol,maxiter)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Least squares estimation of parameters of Haverkamp infiltration      %%
%% equation using the Levenberg-Marquardt algorithm                      %%
%%                                                                       %%
%%  SYNOPSIS: X = least_squares_Haverkamp(d,t_meas,I_meas)               %%
%%            X = least_squares_Haverkamp(d,t_meas,I_meas,app)           %%
%%            X = least_squares_Haverkamp(d,t_meas,I_meas,app,nu)        %%
%%            X = least_squares_Haverkamp(d,t_meas,I_meas,app,nu,N)      %%
%%            X = least_squares_Haverkamp(d,t_meas,I_meas,app,nu,N,xtol) %%
%%            X = least_squares_Haverkamp(d,t_meas,I_meas,app,nu,N, ...  %%
%%                    xtol,maxiter)                                      %%
%% where                                                                 %%
%%  d           [input] # parameters Haverkamp infiltration equation     %%
%%   = 3            3-parameter form                                     %%
%%                      S  [cm/h^0.5]                                    %% 
%%                      Ks [cm/h]                                        %%
%%                      ß  [-]                                           %%
%%   = 4            4-parameter form                                     %%
%%                      S  [cm/h^0.5]                                    %% 
%%                      Ks [cm/h]                                        %%
%%                      ß  [-]                                           %%
%%                      Ki [cm/h]                                        %%
%%  t_meas      [input] nx1 vector with time, t, in hours (h)            %%
%%  I_meas      [input] nx1 vector with cum. inf. I, in centimeters (cm) %%
%%  approach    [input] OPT: Residuals to minimimze         (def: 1)     %%
%%   = 1            Minimize cumulative inf. residuals      default      %%
%%   = 2            Minimize time residuals                              %%
%%  nu          [input] OPT: Multiplicative parameter LM    (def: 2)     %%
%%  N           [input] OPT: # trials (multi-start)         (def: 50)    %%
%%  xtol        [input] OPT: Convergence threshold prmters  (def: 1e-6)  %%
%%  maxiter     [input] OPT: Maximum number iterations      (def: 250)   %%
%%  X_opt       [outpt] Least squares Haverkamp parameters, N trials     %%
%%  F_opt       [outpt] Objective function, N trials                     %%
%%  K_opt       [outpt] # iterations, N trials                           %%
%%                                                                       %%
%%  LITERATURE                                                           %%
%%  Vrugt, J.A., J.W. Hopmans, Y. Gao, M. Rahmati, J. Vanderborght, and  %%
%%      H. Vereecken (2023), The time validity of Philip's two-term      %%
%%      infiltration equation: An elusive theoretical quantity? Vadose   %%
%%      Zone Journal, e20309, pp. 1-25,                                  %%
%%          https://doi.org/10.1002/vzj2.20309                           %%
%%  Vrugt, J.A. and y. Gao (2022), On the three-parameter infiltration   %%
%%      equation of Parlange et al. (1982): Numerical solution,          %%
%%      experimental design, and parameter estimation, Vadose Zone       %%
%%      Journal, 21:e20167, pp. 1-25,                                    %%
%%          https://doi.org/10.1002/vzj2.20167                           %%
%%                                                                       %%
%%  © Written by Jasper A. Vrugt, Jan. 2019                              %%
%%  University of California Irvine                                      %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Time form implies that Ki = 0 --> only works if Ki = 0, thus, d = 3
if d == 4 && app == 2
    error(['least_squares_Haverkamp:Time form of Haverkamp ' ...
        'implies Ki = 0, and, consequently, for d = 4, app ≠ 2. For ' ...
        'app = 2 to work (time minimization), one must use the 3-parameter ' ...
        'variant of Haverkamp''s infiltration equation, then Ki = 0']);
end
if nargin < 8, maxiter = 250; end
if nargin < 7, xtol = 1e-6; end
if nargin < 6, N = 50; end
if nargin < 5, nu = 2; end
if nargin < 4, app = 1; end

% # parameters?
switch d
    case 3       % S   Ks   ß
        x_min = [1e-3 1e-3 1e-3]'; 
        x_max = [ 25   25   2  ]';
    case 4       % S   Ks   ß    Ki
        x_min = [1e-3 1e-3 1e-3 1e-3]'; 
        x_max = [ 25   25   2    25 ]';
    otherwise
        error(['least_squares_Haverkamp:Numbers of parameters ' ...
            'must equal 2, 3 or 4']);
end
% Approach one or two
switch app
    case 1 % Infiltration form: minimizes cum. inf. residuals
        plugin.t = t_meas; f = @(x) Haverkamp_I(x,plugin); y_meas = I_meas; 
    case 2 % Time form: minimizes time residuals
        plugin.I = I_meas; f = @(x) Haverkamp_t(x,plugin); y_meas = t_meas;
end
X = LH_sampling(x_min',x_max',N)';
% Initialize x_opt, K_opt and F_opt
X_opt = nan(d,N); [K_opt,F_opt] = deal(nan(1,N));

for trial = 1:N                                 % Loop over initial guesses
    x = X(1:d,trial);                               % Define initial guess
    Jf = jac_Haverkamp_anal(x,plugin,app);          % nxp Jacobian matrix
    lmd = max(diag(Jf'*Jf));                        % Initial lambda
    dx = inf * ones(d,1);                           % Initlze shift vector
    k = 2;                                          % Initlze loop counter
    while (max(abs(dx)) > xtol) && (k <= maxiter)   % Now while loop
        e = y_meas - f(x(1:d,k-1));                     % nx1 residual vctr
        if k == 2, F_x = sum(e.^2); end                 % objective fnction
        Jf = jac_Haverkamp_anal(x(1:d,k-1),plugin,app); % Jacobian matrix
        D = lmd * diag(diag(Jf'*Jf));                   % dxd damping matrx
        dx = (Jf'*Jf + D)\(Jf'*e);                      % shift vector
        x_try = x(1:d,k-1) + dx;                        % 1xd trial parmtrs
        x_try = check_bound(x_try,x,x_min,x_max,d);     % trial in bound
        e_try = y_meas - f(x_try);                      % nx1 residual vctr
        F_try = sum(e_try.^2);                          % objective fnction
        if F_try < F_x
            x(1:d,k) = x_try; lmd = lmd/nu; F_x = F_try;    % Accept trial
        else
            x(1:d,k) = x(1:d,k-1); lmd = nu*lmd;            % Keep current
        end
        k = k + 1;                                          % Increment ct
    end
    F_opt(1,trial) = F_x;                           % Store objctve fnction
    X_opt(1:d,trial) = x(1:d,k-1);                  % Store parameter vlues
    K_opt(1,trial) = k-1;                           % Store # iterations
end

end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% Secondary function
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% check_bound
function x_try = check_bound(x_try,x,x_min,x_max,d)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Enforces parameter bounds on the trial vector                         %%
%%                                                                       %%
%%  SYNOPSIS: x_try = check_bound(x_try,x,x_min,x_max,d)                 %%
%% where                                                                 %%
%%                                                                       %%
%%  © Written by Jasper A. Vrugt, Jan. 2019                              %%
%%  University of California Irvine                                      %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

var_sigma = 0.05;
for u = 1:d
    if x_try(u) < x_min(u)
        x_try(u) = x_min(u) + var_sigma*(x(u) - x_min(u));
    end
    if x_try(u) > x_max(u)
        x_try(u) = x_max(u) - var_sigma*(x_max(u) - x(u));
    end
end

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
