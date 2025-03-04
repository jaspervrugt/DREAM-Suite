function J = jac_Haverkamp_num(x,plugin,app)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Numerical approximation of the nxd Jacobian matrix of the Haverkamp   %%
%% infiltration equation in time and infiltration form                   %%
%%                                                                       %%
%%  SYNOPSIS: J = jac_Haverkamp_anal(x,plugin,app)                       %%
%% where                                                                 %%
%%  x           [input] 1xd vector of Haverkamp parameters               %%
%%                  3-par form, S [cm/h^.5], Ks [cm/h], ß [-]            %%
%%                  4-par form, S [cm/h^.5], Ks [cm/h], ß [-], Ki [cm/h] %%
%%  plugin      [input] structure time, cum. inf. experimental info      %%
%%  app         [input] Residuals to minimimze                           %%
%%   = 1            Minimize cumulative infiltration residuals           %%
%%   = 2            Minimize time residuals [valid if Ki = 0 -> d = 3]   %%
%%  J           [outpt] nxd Jacobian matrix of time/inf form Haverkamp   %%
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

% Brute-force implementation of differencing. Better methods are available
% but inconsequential for demonstration in VZJ paper

% Check wkether we estimate Ki or not
d = numel(x);
% Check input arguments
if d == 4 && app == 2
    error(['jac_Haverkamp_num:Time form of Haverkamp implies Ki = 0, ' ...
        'and, consequently, for d = 4, app ≠ 2. For app = 2 to work ' ...
        '(time minimization), one must use the 3-parameter variant of ' ...
        'Haverkamp''s infiltration equation, then Ki = 0']);
end
% Multiplier of parameter increments
kappa = 1e-3;           

switch app
    case 1 % Infiltration form
        [dx1,dx2] = deal(nan(1,3));
        % Need numerical solution: Two sided intervals with kappa
        for j = 1:d
            % Define x as default x - and change jth coordinate: -
            x1 = x;
            % Compute delta_eta
            dx1(j) = kappa*x1(j);
            % Update jth parameter
            x1(j) = x1(j) - dx1(j);
            % No check against parameter bounds!!
            % Evaluate implicit solution
            I_x1 = Haverkamp_I(x1,plugin); 
            % Initialize nxd Jacobian matrix
            if j == 1
                n = numel(I_x1); J = nan(n,d);
            end
            % Define x as default x - and change jth coordinate: +
            x2 = x;
            % Compute dx2
            dx2(j) = kappa*x2(j);
            % Update jth parameter
            x2(j) = x2(j) + dx2(j);
            % No check against parameter bounds!!
            % Evaluate implicit solution
            I_x2 = Haverkamp_I(x2,plugin);
            % Now compute Jacobian
            J(1:n,j) = (I_x2 - I_x1)/(dx1(j)+dx2(j));
        end
    case 2 % Time form
        [dx1,dx2] = deal(nan(1,3));
        % Need numerical solution: Two sided intervals with kappa;
        for j = 1:numel(x)
            % Define eta as default eta - and change jth coordinate: min
            x1 = x;
            % Compute delta_eta
            dx1(j) = kappa*x1(j);
            % Update jth parameter
            x1(j) = x1(j) - dx1(j);
            % No check against parameter bounds!!
            % Evaluate implicit solution
            t_x1 = Haverkamp_t(x1,plugin);
            % Initialize nxd Jacobian matrix
            if j == 1
                n = numel(t_x1); J = nan(n,d);
            end
            % Define eta as default eta - and change jth coordinate: min
            x2 = x;
            % Compute delta_eta
            dx2(j) = kappa*x2(j);
            % Update jth parameter
            x2(j) = x2(j) + dx2(j);
            % No check against parameter bounds!!
            % Evaluate implicit solution
            t_x2 = Haverkamp_t(x2,plugin);
            % Now compute Jacobian
            J(1:n,j) = (t_x2 - t_x1)/(dx1(j)+dx2(j));
        end
end

end
