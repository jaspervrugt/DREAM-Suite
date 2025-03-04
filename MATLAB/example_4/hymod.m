function SimRR = hymod(x,plugin,model_code)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Runge Kutta implementation of HYMOD
%%
%% Written by JA Vrugt based on C++ code of Schoups
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 3
    model_code = 3;         % Formulation/language
end                         % 1: MATLAB Runge Kutta implementation of hymod_ode
                            % 2: MATLAB ode45 implementation of hymod_ode
                            % 3: MATLAB Explicit Euler hymod_ode with int_steps steps
                            % 4: MATLAB old implementation
                            % 5: 'C++': Runge Kutta implementation hymod_ode in C++

nv = 6;                     % Initialize number of state variables
ns = numel(plugin.tout);    % Number of time steps
Y = nan(ns,nv);             % Initialize matrix of state variables 
Y(1,1:nv) = plugin.y0(:)';  % Initialize state variables at time 0
        
% Initialization
switch model_code
    case {1,2,3,4}
        Su_max = x(1);              % Maximum storage of unsaturated zone
        beta  = x(2);               % Spatial variability of soil moisture capacity
        alfa  = x(3);               % Flow partitioning coefficient
        Ks    = x(4);               % Recession constant slow reservoir
        Kf    = x(5);               % Recession constant fast reservoir
    case 5
        data.P = plugin.data.P;
        data.Ep = plugin.data.Ep;
        data.m = 1e-2;              % Smoothing parameter
        data.Sumax = x(1);          % Maximum storage of unsaturated zone
        data.beta  = x(2);          % Spatial variability of soil moisture capacity
        data.alfa  = x(3);          % Flow partitioning coefficient
        data.Ks    = x(4);          % Recession constant slow reservoir
        data.Kf    = x(5);          % Recession constant fast reservoir
end

% Execute model
switch model_code

    case 1 %% MATLAB: Runge Kutta implementation (= model similar to hymod_odefcn) 
        hin    = plugin.options.InitialStep;        % Initial time step
        hmax_  = plugin.options.MaxStep;            % Maximum time step
        hmin_  = plugin.options.MinStep;            % minimum time step
        reltol = plugin.options.RelTol;             % Relative tolerance
        abstol = plugin.options.AbsTol;             % Absolute tolerance
        order  = plugin.options.Order;              % Order
        for s = 2:ns                                % Integrate from tprint(1) to tprint(end)
            t1 = plugin.tout(s-1); t2 = plugin.tout(s); % Set start and end times
            h = hin;                                    % Set initial step
            h = max(hmin_,min(h,hmax_));                % Make sure within range
            h = min(h,t2-t1);                           % Determine based on t2-t1
            Y(s,1:nv) = Y(s-1,1:nv);                    % Initial state, u
            t = t1;                                     % Initial time
            %% Integrate from t1 to t2
            while (t < t2)
                ytmp = Y(s,1:nv);                                   % Initial states at t 
                [ytmp,LTE] = rk2(ytmp,h,Su_max,beta,alfa,Ks,Kf,...  % New state at t+h, return int. error
                    plugin.data.P(s-1),plugin.data.Ep(s-1));        
                w = 1 ./ (reltol*abs(ytmp) + abstol);               % Weights
                wrms = sqrt(w*LTE'/nv);                             % Weighted rms of error
                if (h <= hmin_), wrms = 0.5; end                    % Make sure time step increases after acceptance
                if (wrms <= 1) || (h <= hmin_)                      % Accept if error is small enough
                    Y(s,1:nv) = ytmp; t = t + h;
                end
                h = h*max(0.2,min(5.0,0.9*wrms^(-1/order)));        % Compute new integration step
                h = max(hmin_,min(h,hmax_)); h = min(h,t2-t);       % Check it satisfies constraints
            end
        end

    case 2 %% MATLAB: ode45 implementation of hymod_odefcn
        plugin.tout(end) = plugin.tout(end)-1e-10;  % Loop ode goes to maxT and then still one try
                                                    % Then, time index goes out of bound of forcing data
        ode_options = odeset(...                    % Define solver settings
            'InitialStep',plugin.options.InitialStep,...    % initial time-step (d)
            'MaxStep',plugin.options.MaxStep,...            % maximum time-step (d)
            'RelTol',plugin.options.RelTol,...              % relative tolerance
            'AbsTol',plugin.options.AbsTol);                % absolute tolerances (mm)
        [~,Y] = ode23(@hymod_odefcn,plugin.tout,...
            plugin.y0,ode_options,Su_max,beta,...
            alfa,Ks,Kf,plugin);

    case 3 %% MATLAB: Explicit Euler hymod_odefcn with int_steps steps
        int_steps = 10;                             % number of integration steps
        dt = 1/int_steps;                           % time of each int. step                   
        for s = 2:ns                                % Start time loop
            y = Y(s-1,1:nv);                            % Initialize state variables
            for it = 1:int_steps                        % Do integration in int_steps steps
                dydt = hymod_odefcn(s-2,y,...               % Compute dxdt based on current state, par and P, Ep
                    Su_max,beta,alfa,Ks,Kf,plugin);     
                y = y + dydt' * dt;                         % Update states
            end
            Y(s,1:nv) = y;                              % State at t                                
        end                                         % End of time loop

    case 4 %% MATLAB: Old implementation
        bexp = beta;                                % Reminder that bexp is beta in 'old' implementation
        for s = 2:ns                                % Start time loop
            y = Y(s-1,1:nv);                            % Initialize state variables
            [P_e,y(1)] = excess(y(1),Su_max,bexp,...    % Compute excess precipitation and evaporation
                plugin.data.P(s-1,1),...
                plugin.data.Ep(s-1,1));                     
            qf_in = alfa * P_e; qs_in = (1-alfa) * P_e; % Now partition ER between quick and slow flow reservoirs
            [y(2),qs_out] = linres(y(2),qs_in,Ks);      % Route slow flow component with single linear reservoir
            for k = 3:5                                 % Route quick flow component with linear reservoirs
                [y(k),qf_out] = linres(y(k),qf_in,Kf);      % Linear reservoir
                qf_in = qf_out;                             % Outflow equals inflow of next fast reservoir
            end 
            y(6) = y(6) + qf_out + qs_out;              % Total inflow to discharge reservoir             
            Y(s,1:nv) = y;                              % State variables at end of time interval
        end                                        % End of time loop  

    case 5 %% C++: Runge Kutta implementation ( = similar to hymod_odefcn)
        Y = crr_hymod(plugin.tout,plugin.y0,data,plugin.options)';

end

SimRR = diff(Y(plugin.idx,nv));                 % Now extract discharge

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                   Secundairy functions listed below                   %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% 1. Runge Kutta integration
function [y,LTE] = rk2(y,h,Su_max,beta,alfa,Ks,Kf,P,Ep)

dydtE = fRhs(y,Su_max,beta,alfa,Ks,Kf,P,Ep); 
yE = y + h*dydtE;                               % Euler solution
dydtH = fRhs(yE,Su_max,beta,alfa,Ks,Kf,P,Ep);   % Evaluate dx
y = y + 0.5*h*(dydtE + dydtH);                  % Heun solution
LTE = abs(yE - y);                              % Integration error

end

%% 2. HYMOD: conceptual rainfall-runoff model
function [dydt] = fRhs(y,Su_max,beta,alfa,Ks,Kf,P,Ep)

m     = 1e-2;           % Smoothing parameter
% expFlux = @(sr,alfa) (1 - exp(-alfa*sr)) / (1 - exp(-alfa));

Su = y(1);              % Storage in unsaturated zone reservoir
Ss = y(2);              % Storage in slow reservoir
Sf = y(3:5);            % Storage in three fast reservoirs
Sur = min(Su/Su_max,1); % Truncate relative storage, otherwise Perc
                        % This is a brute-force solution, can replace
                        % with an exponential flux relationship, or
                        % alternatively, use present solution
Perc = P * (1 - (1-Sur)^beta );         % Sur should not exceed one
Ea = Ep * Sur * (1 + m)/(Sur + m);      % This is generally OK
%Ea = Ep * expFlux(Su/Sumax,aE);        % Numerically more stable
dydt(1) = P - Perc - Ea;                % Eff. flux into first reservoir

qs_in = (1-alfa) * Perc;                % Flux into slow reservoir
qs_out = Ks * Ss;                       % Flux out of slow reservoir
dydt(2) = qs_in - qs_out;               % Eff. flux into second reservoir

qf_in = alfa * Perc;                    % Flux into first fast reservoir
qf_out = Kf * Sf;                       % Fluxes out of fast reservoirs

dydt(3) = qf_in - qf_out(1);            % Eff. flux into 1st fast reservoir
dydt(4) = qf_out(1) - qf_out(2);        % Eff. flux into 2nd fast reservoir
dydt(5) = qf_out(2) - qf_out(3);        % Eff. flux into 3rd fast reservoir
dydt(6) = qs_out + qf_out(3);           % Eff. flux into infinite reservoir
                                        % --> yield streamflow from delta
end

%% 3. HYMOD: Secondairy function, ODE solver
function dxdt = hymod_odefcn(t,x,Su_max,beta,alfa,Ks,Kf,plugin)

dxdt = nan(6,1);            % Initialize return argument
m = 1e-2;                   % Smoothing coefficient (mm)
Su = x(1);                  % Surface storage (mm)
Ss = x(2);                  % Storage of slow reservoir (mm)
Sf = x(3:5);                % Storage of fast/quick reservoirs (mm)
it = floor(t) + 1;          % Truncate time to current time index
P = plugin.data.P(it,1);    % Get current rainfall (mm/d)
Ep = plugin.data.Ep(it,1);  % Get current Ep (mm/d)
Sur = min(Su/Su_max,1);     % Ratio of surface storage to max storage
qu = P*(1 - (1-Sur)^beta);  % Precipitation converted to flow, qu (mm/d)
Ea = Ep * Sur ...
    * (1 + m)/(Sur + m);    % Actual evaporation (mm/d)
dxdt(1) = P - qu - Ea;      % Flux into/out of unsaturated reservoir (mm/d)

qs = (1-alfa) * qu;         % Inflow to slow reservoir (mm/d)
qs_o = Ks * Ss;             % Flow out of slow reservoir (mm/d)
dxdt(2) = qs - qs_o;        % Flux into/out of slow reservoir (mm/d)

qf = alfa * qu;             % Inflow to first quick reservoir (mm/d)
qf_o = Kf * Sf;             % Flow out of three quick reservoirs (mm/d)
dxdt(3) = qf - qf_o(1);     % Flux into/out of first fast reservoir (mm/d)
dxdt(4) = qf_o(1)-qf_o(2);  % Flux into/out of 2nd fast reservoir (mm/d)
dxdt(5) = qf_o(2)-qf_o(3);  % Flux into/out of 3rd fast reservoir (mm/d)
dxdt(6) = qf_o(3) + qs_o;   % Inflow (mm/d) to infinite reservoir
                            % --> yield streamflow from delta
end

%% 4. Excess precipitation and evaporation from UZ
function [P_e,Su_n] = excess(Su,Su_max,bexp,P,Ep)

Su_old = Su; b1exp = bexp + 1;
prev = Su_max * (1 - ( 1 - b1exp * Su/Su_max )^(1/b1exp));
% Calculate effective rainfall 1 
P_e1 = max(P - Su_max + prev , 0); P = P - P_e1;

dummy = min((prev + P)/Su_max,1);
Su_n = Su_max/b1exp * ( 1 - (1 - dummy)^b1exp );
% Calculate effective rainfall 2
P_e2 = max(P - (Su_n - Su_old) , 0);
% Alternative approach: ET linearly related to soil moisture storage
evap = Ep * (1 - (Su_max/b1exp - Su_n) / (Su_max / b1exp) ); 
Su_n = max ( Su_n - evap , 0); % update state
%evap = min(xn,PETval); xn = xn-evap;
% Compute effective rainfall
P_e = P_e1 + P_e2;

end

%% 5. Linear reservoir
function [Sn,q_out] = linres(S,q_in,K)
% Linear reservoir
% S:    state of reservoir (L) 
% q_in  inflow (L/T)
% K     recession constant (1/L)
% with imaginary time step of 1 unit of L
% itu = 1;

% Unfamiliar implementation: Implicit for fixed time step of 1 unit 
Sn = (1-K)*S + (1-K)*q_in;    % New state of reservoir (L)
q_out = (K/(1-K))*Sn;         % Flux out of reservoir (L)
% --> adaptation is numerically more appropriate when itu is large

% Conventional method: instaneous flux integrated over itu units
% q_out = K * S;                  % Flux out of reservoir (L/T)
% Sn = S + (q_in - q_out)*itu;    % New state of reservoir (L)
% Should yield similar solution as above if we repeat n times per itu

end