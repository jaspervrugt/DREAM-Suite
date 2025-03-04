function SimRR = hmodel(x,plugin,model_code)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Runge Kutta implementation of Hmodel
%%
%% Written by JA Vrugt based on C++ code of Schoups
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 3
    model_code = 3;         % Formulation/language
end                         % 1: MATLAB Runge Kutta implementation of hmodel_ode
                            % 2: MATLAB ode45 implementation of hmodel_ode
                            % 3: MATLAB Explicit Euler hmodel_ode with int_steps steps
                            % 5: 'C++': Runge Kutta implementation hmodel_ode in C++

nv = 5;                     % Initialize number of state variables
ns = numel(plugin.tout);    % Number of time steps
Y = nan(ns,nv);             % Initialize matrix of state variables 
Y(1,1:nv) = plugin.y0(:)';  % Initialize state variables at time 0

% Initialization
switch model_code
    case {1,2,3}
        I_max  = x(1);      % Interception storage capacity (mm)
        Su_max = x(2);      % Unsaturated zone storage capacity (mm)
        Qs_max = x(3);      % Maximum percolation rate (mm/d)
        aE    = x(4);       % Evaporation coefficient
        aF    = x(5);       % Runoff coefficient
        aS    = 1e-6;       % Percolation coefficient
        Kf    = x(6);       % Fast-flow response time (d)
        Ks    = x(7);       % Slow-flow response time (d)
    case 5
        data.P = plugin.data.P';
        data.Ep = plugin.data.Ep';
        data.I_max  = x(1);     % interception storage capacity (mm)
        data.Su_max = x(2);     % unsaturated zone storage capacity (mm)
        data.Qs_max = x(3);     % maximum percolation rate (mm/d)
        data.aE    = x(4);      % evaporation coefficient
        data.aF    = x(5);      % runoff coefficient
        data.aS    = 1e-6;      % percolation coefficient
        data.Kf    = x(6);      % fast-flow response time (d)
        data.Ks    = x(7);      % slow-flow response time (d)

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
                ytmp = Y(s,1:nv);                                   % Advance solution by step h
                [ytmp,LTE] = rk2(ytmp,h,I_max,Su_max,Qs_max,aE,aF,aS, ...
                    Kf,Ks,plugin.data.P(s-1),plugin.data.Ep(s-1));  % Evaluate rk2
                w = 1 ./ (reltol*abs(ytmp) + abstol);               % Weights
                wrms = sqrt(sum((w.*LTE).^2)/nv);                   % CORRECTED, Nov. 2022
%                wrms = sqrt(w*LTE'/nv);                             % Weighted rms of error
                if (h <= hmin_), wrms = 0.5; end                    % Make sure time step increases after acceptance
                if (wrms <= 1) || (h <= hmin_)                      % Accept if error is small enough
                    Y(s,1:nv) = ytmp; t = t + h;
                end
                h = h*max(0.2,min(5.0,0.9*wrms^(-1/order)));        % Compute new step
                h = max(hmin_,min(h,hmax_)); h = min(h,t2-t);       % Another turn
            end
        end

    case 2 %% MATLAB: ode45 implementation of hymod_odefcn
        plugin.tout(end) = plugin.tout(end)-1e-10;  % Loop ode goes to maxT and then still one try
                                                    % Then, time index goes out of bound of forcing data
        ode_options = odeset('InitialStep',1,...    % initial time-step (d)
            'MaxStep',1,...                         % maximum time-step (d)
            'RelTol',1e-3,...                       % relative tolerance
            'AbsTol',1e-3);                         % absolute tolerances (mm)
        [~,Y] = ode23(@hmodel_odefcn,plugin.tout,...
            plugin.y0,ode_options,I_max,Su_max, ...
            Qs_max,aE,aF,aS,Kf,Ks,plugin);

    case 3 %% MATLAB: Explicit Euler hymod_odefcn with int_steps steps
        int_steps = 10;                             % number of integration steps
        dt = 1/int_steps;                           % time of each int. step                   
        for s = 2:ns                                % Start time loop
            y = Y(s-1,1:nv);                            % Initialize state variables
            for it = 1:int_steps                        % Do integration in int_steps steps
                dydt = hmodel_odefcn(s-2,y,...              % Compute dxdt based on current state, par and P, Ep
                    I_max,Su_max,Qs_max,aE,aF, ...
                    aS,Kf,Ks,plugin);     
                y = y + dydt' * dt;                         % Update states
            end
            Y(s,1:nv) = y;                              % State at t                                
        end                                         % End of time loop

    case 5 %% C++: Runge Kutta implementation ( = similar to hymod_odefcn)
        Y = crr_hmodel(plugin.tout,plugin.y0,data,plugin.options)';

end
SimRR = diff(Y(plugin.idx,nv));                 % Now extract discharge

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                   Secundairy functions listed below                   %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% 1. Runge Kutta integration
function [y,LTE] = rk2(y,h,I_max,Su_max,Qs_max,aE,aF,aS,Kf,Ks,P,Ep)

dydtE = fRhs(y,I_max,Su_max,Qs_max,aE,aF,aS,Kf,Ks,P,Ep); 
yE = y + h*dydtE;                                           % Euler solution
dydtH = fRhs(yE,I_max,Su_max,Qs_max,aE,aF,aS,Kf,Ks,P,Ep);   % Evaluate dx
y = y + 0.5*h*(dydtE + dydtH);                              % Heun solution
LTE = abs(yE - y);                                          % Compute estimate of LTE

end

%% 2. HMODEL: conceptual rainfall-runoff model
function [dydt] = fRhs(y,I_max,Su_max,Qs_max,aE,aF,aS,Kf,Ks,P,Ep)

% expFlux = @(sr,alfa) (1 - exp(-alfa*sr)) / (1 - exp(-alfa));

Si = y(1);              % Interception storage (mm)
Su = y(2);              % Surface storage (mm)
Sf = y(3);              % Storage of fast/quick reservoirs (mm)
Ss = y(4);              % Storage of slow reservoir (mm)
Sur = min(Su/Su_max,1); % Truncate relative storage, otherwise Perc
                        % This is a brute-force solution, can replace
                        % with an exponential flux relationship, or
                        % alternatively, use present solution
if I_max > 0
    Sir = Si/I_max;                 % Ratio of interception storage to max storage
    EvapI  = Ep * expFlux(Sir,50);  % Interception evaporation (mm/d)
    P_e    = P * expFlux(Sir,-50);  % Effective precipitation (mm/d)
    Ep_e   = max(0.,Ep - EvapI);    % Effective evaporation (mm/d)
else
    EvapI = 0;                      % Interception evaporation (mm/d)
    P_e = P;                        % Effective precipitation (mm/d)
    Ep_e = Ep;                      % Effective evaporation (mm/d)
end
Ea = Ep_e * expFlux(Sur,aE);        % Actual evaporation (mm/d)
prc = Qs_max * expFlux(Sur,aS);     % Percolation (mm/d)
rnf = P_e * expFlux(Sur,aF);        % Surface runoff (mm/d)
qf  = Sf/Kf;                        % Outflow of fast reservoir (mm/d)
qs  = Ss/Ks;                        % Outflow of slow reservoir (mm/d)

dydt(1) = P - EvapI - P_e;          % Flux into/out of interception reservoir (mm/d)
dydt(2) = P_e - Ea - prc - rnf;     % Flux into/out of unsaturated reservoir (mm/d)    
dydt(3) = rnf - qf;                 % Flux into/out of fast reservoir (mm/d)
dydt(4) = prc - qs;                 % Flux into/out of slow reservoir (mm/d)
dydt(5) = qf + qs;                  % Inflow (mm/d) to infinite reservoir
                                    % --> yield streamflow from delta
end

%% 3. HMODEL: Secondairy function, ODE solver
function dxdt = hmodel_odefcn(t,x,I_max,Su_max,Qs_max,aE,aF,aS,Kf,Ks,plugin)

dxdt = nan(5,1);            % Initialize return argument
Si = x(1);                  % Interception storage (mm)
Su = x(2);                  % Surface storage (mm)
Sf = x(3);                  % Storage of fast/quick reservoirs (mm)
Ss = x(4);                  % Storage of slow reservoir (mm)
it = floor(t) + 1;          % Truncate time to current time index
P = plugin.data.P(it,1);    % Get current rainfall (mm/d)
Ep = plugin.data.Ep(it,1);  % Get current Ep (mm/d)
Sur = Su/Su_max;            % Ratio of unsaturated storage to max storage

if I_max > 0
    Sir = Si/I_max;                 % Ratio of interception storage to max storage
    EvapI  = Ep * expFlux(Sir,50);  % Interception evaporation (mm/d)
    P_e    = P * expFlux(Sir,-50);  % Effective precipitation (mm/d)
    Ep_e   = max(0.,Ep - EvapI);    % Effective evaporation (mm/d)
else
    EvapI = 0;                      % Interception evaporation (mm/d)
    P_e = P;                        % Effective precipitation (mm/d)
    Ep_e = Ep;                      % Effective evaporation (mm/d)
end
Ea = Ep_e * expFlux(Sur,aE);        % Actual evaporation (mm/d)
prc = Qs_max * expFlux(Sur,aS);     % Percolation (mm/d)
rnf = P_e * expFlux(Sur,aF);        % Surface runoff (mm/d)
qf  = Sf/Kf;                        % Outflow of fast reservoir (mm/d)
qs  = Ss/Ks;                        % Outflow of slow reservoir (mm/d)

dxdt(1) = P - EvapI - P_e;          % Flux into/out of interception reservoir (mm/d)
dxdt(2) = P_e - Ea - prc - rnf;     % Flux into/out of unsaturated reservoir (mm/d)    
dxdt(3) = rnf - qf;                 % Flux into/out of fast reservoir (mm/d)
dxdt(4) = prc - qs;                 % Flux into/out of slow reservoir (mm/d)
dxdt(5) = qf + qs;                  % Inflow (mm/d) to infinite reservoir
                                    % --> yield streamflow from delta
end

%% Exponential flux
function Qr = expFlux(Sr,a)
% Relative flux from exponential storage-flux relation
Sr = max(0.0,min(1.0,Sr));
if abs(a) < 1e-6
    Qr = Sr;    % approximately linear
else
    Qr = (1.-exponen(-a*Sr))/(1.-exponen(-a));
end

%% Exponent
function f = exponen(x)
% Exponential function with protection against overflow
f = exp(min(300,x));
end
end