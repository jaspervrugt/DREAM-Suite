function SimRR = sacsma(x,plugin,model_code)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% RK2 implementation of SAC-SMA
%%
%% Written by JA Vrugt
%% Flux equations based on FUSE paper of Clark et al. (WRR, 2008)
%% Integration based on hmodel of Schoups (2010)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 3
    model_code = 3;         % Formulation/language
end                         % 1: MATLAB Runge Kutta implementation of sacsma_ode
                            % 2: MATLAB ode45 implementation of sacsma_ode
                            % 3: MATLAB Explicit Euler sacsma_ode with int_steps steps
                            % 4: 'C++': Runge Kutta implementation sacsma_ode in C++

nv = 9;                     % Initialize number of state variables
ns = numel(plugin.tout);    % Number of time steps
Y = nan(ns,nv);             % Initialize matrix of state variables
Y(1,1:nv) = plugin.y0(:)';  % Initialize state variables at time 0

% Initialization
switch model_code
    case {1,2,3}
        %% Unpack parameters
        uzfwm = x(1);               % Maximum free water storage of upper zone (mm)
        uztwm = x(2);               % Maximum tension water storage of upper zone (mm)
        lzfpm = x(3);               % Maximum free water storage of lower zone primary (mm)
        lzfsm = x(4);               % Maximum free water storage of lower zone secundary (mm)
        lztwm = x(5);               % Maximum tension water storage of lower zone (mm)
        zperc = x(6);               % Multiplier percolation function (-)
        rexp = x(7);                % Power of percolation function (-)
        uzk = x(8);                 % Interflow rate (1/d)
        pfree = x(9);               % Fraction of percolation to tension storage in lower layer (-)
        lzpk = x(10);               % Base flow depletion rate for primary reservoir (1/d)
        lzsk = x(11);               % Base flow depletion rate for secondary reservoir (1/d)
        acm = x(12);                % Maximum fraction of saturated area (-)
        kf = x(13);                 % Recession constant fast reservoir (1/d)
                                    % --> to compute channel inflow
    case 4
        data.P = plugin.data.P;
        data.Ep = plugin.data.Ep;
        data.uzfwm = x(1);          % Maximum free water storage of upper zone (mm)
        data.uztwm = x(2);          % Maximum tension water storage of upper zone (mm)
        data.lzfpm = x(3);          % Maximum free water storage of lower zone primary (mm)
        data.lzfsm = x(4);          % Maximum free water storage of lower zone secundary (mm)
        data.lztwm = x(5);          % Maximum tension water storage of lower zone (mm)
        data.zperc = x(6);          % Multiplier percolation function (-)
        data.rexp = x(7);           % Power of percolation function (-)
        data.uzk = x(8);            % Interflow rate (1/d)
        data.pfree = x(9);          % Fraction of percolation to tension storage in lower layer (-)
        data.lzpk = x(10);          % Base flow depletion rate for primary reservoir (1/d)
        data.lzsk = x(11);          % Base flow depletion rate for secondary reservoir (1/d)
        data.acm = x(12);           % Maximum fraction of saturated area (-)
        data.kf = x(13);            % Recession constant fast reservoir (1/d)
                                    % --> to compute channel inflow
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
                [ytmp,LTE] = rk2(ytmp,h,uzfwm,uztwm,lzfpm,lzfsm,...
                    lztwm,zperc,rexp,uzk,pfree,lzpk,lzsk,acm,kf,...
                    plugin.data.P(s-1),plugin.data.Ep(s-1));        % Evaluate rk2
                w = 1 ./ (reltol*abs(ytmp) + abstol);               % Weights
                wrms = sqrt(sum((w.*LTE).^2)/nv);                   % CORRECTED, Nov. 2022
                % wrms = sqrt(w*LTE'/nv);                           % Weighted rms of error
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
        ode_options = odeset('InitialStep',1e-5,...    % initial time-step (d)
            'MaxStep',1,...                         % maximum time-step (d)
            'RelTol',1e-3,...                       % relative tolerance
            'AbsTol',1e-3);                         % absolute tolerances (mm)
        [~,Y] = ode45(@sacsma_odefcn,plugin.tout,...
            plugin.y0,ode_options,uzfwm,uztwm,...
            lzfpm,lzfsm,lztwm,zperc,rexp,uzk,...
            pfree,lzpk,lzsk,acm,kf,plugin);

    case 3 %% MATLAB: Explicit Euler sacsma_odefcn with int_steps steps
        int_steps = 20;                             % number of integration steps
        dt = 1/int_steps;                           % time of each int. step
        for s = 2:ns                                % Start time loop
            y = Y(s-1,1:nv);                            % Initialize state variables
            for it = 1:int_steps                        % Do integration in int_steps steps
                dydt = sacsma_odefcn(s-2,y,...              % Compute dxdt based on current state, par and P, Ep
                    uzfwm,uztwm,lzfpm,lzfsm,...
                    lztwm,zperc,rexp,uzk,pfree,...
                    lzpk,lzsk,acm,kf,plugin);
                y = y + dydt' * dt;                         % Update states
            end
            Y(s,1:nv) = y;                              % State at t
        end                                         % End of time loop

    case 4 %% C++: Runge Kutta implementation ( = similar to hymod_odefcn)
        Y = crr_sacsma(plugin.tout,plugin.y0,data,plugin.options)';

end

SimRR = diff(Y(plugin.idx,nv));                 % Now extract discharge

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                   Secundairy functions listed below                   %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% 1. Runge Kutta integration
function [y,LTE] = rk2(y,h,uzfwm,uztwm,lzfpm,lzfsm,lztwm,zperc,...
    rexp,uzk,pfree,lzpk,lzsk,acm,kf,P,Ep)

dydtE = fRhs(y,uzfwm,uztwm,lzfpm,lzfsm,lztwm,zperc,...
    rexp,uzk,pfree,lzpk,lzsk,acm,kf,P,Ep);      % Euler solution
yE = y + h*dydtE;
dydtH = fRhs(yE,uzfwm,uztwm,lzfpm,lzfsm,lztwm,zperc,...
    rexp,uzk,pfree,lzpk,lzsk,acm,kf,P,Ep);      % Heun solution
y = y + 0.5*h*(dydtE + dydtH);                  % New y
LTE = abs(yE - y);                              % Compute estimate of LTE

end

%% 2. SACSMA: conceptual rainfall-runoff model
function [dydt] = fRhs(y,uzfwm,uztwm,lzfpm,lzfsm,lztwm,zperc,rexp,uzk,pfree,lzpk,lzsk,acm,kf,P,Ep)

dydt = nan(1,9);            % Initialize return argument

%% Define smoothing function
eps = 5; rho = 1e-2;        % Dimensionless smoothing coefficients
% WRONG IMPLEMENTATION ACCORDING TO FUSE PAPER (Eq. 12h)
%lf = @(S,Smax) 1 / ( 1 + exp( (S - Smax - rho * Smax * eps)/(rho * Smax) ) );
% CORRECT IMPLEMENTATION ACCORDING TO FUSE SOURCE CODE (GitHub) 
lf = @(S,Smax) 1 / ( 1 + exp( - ( S - (Smax - rho * Smax * eps) )/(rho * Smax) ) );

Lzm = lzfpm + lzfsm + lztwm;   % Maximum storage of lower zone (mm)

%% Unpack states and forcing
Uztw = y(1);                   % Tension water storage upper zone (mm)
Uzfw = y(2);                   % Free water storage upper zone (mm)
Lztw = y(3);                   % Tension water storage lower zone (mm)
Lzps = y(4);                   % Free water storage of lower zone primary (mm)
Lzfs = y(5);                   % Free water storage of lower zone secundary (mm)
Sf = y(6:8);                   % Storage of 1st, 2nd and 3rd fast reservoir (routing, mm)
Lztot = Lztw + Lzps + Lzfs;    % Total lower zone storage (mm)

%% Compute fluxes between layers
e_1 = Ep * min(Uztw,uztwm)/uztwm;       % Evaporation from upper soil layer (mm/d)
e_2 = (Ep-e_1) * min(Lztw,lztwm)/lztwm; % Evaporation from lower soil layer (mm/d)

if e_2 < 0, e_2 = 0; end                % Residual evporation cannot be negative
q_0 = lzpk*lzfpm + lzsk*lzfsm;          % Baseflow at saturation (mm/d)
dlz = 1 + zperc*(Lztot/Lzm)^rexp;       % Lower-zone percolation demand (-)
q_12 = q_0*dlz*(Uzfw/uzfwm);            % Percolation from UZ to LZ (mm/d)
q_if  = uzk*(Uzfw/uzfwm);               % Interflow (mm/d)
q_bp = lzpk*Lzps;                       % Primary baseflow (mm/d)
q_bs = lzsk*Lzfs;                       % Secundary baseflow (mm/d)
ac = Uztw/uztwm * acm;                  % Saturated area (-)
q_sx = ac * P;                          % Saturation excess (mm/d)
q_utof = (P - q_sx) * lf(Uztw,uztwm);   % Overflow of water from tension storage in upper soil layer (mm/d)
q_ufof = q_utof * lf(Uzfw,uzfwm);       % Overflow of water from free storage in the upper soil layer (mm/d)

q_stof = pfree*q_12*lf(Lztw,lztwm);     % Overflow of water from tension storage in the lower soil layer (mm/d)
prc_s = (1-pfree)*q_12/2 + q_stof/2;    % Percolation and q_stof (mm/d)
q_sfofp = prc_s * lf(Lzps,lzfpm);       % Overflow of water from primary base flow storage in lower soil layer (mm/d)
q_sfofs = prc_s * lf(Lzfs,lzfsm);       % Overflow of water from secondary base flow storage in lower soil layer (mm/d)

%% Compute fluxes out of fast reservoirs
q_fout = kf * Sf;                       % Channel inflow is q_fout(3)

%% Compute net flux into/out of reservoir
dydt(1) = P - q_sx - e_1 - q_utof;              % Net flux into/out of S1T  (mm/d)
dydt(2) = q_utof - q_12 - q_if - q_ufof;        % Net flux into/out of S1F  (mm/d)
dydt(3) = pfree*q_12 - e_2 - q_stof;            % Net flux into/out of S2T  (mm/d)
dydt(4) = prc_s - q_bp - q_sfofp;               % Net flux into/out of S2Fa (mm/d)
dydt(5) = prc_s - q_bs - q_sfofs;               % Net flux into/out of S2Fb (mm/d)
dydt(6) = q_if + q_sx + q_ufof + q_sfofp ...
    + q_sfofs - q_fout(1);                      % Net flux into/out of fast reservoir 1 (mm/d)
                                                % NOTE: q_bp and q_bs are not routed (JHM, 2006 paper) 
                                                % But can route q_bp and q_bs together with other fluxes
                                                % Then, add to dxdt(6) and remove from dxdt(9)
dydt(7) = q_fout(1) - q_fout(2);                % Net flux into/out of fast reservoir 2 (mm/d)
dydt(8) = q_fout(2) - q_fout(3);                % Net flux into/out of fast reservoir 3 (mm/d)
dydt(9) = q_bp + q_bs + q_fout(3);              % Inflow (mm/d) to infinite reservoir
                                                % --> yield streamflow from delta                                              
end

%% 3. SACSMA: Secondairy function, ODE solver
function dxdt = sacsma_odefcn(t,x,uzfwm,uztwm,lzfpm,lzfsm,lztwm,...
    zperc,rexp,uzk,pfree,lzpk,lzsk,acm,kf,plugin)

%% Define return argument
dxdt = nan(9,1);            % Initialize return argument

%% Define smoothing function
eps = 5; rho = 0.01;        % Dimensionless smoothing coefficients
% WRONG IMPLEMENTATION ACCORDING TO FUSE PAPER (Eq. 12h)
%lf = @(S,Smax) 1 / ( 1 + exp( (S - Smax - rho * Smax * eps)/(rho * Smax) ) );
% CORRECT IMPLEMENTATION ACCORDING TO FUSE SOURCE CODE (GitHub) 
lf = @(S,Smax) 1 / ( 1 + exp( - ( S - (Smax - rho * Smax * eps) )/(rho * Smax) ) );

Lzm = lzfpm + lzfsm + lztwm;    % Maximum storage of lower zone (mm)

%% Unpack states and forcing
Uztw = x(1);                   % Tension water storage upper zone (mm)
Uzfw = x(2);                   % Free water storage upper zone (mm)
Lztw = x(3);                   % Tension water storage lower zone (mm)
Lzps = x(4);                   % Free water storage of lower zone primary (mm)
Lzfs = x(5);                   % Free water storage of lower zone secundary (mm)
Sf = x(6:8);                   % Storage of 1st, 2nd and 3rd fast reservoir (routing, mm)
Lztot = Lztw + Lzps + Lzfs;    % Total lower zone storage (mm)
it = floor(t) + 1;             % Truncate time to current time index
P = plugin.data.P(it,1);       % Get current rainfall (mm/d)
Ep = plugin.data.Ep(it,1);     % Get current Ep (mm/d)

%% Compute fluxes between layers
e_1 = Ep * min(Uztw,uztwm)/uztwm;       % Evaporation from upper soil layer (mm/d)
e_2 = (Ep-e_1) * min(Lztw,lztwm)/lztwm; % Evaporation from lower soil layer (mm/d)
q_0 = lzpk*lzfpm + lzsk*lzfsm;          % Baseflow at saturation (mm/d)
dlz = 1 + zperc*(Lztot/Lzm)^rexp;       % Lower-zone percolation demand (-)
q_12 = q_0*dlz*(Uzfw/uzfwm);            % Percolation from UZ to LZ (mm/d)
q_if  = uzk*(Uzfw/uzfwm);               % Interflow (mm/d)
q_bp = lzpk*Lzps;                       % Primary baseflow (mm/d)
q_bs = lzsk*Lzfs;                       % Secundary baseflow (mm/d)
ac = Uztw/uztwm * acm;                  % Saturated area (-)
q_sx = ac * P;                          % Saturation excess (mm/d)
q_utof = (P - q_sx) * lf(Uztw,uztwm);   % Overflow of water from tension storage in upper soil layer (mm/d)
q_ufof = q_utof * lf(Uzfw,uzfwm);       % Overflow of water from free storage in the upper soil layer (mm/d)

if e_2 < 0, e_2 = 0; end                % Residual evporation cannot be negative
q_stof = pfree*q_12*lf(Lztw,lztwm);     % Overflow of water from tension storage in the lower soil layer (mm/d)
prc_s = (1-pfree)*q_12/2 + q_stof/2;    % Percolation and q_stof (mm/d)
q_sfofp = prc_s * lf(Lzps,lzfpm);       % Overflow of water from primary base flow storage in lower soil layer (mm/d)
q_sfofs = prc_s * lf(Lzfs,lzfsm);       % Overflow of water from secondary base flow storage in lower soil layer (mm/d)

%% Compute fluxes out of fast reservoirs
q_fout = kf * Sf;                       % Channel inflow is q_fout(3)

%% Compute net flux into/out of reservoir
dxdt(1) = P - q_sx - e_1 - q_utof;               % Net flux into/out of S1T  (mm/d)
dxdt(2) = q_utof - q_12 - q_if - q_ufof;        % Net flux into/out of S1F  (mm/d)
dxdt(3) = pfree*q_12 - e_2 - q_stof;             % Net flux into/out of S2T  (mm/d)
dxdt(4) = prc_s - q_bp - q_sfofp;               % Net flux into/out of S2Fa (mm/d)
dxdt(5) = prc_s - q_bs - q_sfofs;               % Net flux into/out of S2Fb (mm/d)
dxdt(6) = q_if + q_sx + q_ufof + q_sfofp ...
    + q_sfofs - q_fout(1);                      % Net flux into/out of fast reservoir 1 (mm/d)
                                                % NOTE: q_bp and q_bs are not routed (JHM, 2006 paper) 
                                                % But can route q_bp and q_bs together with other fluxes
                                                % Then, add to dxdt(6) and remove from dxdt(9)
dxdt(7) = q_fout(1) - q_fout(2);                % Net flux into/out of fast reservoir 2 (mm/d)
dxdt(8) = q_fout(2) - q_fout(3);                % Net flux into/out of fast reservoir 3 (mm/d)
dxdt(9) = q_bp + q_bs + q_fout(3);              % Inflow (mm/d) to infinite reservoir
                                                % --> yield streamflow from delta
end

% QUESTIONS: 1. Do we have to account for evaporation of channel inflow?
%            2. Do we have to treat impervious fraction ( = state variable)
%               instead of saturation excess?
% NOTE:      1. Current formulation has 13 parameters instead of 14 as
%               impervious fraction is not modeled. This has a state
%               variable, evaporation, impervious area and direct runoff
%            2. q_if + q_sx + q_ufof + q_sfofa + q_sfofb is routed through 
%               three linear fast reservoirs with same recession constant
