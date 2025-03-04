import numpy as np
from scipy.integrate import solve_ivp

############################
### Case study 23
############################
def sacsma(x, plugin, model_code = 3):
    # ####################################################################### #
    # SACSMA Different numerical solutions of SACSMA model of Burnash, 1973   #
    # This forumulation has been presented in Vrugt, 2024                     #
    # Formulation/language                                                    #
    #   1: Python Runge Kutta implementation of sacsma_ode                    #
    #   2: Python ode45 implementation of sacsma_ode                          #
    #   3: Python MATLAB Explicit Euler sacsma_ode with int_steps steps       #
    #   4: 'C++': Runge Kutta implementation sacsma_ode in C++                #
    # ####################################################################### #

    nv = 9                          # Number of state variables
    ns = len(plugin['tout'])        # Number of time steps
    Y = np.full((ns, nv), np.nan)   # Matrix for state variables
    Y[0, :nv] = plugin['y0'].copy() # Initialize state variables at time 0

    # Unpacking parameters
    if model_code in [1, 2, 3]:
        S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, nu_p, nu_s, a_cmax, kf = x[:13]
    elif model_code == 4:
        data = plugin['data']
        data['P'] = data['P']
        data['Ep'] = data['Ep']
        data.update({               # Unpack parameters into the data dict
            'uzfwm': x[0],          # Maximum free water storage of upper zone (mm)
            'uztwm': x[1],          # Maximum tension water storage of upper zone (mm)
            'lzfpm': x[2],          # Maximum free water storage of lower zone primary (mm)
            'lzfsm': x[3],          # Maximum free water storage of lower zone secundary (mm)
            'lztwm': x[4],          # Maximum tension water storage of lower zone (mm)
            'zperc': x[5],          # Multiplier percolation function (-)
            'rexp': x[6],           # Power of percolation function (-)    
            'uzk': x[7],            # Interflow rate (1/d)    
            'pfree': x[8],          # Fraction of percolation to tension storage in lower layer (-)
            'lzpk': x[9],           # Base flow depletion rate for primary reservoir (1/d)
            'lzsk': x[10],          # Base flow depletion rate for secondary reservoir (1/d)
            'acm': x[11],           # Maximum fraction of saturated area (-)
            'kf': x[12]             # Recession constant fast reservoir (1/d)
        })

    # Execute the model based on the specified model_code
    if model_code == 1:  # Runge-Kutta method
        hin = plugin['options']['InitialStep']
        hmax_ = plugin['options']['MaxStep']
        hmin_ = plugin['options']['MinStep']
        reltol = plugin['options']['RelTol']
        abstol = plugin['options']['AbsTol']
        order = plugin['options']['Order']
        
        for s in range(1, ns):
            t1, t2 = plugin['tout'][s-1], plugin['tout'][s]
            h = hin
            h = max(hmin_, min(h, hmax_))
            h = min(h, t2 - t1)
            Y[s, :nv] = Y[s-1, :nv].copy()      # Set initial state
            t = t1
            
            # Integrate from t1 to t2
            while t < t2:
                ytmp = Y[s, :nv]
                ytmp, LTE = rk2_sac(ytmp, h, S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, 
                                nu_p, nu_s, a_cmax, kf, plugin['data']['P'][s-1], plugin['data']['Ep'][s-1])
                w = 1 / (reltol * np.abs(ytmp) + abstol)
                wrms = np.sqrt(np.sum((w * LTE) ** 2) / nv)
                if h <= hmin_:
                    wrms = 0.5
                if wrms <= 1 or h <= hmin_:
                    Y[s, :nv] = ytmp.copy()
                    t += h
                h = h * max(0.2, min(5.0, 0.9 * wrms ** (-1 / order)))
                h = max(hmin_, min(h, hmax_))
                h = min(h, t2 - t)
        
    elif model_code == 2:  # ode45 implementation
        plugin['tout'][-1] -= 1e-10
        ode_options = { 'rtol': plugin['options']['RelTol'],        # Relative tolerance
                        'atol': plugin['options']['AbsTol'],        # Absolute tolerance
                        'max_step': plugin['options']['MaxStep'],   # Maximum step size
                        'min_step': plugin['options']['MinStep']}   # Minimum step size (optional)
        
        def ode_func(t, y):
            return sacsma_odefcn(t, y, S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, 
                                  nu_p, nu_s, a_cmax, kf, plugin)

        solution = solve_ivp(ode_func, (plugin['tout'][0], plugin['tout'][-1]), plugin['y0'], t_eval=plugin['tout'], **ode_options)
        Y = solution.y.T
        
    elif model_code == 3:  # Explicit Euler method - 20 integration steps for Δt = 1
        int_steps = 20
        dt = 1 / int_steps
        for s in range(1, ns):
            y = Y[s-1, :nv].copy()
            for _ in range(int_steps):
                dydt = sacsma_odefcn(s-1, y, S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, 
                                      nu_p, nu_s, a_cmax, kf, plugin)
                y = y + dydt * dt
            Y[s, :nv] = y.copy()

    # elif model_code == 4:  # C++ Runge-Kutta in C++ using dll [= fix later]
    #    Y = crr_sacsma(plugin['tout'], plugin['y0'], data, plugin['options']).T
    SimRR = np.diff(Y[plugin['idx'], nv-1])  # Discharge (inferred output)

    return SimRR


# Runge-Kutta integration
def rk2_sac(y, h, S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, nu_p, nu_s, a_cmax, kf, P, Ep):
    dydtE = fRhs_sac(y, S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, nu_p, nu_s, a_cmax, kf, P, Ep)
    yE = y + h * dydtE
    dydtH = fRhs_sac(yE, S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, nu_p, nu_s, a_cmax, kf, P, Ep)
    y = y + 0.5 * h * (dydtE + dydtH)
    LTE = np.abs(yE - y)

    return y, LTE


# SACSMA conceptual model flux calculation
def fRhs_sac(y, S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, nu_p, nu_s, a_cmax, kf, P, Ep):
    # Define smoothing function (corrected)
    eps, rho = 5, 1e-2
    lf = lambda S, Smax: 1 / (1 + np.exp(-(S - (Smax - rho * Smax * eps)) / (rho * Smax)))
    
    Lzm = S2_FPmax + S2_FSmax + S2_Tmax
    Uztw, Uzfw, Lztw, Lzps, Lzfs = y[:5]
    Sf = y[5:8]
    Lztot = Lztw + Lzps + Lzfs
    
    e_1 = Ep * min(Uztw, S1_Tmax) / S1_Tmax
    e_2 = (Ep - e_1) * min(Lztw, S2_Tmax) / S2_Tmax
    if e_2 < 0:
        e_2 = 0
    
    q_0 = nu_p * S2_FPmax + nu_s * S2_FSmax
    dlz = 1 + alfa * (Lztot / Lzm) ** psi
    q_12 = q_0 * dlz * (Uzfw / S1_Fmax)
    q_if = ki * (Uzfw / S1_Fmax)
    q_bp = nu_p * Lzps
    q_bs = nu_s * Lzfs
    ac = Uztw / S1_Tmax * a_cmax
    q_sx = ac * P
    q_utof = (P - q_sx) * lf(Uztw, S1_Tmax)
    q_ufof = q_utof * lf(Uzfw, S1_Fmax)
    
    q_stof = kappa * q_12 * lf(Lztw, S2_Tmax)
    prc_s = (1 - kappa) * q_12 / 2 + q_stof / 2
    q_sfofp = prc_s * lf(Lzps, S2_FPmax)
    q_sfofs = prc_s * lf(Lzfs, S2_FSmax)
    
    q_fout = kf * Sf
    
    dydt = np.zeros(9)
    dydt[0] = P - q_sx - e_1 - q_utof
    dydt[1] = q_utof - q_12 - q_if - q_ufof
    dydt[2] = kappa * q_12 - e_2 - q_stof
    dydt[3] = prc_s - q_bp - q_sfofp
    dydt[4] = prc_s - q_bs - q_sfofs
    dydt[5] = q_if + q_sx + q_ufof + q_sfofp + q_sfofs - q_fout[0]
    dydt[6] = q_fout[0] - q_fout[1]
    dydt[7] = q_fout[1] - q_fout[2]
    dydt[8] = q_bp + q_bs + q_fout[2]
    
    return dydt


def sacsma_odefcn(t, x, S1_Fmax, S1_Tmax, S2_FPmax, S2_FSmax, S2_Tmax, alfa, psi, ki, kappa, nu_p, nu_s, a_cmax, kf, plugin):
    # Define return argument (dxdt as a numpy array)
    dxdt = np.full(9, np.nan)  # Initialize dxdt array with NaN values
    
    # Define smoothing function
    eps = 5
    rho = 0.01  # Dimensionless smoothing coefficients
    
    # Correct implementation based on FUSE source code (GitHub)
    lf = lambda S, Smax: 1 / (1 + np.exp(-(S - (Smax - rho * Smax * eps)) / (rho * Smax)))
    
    # Maximum storage of lower zone (mm)
    Lzm = S2_FPmax + S2_FSmax + S2_Tmax
    
    # Unpack states and forcing
    Uztw = x[0]                                         # Tension water storage upper zone (mm)
    Uzfw = x[1]                                         # Free water storage upper zone (mm)
    Lztw = x[2]                                         # Tension water storage lower zone (mm)
    Lzps = x[3]                                         # Free water storage of lower zone primary (mm)
    Lzfs = x[4]                                         # Free water storage of lower zone secondary (mm)
    Sf = x[5:8]                                         # Storage of 1st, 2nd, and 3rd fast reservoirs (routing, mm)
    Lztot = Lztw + Lzps + Lzfs                          # Total lower zone storage (mm)
    
    it = int(np.floor(t))                               # Truncate time to current time index
    P = plugin['data']['P'][it]                         # Current rainfall (mm/d)
    Ep = plugin['data']['Ep'][it]                       # Current Ep (mm/d)
    
    # Compute fluxes between layers
    e_1 = Ep * min(Uztw, S1_Tmax) / S1_Tmax             # Evaporation from upper soil layer (mm/d)
    e_2 = (Ep - e_1) * min(Lztw, S2_Tmax) / S2_Tmax     # Evaporation from lower soil layer (mm/d)
    q_0 = nu_p * S2_FPmax + nu_s * S2_FSmax             # Baseflow at saturation (mm/d)
    # Handle negative values (= futile but important)
    if (Lztot / Lzm) >= 0:
        dlz = 1 + alfa * (Lztot / Lzm) ** psi           # Lower-zone percolation demand (-)
    else:
        dlz = 1

    q_12 = q_0 * dlz * (Uzfw / S1_Fmax)                 # Percolation from UZ to LZ (mm/d)
    q_if = ki * (Uzfw / S1_Fmax)                        # Interflow (mm/d)
    q_bp = nu_p * Lzps                                  # Primary baseflow (mm/d)
    q_bs = nu_s * Lzfs                                  # Secondary baseflow (mm/d)
    ac = Uztw / S1_Tmax * a_cmax                        # Saturated area (-)
    q_sx = ac * P                                       # Saturation excess (mm/d)
    q_utof = (P - q_sx) * lf(Uztw, S1_Tmax)             # Overflow of water from tension storage in upper soil layer (mm/d)
    q_ufof = q_utof * lf(Uzfw, S1_Fmax)                 # Overflow of water from free storage in the upper soil layer (mm/d)
    
    # Ensure that residual evaporation cannot be negative
    e_2 = max(e_2, 0)
    
    q_stof = kappa * q_12 * lf(Lztw, S2_Tmax)           # Overflow of water from tension storage in the lower soil layer (mm/d)
    prc_s = (1 - kappa) * q_12 / 2 + q_stof / 2         # Percolation and q_stof (mm/d)
    q_sfofp = prc_s * lf(Lzps, S2_FPmax)                # Overflow of water from primary base flow storage in lower soil layer (mm/d)
    q_sfofs = prc_s * lf(Lzfs, S2_FSmax)                # Overflow of water from secondary base flow storage in lower soil layer (mm/d)
    
    # Compute fluxes out of fast reservoirs
    q_fout = kf * Sf                                    # Channel inflow is q_fout[3]
    
    # Compute net flux into/out of reservoir
    dxdt[0] = P - q_sx - e_1 - q_utof                   # Net flux into/out of S1T (mm/d)
    dxdt[1] = q_utof - q_12 - q_if - q_ufof             # Net flux into/out of S1F (mm/d)
    dxdt[2] = kappa * q_12 - e_2 - q_stof               # Net flux into/out of S2T (mm/d)
    dxdt[3] = prc_s - q_bp - q_sfofp                    # Net flux into/out of S2Fa (mm/d)
    dxdt[4] = prc_s - q_bs - q_sfofs                    # Net flux into/out of S2Fb (mm/d)
    dxdt[5] = q_if + q_sx + q_ufof + q_sfofp + q_sfofs - q_fout[0]  # Net flux into/out of fast reservoir 1 (mm/d)
    dxdt[6] = q_fout[0] - q_fout[1]                     # Net flux into/out of fast reservoir 2 (mm/d)
    dxdt[7] = q_fout[1] - q_fout[2]                     # Net flux into/out of fast reservoir 3 (mm/d)
    dxdt[8] = q_bp + q_bs + q_fout[2]                   # Inflow (mm/d) to infinite reservoir (streamflow from delta)
    
    return dxdt
