import numpy as np
from scipy.integrate import odeint    

#from scipy.linalg import inv, det

############################
### Case study 4, 17, 32
############################
def hymod(x, plugin, model_code = 1):
    # ####################################################################### #
    # HYMOD Different numerical solutions of Hydrologic Model of Boyle (2001) #
    # Formulation/language                                                    #  
    #   1: Python Runge Kutta implementation of hymod_ode                     #  
    #   2: Python ode45 implementation of hymod_ode                           #
    #   3: Python MATLAB Explicit Euler hymod_ode with int_steps steps        #
    #   4: Python old implementation                                          #
    #   5: 'C++': Runge Kutta implementation hymod_ode in C++                 #
    # ####################################################################### #

    nv = 6                              # number of state variables
    ns = len(plugin['tout'])            # number of time steps
    Y = np.full((ns, nv), np.nan)       # initialize matrix of state variables
    Y[0, :nv] = np.array(plugin['y0'])  # initialize state variables at time 0

    # Initialization
    if model_code in [1, 2, 3, 4]:
        Su_max = x[0]                   # Maximum storage of unsaturated zone
        beta_ = x[1]                    # Spatial variability of soil moisture capacity
        alfa = x[2]                     # Flow partitioning coefficient
        Ks = x[3]                       # Recession constant slow reservoir
        Kf = x[4]                       # Recession constant fast reservoir
    elif model_code == 5:
        data = {
            'P': plugin['data']['P'],
            'Ep': plugin['data']['Ep'],
            'm': 1e-2,                  # smoothing parameter
            'Sumax': x[0],              # Maximum storage of unsaturated zone
            'beta': x[1],               # Spatial variability of soil moisture capacity
            'alfa': x[2],               # Flow partitioning coefficient
            'Ks': x[3],                 # Recession constant slow reservoir
            'Kf': x[4]                  # Recession constant fast reservoir
        }

    # Execute model
    if model_code == 1:  # Runge Kutta implementation
        hin = plugin['options']['InitialStep']
        hmax_ = plugin['options']['MaxStep']
        hmin_ = plugin['options']['MinStep']
        reltol = plugin['options']['RelTol']
        abstol = plugin['options']['AbsTol']
        order = plugin['options']['Order']
        
        for s in range(1, ns):      ## Integrate from tprint(1) to tprint(end)
            t1 = plugin['tout'][s - 1]
            t2 = plugin['tout'][s]
            h = hin
            h = max(hmin_, min(h, hmax_))
            h = min(h, t2 - t1)
            Y[s, :nv] = Y[s - 1, :nv].copy()
            t = t1
            while t < t2:
                ytmp = Y[s, :nv]    # Initial state at t
                ytmp, LTE = rk2_hymod(ytmp, h, Su_max, beta_, alfa, Ks, Kf, plugin['data']['P'][s - 1], plugin['data']['Ep'][s - 1])
                w = 1. / (reltol * np.abs(ytmp) + abstol)
                wrms = np.sqrt(np.dot(w, LTE) / nv)
                if h <= hmin_:
                    wrms = 0.5
                if wrms <= 1 or h <= hmin_:
                    Y[s, :nv] = ytmp.copy()
                    t += h
                h = h * max(0.2, min(5.0, 0.9 * wrms ** (-1 / order)))
                h = max(hmin_, min(h, hmax_))
                h = min(h, t2 - t)
    elif model_code == 2:           ## ode45 implementation
        plugin['tout'][-1] -= 1e-10
        ode_options = { 'rtol': plugin['options']['RelTol'],        # Relative tolerance
                        'atol': plugin['options']['AbsTol'],        # Absolute tolerance
                        'max_step': plugin['options']['MaxStep'],   # Maximum step size
                        'min_step': plugin['options']['MinStep']}   # Minimum step size (optional)
        # t, Y = odeint(hymod_odefcn, plugin['y0'], plugin['tout'], args=(Su_max, beta_, alfa, Ks, Kf, plugin), tfirst=True)
        
        def ode_func(t, y):
            return hymod_odefcn(t, y, Su_max, beta_, alfa, Ks, Kf, plugin, plugin)

        solution = solve_ivp(ode_func, (plugin['tout'][0], plugin['tout'][-1]), plugin['y0'], t_eval=plugin['tout'], **ode_options)
        Y = solution.y.T

    elif model_code == 3:  ## Explicit Euler implementation
        int_steps = 10
        dt = 1 / int_steps
        for s in range(1, ns):      # MATLAB: s = 2,...,ns    Python: s = 1,...,ns-1
            y = Y[s - 1, :nv].copy()
            for it in range(int_steps):
                dydt = hymod_odefcn(s - 1, y, Su_max, beta_, alfa, Ks, Kf, plugin)
                y += dydt * dt
            Y[s, :nv] = y.copy()

    elif model_code == 4:  ## Old implementation
        bexp = beta_
        for s in range(1, ns):
            y = Y[s - 1, :nv].copy()
            P_e, y[0] = excess(y[0], Su_max, bexp, plugin['data']['P'][s - 1], plugin['data']['Ep'][s - 1])
            qf_in = alfa * P_e
            qs_in = (1 - alfa) * P_e
            y[1], qs_out = linres(y[1], qs_in, Ks)
            for k in range(2, 5):
                y[k], qf_out = linres(y[k], qf_in, Kf)
                qf_in = qf_out
            y[5] = y[5] + qf_out + qs_out
            Y[s, :nv] = y.copy()
    # elif model_code == 5:  ## C++ implementation
    #     Y = crr_hymod(plugin['tout'], plugin['y0'], data, plugin['options']).T

    idx = np.array(plugin['idx']).astype(int)
    SimRR = np.diff(Y[idx, nv-1])                       # extract discharge in mm/day   
    # Be careful with indices! Fixed in input file [re: MATLAB/Python]i

    return SimRR


# Runge-Kutta integration function
def rk2_hymod(y, h, Su_max, beta_, alfa, Ks, Kf, P, Ep):

    dydtE = fRhs_hymod(y, Su_max, beta_, alfa, Ks, Kf, P, Ep)
    yE = y + h * dydtE                  # Euler solution
    dydtH = fRhs_hymod(yE, Su_max, beta_, alfa, Ks, Kf, P, Ep)
    y = y + 0.5 * h * (dydtE + dydtH)   # Heun solution
    LTE = np.abs(yE - y)                # Integration error
    return y, LTE


# Right-hand side function for the HYMOD model
def fRhs_hymod(y, Su_max, beta_, alfa, Ks, Kf, P, Ep):

    m = 1e-2  # smoothing parameter
    Su = y[0]
    Ss = y[1]
    Sf = y[2:5]
    Sur = min(Su / Su_max, 1)
    Perc = P * (1 - (1 - Sur) ** beta_)
    Ea = Ep * Sur * (1 + m) / (Sur + m)
    dydt = np.zeros(6)
    dydt[0] = P - Perc - Ea
    qs_in = (1 - alfa) * Perc
    qs_out = Ks * Ss
    dydt[1] = qs_in - qs_out
    qf_in = alfa * Perc
    qf_out = Kf * Sf
    dydt[2] = qf_in - qf_out[0]
    dydt[3] = qf_out[0] - qf_out[1]
    dydt[4] = qf_out[1] - qf_out[2]
    dydt[5] = qs_out + qf_out[2]
    return dydt


def hymod_odefcn(t, x, Su_max, beta_, alfa, Ks, Kf, plugin):

    dxdt = np.nan * np.ones(6)              # Initialize return argument
    m = 1e-2                                # Smoothing coefficient (mm)
    Su = x[0]                               # Surface storage (mm)
    Ss = x[1]                               # Storage of slow reservoir (mm)
    Sf = x[2:5]                             # Storage of fast/quick reservoirs (mm)
    it = int(np.floor(t))                   # Truncate time to current time index
    P = plugin['data']['P'][it]             # Get current rainfall (mm/d)
    Ep = plugin['data']['Ep'][it]           # Get current Ep (mm/d)
    Sur = min(Su / Su_max, 1)               # Ratio of surface storage to max storage
    qu = P * (1 - (1 - Sur) ** beta_)       # Precipitation converted to flow, qu (mm/d)
    Ea = Ep * Sur * (1 + m) / (Sur + m)     # Actual evaporation (mm/d)
    dxdt[0] = P - qu - Ea                   # Flux into/out of unsaturated reservoir (mm/d)
    qs = (1 - alfa) * qu                    # Inflow to slow reservoir (mm/d)
    qs_o = Ks * Ss                          # Flow out of slow reservoir (mm/d)
    dxdt[1] = qs - qs_o                     # Flux into/out of slow reservoir (mm/d)
    qf = alfa * qu                          # Inflow to first quick reservoir (mm/d)
    qf_o = Kf * Sf                          # Flow out of three quick reservoirs (mm/d)
    dxdt[2] = qf - qf_o[0]                  # Flux into/out of first fast reservoir (mm/d)
    dxdt[3] = qf_o[0] - qf_o[1]             # Flux into/out of 2nd fast reservoir (mm/d)
    dxdt[4] = qf_o[1] - qf_o[2]             # Flux into/out of 3rd fast reservoir (mm/d)
    dxdt[5] = qf_o[2] + qs_o                # Inflow (mm/d) to infinite reservoir
    
    return dxdt


def excess(Su, Su_max, bexp, P, Ep):
    
    Su_old = Su
    b1exp = bexp + 1
    prev = Su_max * (1 - (1 - b1exp * Su / Su_max) ** (1 / b1exp))
    
    # Calculate effective rainfall 1
    P_e1 = max(P - Su_max + prev, 0)
    P -= P_e1
    dummy = min((prev + P) / Su_max, 1)
    Su_n = Su_max / b1exp * (1 - (1 - dummy) ** b1exp)
    # Calculate effective rainfall 2
    P_e2 = max(P - (Su_n - Su_old), 0)
    # Alternative approach: ET linearly related to soil moisture storage
    evap = Ep * (1 - (Su_max / b1exp - Su_n) / (Su_max / b1exp))
    Su_n = max(Su_n - evap, 0)  # Update state
    # Compute effective rainfall
    P_e = P_e1 + P_e2
    
    return P_e, Su_n


def linres(S, q_in, K):
    # Linear reservoir
    
    # S: state of reservoir (L)
    # q_in: inflow (L/T)
    # K: recession constant (1/L)
    
    # Unfamiliar implementation: Implicit for fixed time step of 1 unit
    Sn = (1 - K) * S + (1 - K) * q_in  	# New state of reservoir (L)
    q_out = (K / (1 - K)) * Sn  	# Flux out of reservoir (L)
    
    # Conventional method: instantaneous flux integrated over time units
    # q_out = K * S  # Flux out of reservoir (L/T)
    # Sn = S + (q_in - q_out) * itu  	# New state of reservoir (L)
    
    return Sn, q_out