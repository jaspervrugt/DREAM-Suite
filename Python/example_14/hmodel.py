import numpy as np
from scipy.integrate import odeint   
from scipy.integrate import solve_ivp

############################
### Case study 6, 9, 14, 28
############################
def hmodel(x, plugin, model_code = 3):
    # ####################################################################### #
    # HMODEL Different numerical solutions of HMODEL of Schoups & Vrugt, 2010 #
    # Formulation/language                                                    #
    #   1: Python Runge Kutta implementation of hmodel_ode                    #
    #   2: Python ode45 implementation of hmodel_ode                          #
    #   3: Python MATLAB Explicit Euler hmodel_ode with int_steps steps       #
    #   4: 'C++': Runge Kutta implementation hmodel_ode in C++                #
    # ####################################################################### #

    nv = 5                          # Initialize number of state variables
    ns = len(plugin['tout'])        # Number of time steps
    Y = np.full((ns, nv), np.nan)   # Initialize matrix of state variables
    Y[0, :] = plugin['y0'].copy()   # Initialize state variables at time 0

    # Initialize model parameters
    if model_code in [1, 2, 3]:
        I_max = x[0]                # Interception storage capacity (mm)
        Su_max = x[1]               # Unsaturated zone storage capacity (mm)
        Qs_max = x[2]               # Maximum percolation rate (mm/d)
        aE = x[3]                   # Evaporation coefficient
        aF = x[4]                   # Runoff coefficient
        aS = 1e-6                   # Percolation coefficient
        Kf = x[5]                   # Fast-flow response time (d)
        Ks = x[6]                   # Slow-flow response time (d)
    elif model_code == 5:
        data = {'P': plugin['data']['P'].T,
                'Ep': plugin['data']['Ep'].T,
                'I_max': x[0],
                'Su_max': x[1],
                'Qs_max': x[2],
                'aE': x[3],
                'aF': x[4],
                'aS': 1e-6,
                'Kf': x[5],
                'Ks': x[6]}
    
    # Execute model based on model_code
    if model_code == 1:             # MATLAB Runge Kutta implementation
        hin = plugin['options']['InitialStep']
        hmax_ = plugin['options']['MaxStep']
        hmin_ = plugin['options']['MinStep']
        reltol = plugin['options']['RelTol']
        abstol = plugin['options']['AbsTol']
        order = plugin['options']['Order']
        
        for s in range(1, ns):
            t1 = plugin['tout'][s-1]
            t2 = plugin['tout'][s]
            h = hin
            h = max(hmin_, min(h, hmax_))
            h = min(h, t2 - t1)
            Y[s, :] = Y[s-1, :].copy()
            t = t1

            while t < t2:
                ytmp = Y[s, :]
                ytmp, LTE = rk2_hmodel(ytmp, h, I_max, Su_max, Qs_max, aE, aF, aS, Kf, Ks,
                                plugin['data']['P'][s-1], plugin['data']['Ep'][s-1])
                w = 1 / (reltol * np.abs(ytmp) + abstol)
                wrms = np.sqrt(np.sum((w * LTE) ** 2) / nv)
                if h <= hmin_:
                    wrms = 0.5
                if wrms <= 1 or h <= hmin_:
                    Y[s, :] = ytmp.copy()
                    t += h
                h = h * max(0.2, min(5.0, 0.9 * wrms ** (-1 / order)))
                h = max(hmin_, min(h, hmax_))
                h = min(h, t2 - t)

    elif model_code == 2:  # ode45 implementation of hmodel_odefcn
        plugin['tout'][-1] -= 1e-10
        ode_options = { 'rtol': plugin['options']['RelTol'],        # Relative tolerance
                        'atol': plugin['options']['AbsTol'],        # Absolute tolerance
                        'max_step': plugin['options']['MaxStep'],   # Maximum step size
                        'min_step': plugin['options']['MinStep']}   # Minimum step size (optional)
        #t, Y = odeint(hmodel_odefcn, plugin['y0'], plugin['tout'], args=(I_max, Su_max, Qs_max, aE, aF, aS, Kf, Ks, plugin), tfirst=True)
        def ode_func(t, y):
            return hmodel_odefcn(t, y, I_max, Su_max, Qs_max, aE, aF, aS, Kf, Ks, plugin)

        solution = solve_ivp(ode_func, (plugin['tout'][0], plugin['tout'][-1]), plugin['y0'], t_eval=plugin['tout'], **ode_options)
        Y = solution.y.T

    elif model_code == 3:  # Explicit Euler implementation
        int_steps = 10
        dt = 1 / int_steps
        for s in range(1, ns):      # MATLAB: s = 2,...,ns    Python: s = 1,...,ns-1
            y = Y[s-1, :].copy()
            for it in range(int_steps):
                dydt = hmodel_odefcn(s-1, y, I_max, Su_max, Qs_max, aE, aF, aS, Kf, Ks, plugin)
                y += dydt * dt
            Y[s, :] = y.copy()

    # elif model_code == 5:  # C++ implementation (with dll)
    #     Y = crr_hmodel(plugin['tout'], plugin['y0'], data, plugin['options']).T
    SimRR = np.diff(Y[plugin['idx'], nv-1])  # Extract discharge

    return SimRR


# Runge-Kutta integration function
def rk2_hmodel(y, h, I_max, Su_max, Qs_max, aE, aF, aS, Kf, Ks, P, Ep):
    
    dydtE = fRhs_hmodel(y, I_max, Su_max, Qs_max, aE, aF, aS, Kf, Ks, P, Ep)
    yE = y + h * dydtE                  # Euler solution
    dydtH = fRhs_hmodel(yE, I_max, Su_max, Qs_max, aE, aF, aS, Kf, Ks, P, Ep)       # Evaluate dx
    y = y + 0.5 * h * (dydtE + dydtH)   # Heun solution
    LTE = np.abs(yE - y)                # Compute estimate of LTE
    return y, LTE


# Right-hand-side for hmodel function
def fRhs_hmodel(y, I_max, Su_max, Qs_max, aE, aF, aS, Kf, Ks, P, Ep):
    
    Si = y[0]  # Interception storage (mm)
    Su = y[1]  # Surface storage (mm)
    Sf = y[2]  # Storage of fast/quick reservoirs (mm)
    Ss = y[3]  # Storage of slow reservoir (mm)
    Sur = min(Su / Su_max, 1)
    
    if I_max > 0:
        Sir = Si / I_max
        EvapI = Ep * expFlux(Sir, 50)
        P_e = P * expFlux(Sir, -50)
        Ep_e = max(0., Ep - EvapI)
    else:
        EvapI = 0
        P_e = P
        Ep_e = Ep

    Ea = Ep_e * expFlux(Sur, aE)
    prc = Qs_max * expFlux(Sur, aS)
    rnf = P_e * expFlux(Sur, aF)
    qf = Sf / Kf
    qs = Ss / Ks

    dydt = np.zeros(5)
    dydt[0] = P - EvapI - P_e
    dydt[1] = P_e - Ea - prc - rnf
    dydt[2] = rnf - qf
    dydt[3] = prc - qs
    dydt[4] = qf + qs
    return dydt

def expFlux(Sr, a):
    Sr = max(0.0, min(1.0, Sr))
    if abs(a) < 1e-6:
        return Sr  # Approximately linear
    else:
        return (1. - np.exp(-a * Sr)) / (1. - np.exp(-a))

def hmodel_odefcn(t, x, I_max, Su_max, Qs_max, aE, aF, aS, Kf, Ks, plugin):
    dxdt = np.full(5, np.nan)
    Si = x[0]
    Su = x[1]
    Sf = x[2]
    Ss = x[3]
    it = int(np.floor(t))
    P = plugin['data']['P'][it]
    Ep = plugin['data']['Ep'][it]
    Sur = Su / Su_max
    
    if I_max > 0:
        Sir = Si / I_max
        EvapI = Ep * expFlux(Sir, 50)
        P_e = P * expFlux(Sir, -50)
        Ep_e = max(0., Ep - EvapI)
    else:
        EvapI = 0
        P_e = P
        Ep_e = Ep

    Ea = Ep_e * expFlux(Sur, aE)
    prc = Qs_max * expFlux(Sur, aS)
    rnf = P_e * expFlux(Sur, aF)
    qf = Sf / Kf
    qs = Ss / Ks

    dxdt[0] = P - EvapI - P_e
    dxdt[1] = P_e - Ea - prc - rnf
    dxdt[2] = rnf - qf
    dxdt[3] = prc - qs
    dxdt[4] = qf + qs
    return dxdt


# Exponential flux function
def exponen(x):
    return np.exp(np.minimum(300, x))
