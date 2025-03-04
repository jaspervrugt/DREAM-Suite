import numpy as np
from scipy.interpolate import interp1d

############################
### Case study 34
############################
def Haverkamp_I_patch(eta, plugin, rtol = 1e-12, kmax = 20, itol = 1e-10):
    """
    Solve for the I(t) relationship of Haverkamp using Newton's method.
    
    Parameters:
    eta : array-like (4,)
        Parameter values for Haverkamp model:
        - S: [cm/h^0.5]
        - Ks: [cm/h]
        - ß: [-]
        - Ki: [cm/h]
    plugin : object
        Structure with time data and model setup
    rtol : float, optional (default=1e-12)
        Tolerance on function value root
    kmax : int, optional (default=20)
        Maximum number of Newton iterations
    itol : float, optional (default=1e-10)
        Tolerance constant rate assumption
    
    Returns:
    I : array
        Cumulative infiltration (cm)
    i : array
        Infiltration rate (cm/h)
    flag : int
        Exit flag:
        - 1: Exact solution
        - 2: Approximate solution
    """
    
    # Unpack parameter values
    S = eta[0]          # cm/h^(1/2)
    Ks = eta[1]         # cm/h
    if plugin['model_setup'] == 1:
        B = eta[2]      # Unitless
        Ki = 0          # cm/h
    elif plugin['model_setup'] == 2:
        Ki = eta[2]     # cm/h
        B = eta[3]      # Unitless

    # Initialization
    dK = Ks - Ki
    xi = dK / S**2
    t = np.array(plugin['t']).flatten()             # time vector
    n = len(t)
    I = np.full(n, np.nan)                          # Cumulative infiltration
    i = np.full(n, np.nan)                          # Infiltration rate
    flag = 1                                        # Initialize flag

    # Set I(1) to zero if t[0] is zero
    if t[0] == 0:
        I[0] = 0
        j = 1
    else:
        j = 0
    # j = 1 if t[0] != 0 else 2
    dt = np.concatenate(([t[0]], np.diff(t)))       # Time step differences

    # Residual and derivative functions
    r = lambda I, t: I - 1 / (2 * xi) * np.log(np.exp(2 * xi * B * (I - Ki * t)) / B + (B - 1) / B) - (dK * (1 - B) + Ki) * t
    dr = lambda I, t: (1 - (B * np.exp(2 * xi * B * (I - Ki * t))) / (np.exp(2 * xi * B * (I - Ki * t)) + B - 1))

    # Dynamic part
    I_up = S * np.sqrt(t) + t * Ks
    I_low = np.maximum(I_up[0] / 2, 0.1)
    gamma = 2 * xi * B * (I_up - Ki * t)

    for m in range(j, n):
        if gamma[m] < 709.783:
            y = [I_low]
            k = 0
            while abs(r(y[k], t[m])) > rtol and k < kmax:
                y.append(y[k] - r(y[k], t[m]) / dr(y[k], t[m]))  # Next iterate
                k += 1

            I[m] = y[k]     # Cumulative infiltration at time t[m]
            I_low = y[k]    # Update lower bound for next iteration
            # Calculate infiltration rate
            A = np.exp(2 * B * xi * (I[m] - Ki * t[m]))
            i[m] = (((1 - B) * Ks + B * Ki) * (A + B - 1) - A * B * Ki) / (A + B - A * B - 1)
        elif m > 0:
            if abs(i[m - 1] - Ks) < itol:
                i[m] = Ks
                I[m] = I[m - 1] + dt[m] * i[m]  # Constant infiltration rate
            else:
                flag = 2
                break       # Approximate solution
        else:
            flag = 2
            break           # Approximate solution

    return I, i, flag
