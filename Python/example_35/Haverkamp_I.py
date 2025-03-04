import numpy as np

############################
### Case study 35
############################
def Haverkamp_I(eta, plugin, rtol = 1e-12, kmax = 20):
    """
    This function solves for the I(t) relationship of Haverkamp using Newton's method.
    
    Parameters:
    eta (array-like): 4x1 vector with S [cm/h^.5], Ks [cm/h], ß [-], Ki [cm/h]
    plugin (dict): A dictionary with 't' (time vector) as a key, where 't' is a numpy array of time values.
    rtol (float, optional): Tolerance on function value at root (default: 1e-12).
    kmax (int, optional): Maximum number of Newton iterations (default: 20).
    
    Returns:
    I (numpy array): Cumulative infiltration (cm) as a function of time.
    i (numpy array): Infiltration rate (cm/h) as a function of time.
    flag (int): Exit flag: [1] exact solution, [2] approximate solution.
    """

    # Unpack parameters
    if len(eta) < 3:
        eta = np.append(eta, plugin['B'])
    
    if len(eta) < 4:
        eta = np.append(eta, 0)

    S, Ks, B, Ki = eta
    dK = Ks - Ki
    xi = dK / S**2
    t = np.asarray(plugin['t'])  # Ensure t is a numpy array
    n = len(t)
    
    # Initialize variables
    I = np.full(n, np.inf)
    flag = 1  # Assume success (exact solution)
    
    if dK <= 0:
        return np.inf * np.ones(n), np.inf * np.ones(n), 0  # If dK <= 0, return inf and flag as 0
    
    # Initial condition for I(1)
    if t[0] == 0:
        I[0] = 0
        j = 1
    else:
        j = 0
    
    # Define the residual and derivative functions
    def r_func(I, t):
        return I - (1 / (2 * xi)) * np.log(np.exp(2 * xi * B * (I - Ki * t)) / B + (B - 1) / B) - (dK * (1 - B) + Ki) * t

    def dr_func(I, t):
        exp_term = np.exp(2 * xi * B * (I - Ki * t))
        return 1 - (B * exp_term) / (exp_term + B - 1)

    # Dynamic part of the model
    I_up = S * np.sqrt(t) + t * Ks
    I_low = np.maximum(I_up[0] / 2, 0.1)
    gamma = 2 * xi * B * (I_up - Ki * t)
        
    for m in range(j, n):       # Iterate over each time step
        if gamma[m] < 709.783:  # Check overflow condition
            y = [I_low]         # Initial guess for the root
            k = 0               # Iteration counter
            # Newton's method iteration
            while abs(r_func(y[k], t[m])) > rtol and k < kmax:
                y.append(y[k] - r_func(y[k], t[m]) / dr_func(y[k], t[m]))
                k += 1
            
            # Assign the root (I(m))
            I[m] = I_low = y[k]
        else:
            flag = 2            # Set flag to 2 if overflow occurs
            break

    # Compute alpha (A) for the infiltration rate
    A = np.exp(2 * B * xi * (I - Ki * t))
    
    # Compute infiltration rate (i)
    i = (((1 - B) * Ks + B * Ki) * (A + B - 1) - A * B * Ki) / (A + B - A * B - 1)

    return I, i, flag
