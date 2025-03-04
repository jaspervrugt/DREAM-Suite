import numpy as np

############################
### Case study 35
############################
def Haverkamp_t(eta, plugin):
    """
    This function solves for the t(I) relationship of Haverkamp.
    
    Parameters:
    eta (array-like): 3x1 vector with S [cm/h^0.5], Ks [cm/h], and beta [-]
    plugin (dict): Dictionary with 'I' as input (cumulative infiltration values in cm)
    
    Returns:
    tuple: (t, dtdI, flag)
        - t (numpy.ndarray): time, t, in hours corresponding to I (cm)
        - dtdI (numpy.ndarray): derivative of time with respect to I, in cm
        - flag (int): exit flag: [1] finite, [2] infinite time values
    """
    # Unpack eta and plugin values
    S = eta[0]
    Ks = eta[1]
    Ki = 0                      # Ki is set to zero as per the original function
    B = eta[2] if len(eta) >= 3 else plugin['B']
    
    # Compute dK and xi
    dK = Ks - Ki
    xi = dK / S**2
    I = np.array(plugin['I'])   # Cumulative infiltration
    n = len(I)
    
    t = np.full(n, np.nan)      # Initialize time array
    flag = 1                    # Initialize flag
    
    # Calculate gamma (expression for later use)
    gamma = 2 * B * xi * I
    
    if I[0] == 0:
        t[0] = 0
        j = 1
    else:
        j = 0
    
    # Dynamic part: Evaluate time expression
    for m in range(j, n):
        if gamma[m] < 709.783:  # Avoid numerical overflow
            t[m] = 1 / (dK * (B - 1) * xi) * (1 / 2 * np.log(np.exp(gamma[m]) / B + (B - 1) / B) - xi * I[m])
        else:
            flag = 2
            break
    
    # Calculate the derivative of time with respect to I (dtdI)
    dtdI = (np.exp((2 * B * I * Ks) / S**2)) / ((B - 1) * (np.exp((2 * B * I * Ks) / S**2) / B + (B - 1) / B)) - 1 / (Ks * (B - 1))
    
    return t, dtdI, flag
