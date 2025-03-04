import numpy as np
from Haverkamp_I import Haverkamp_I
from Haverkamp_t import Haverkamp_t

############################
### Case study 35
############################
def Haverkamp_It(eta, plugin):
    """
    Evaluates Haverkamp_I and Haverkamp_t and combines their infiltration times
    and corresponding cumulative infiltration values.
    
    Parameters:
    eta (array-like): 4x1 vector with S [cm/h^0.5], Ks [cm/h], beta [-], Ki [cm/h]
    plugin (dict): A dictionary containing input variables for Haverkamp model.
    
    Returns:
    numpy.ndarray: 2nx1 vector of cumulative infiltration (I, cm) and time (t, h).
    """
    # First, calculate cumulative infiltration I using Haverkamp_I function
    I, _, _ = Haverkamp_I(eta, plugin)     # Assume Haverkamp_I returns (I, i, flag)
    
    # Then, calculate time using Haverkamp_t function (this function must exist in your code)
    t = Haverkamp_t(eta, plugin)[0]        # Assume Haverkamp_t returns time as an array
    
    # Combine time and cumulative infiltration into one array
    y = np.concatenate((t, I))             # Stack t and I vertically to form a 2nx1 array

    return y