import numpy as np
import scipy.special as sp


############################
### Case study 22
############################
def Nash_Cascade_2(x):
    # ####################################################################### #
    # NASH_CASCADE Nash-Cascade unit hydrograph: serie of 3 linear reservoirs #
    # ####################################################################### #
    
    if not hasattr(Nash_Cascade_2, "initialized"):                      # Store local variables in memory
        Nash_Cascade_2.maxT = 25                                        # Maximum time       
        Nash_Cascade_2.P =  [ 10, 25, 8, 2, 0, 0, 0, 0, 0, 0, 0, 0, \
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]                      # Precipitation
        Nash_Cascade_2.Time = np.arange(1, Nash_Cascade_2.maxT + 1)     # Define Time
        Nash_Cascade_2.initialized = True                               # Flag to indicate that initialization is complete        
    
    # ------------
    # Model script
    # ------------
    k = x[0]
    n = x[1]
    if k < 1:                                                           # Write to screen
        print('Nash_Cascade_2: Recession constant < 1 day --> numerical errors possible')

    A = np.zeros((Nash_Cascade_2.maxT, Nash_Cascade_2.maxT))            # Define help matrix
    IUH = 1 / (k * sp.gamma(n)) * (Nash_Cascade_2.Time / k) \
        ** (n - 1) * np.exp(-Nash_Cascade_2.Time / k)                   # Instantaneous unit hydrograph
    for t in range(Nash_Cascade_2.maxT):                                # Loop over time
        id = np.arange(0, Nash_Cascade_2.maxT - t)                      # Define id
        A[t, t:Nash_Cascade_2.maxT] = Nash_Cascade_2.P[t] * IUH[id]     # Calculate flow

    sim_Q = np.sum(A, axis = 0)                                         # Now determine total flow
    
    return sim_Q
