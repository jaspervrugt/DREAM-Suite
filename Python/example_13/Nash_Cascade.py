import numpy as np
import scipy.special as sp  

############################
### Case study 13
############################
def Nash_Cascade(k):
    # Nash-Cascade unit hydrograph -- series of three linear reservoirs
    
    if not hasattr(Nash_Cascade, "initialized"):                # Store local variables in memory
        daily_data = np.loadtxt('03451500.dly')                 # Load French Broad data
        Nash_Cascade.maxT = 365                                 # Maximum time       
        Nash_Cascade.P = daily_data[0:Nash_Cascade.maxT, 3]     # Precipitation
        Nash_Cascade.F = 767 * (1000) / (60 * 60 * 24)          # Area factor: from mm/d to m3/s
        Nash_Cascade.n = 3                                      # Define number of linear reservoirs
        Nash_Cascade.Time = np.arange(1, Nash_Cascade.maxT + 1) # Define Time
        Nash_Cascade.initialized = True                         # Flag to indicate that initialization is complete        
    
    # ------------
    # Model script
    # ------------
    if k < 1:                                                   # Write to screen
        print('Nash_Cascade: Recession constant < 1 day --> numerical errors possible')

    A = np.zeros((Nash_Cascade.maxT, Nash_Cascade.maxT))            # Define help matrix
    IUH = 1 / (k * sp.gamma(Nash_Cascade.n)) * (Nash_Cascade.Time / k) \
        ** (Nash_Cascade.n - 1) * np.exp(-Nash_Cascade.Time / k)    # Instantaneous unit hydrograph
    for t in range(Nash_Cascade.maxT):                              # Loop over time
        id = np.arange(0, Nash_Cascade.maxT - t)                    # Define id
        A[t, t:Nash_Cascade.maxT] = Nash_Cascade.P[t] * IUH[id]     # Calculate flow

    sim_Q = Nash_Cascade.F * np.sum(A, axis = 0)                    # Now determine total flow
    
    return sim_Q

