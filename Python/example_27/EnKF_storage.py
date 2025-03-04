import numpy as np
from scipy.integrate import solve_ivp

############################
### Case study 27
############################
def EnKF_storage(par, plugin):
    # Vrugt, J. A., S.C. Dekker, and W. Bouten, Identification of rainfall 
    #   interception model parameters from measurements of throughfall and 
    #   forest canopy storage, Water Resources Research, 39(9), 1251,
    #   doi:10.1029/2003WR002013, 2003.   
    # Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten
    #   (2005), Improved treatment of uncertainty in hydrologic modeling: 
    #   Combining the strengths of global optimization and data assimilation, 
    #   Water Resources Research, 41, W01017, doi:10.1029/2004WR003059
    if not hasattr(EnKF_storage, "initialized"):    # Store local variables in memory
        EnKF_storage.N = plugin['N']                # ensemble size EnKF
        EnKF_storage.m = plugin['m']                # number of measurements at given time
        EnKF_storage.R = plugin['R']                # m x m measurement error covariance matrix
        EnKF_storage.C = plugin['C']                # m x m model error covariance matrix at given time
        EnKF_storage.obsTime = plugin['obsTime']    # nx1 vector of observation times
        EnKF_storage.Precip = plugin['Precip']      # nx1 vector of measured precipitation rates: mm/d
        EnKF_storage.Potevap = plugin['Potevap']    # nx1 vector of measured potential evapo rates: mm/d
        def InterceptionModelEquations(t, S, P_int, E0_int, a, b, c, d):
            """
            Interception model equations.
            Args:
            - t: Current time
            - S: Current storage
            - P_int: Interpolated rainfall value at time t
            - E0_int: Interpolated potential evapotranspiration value at time t
            - a, b, c, d: Model parameters
            Returns:
            - dSdt: Rate of change of storage (mm/day)
            """
            I = a * P_int                       # Interception in mm/day
            if S > c:                           # Calculate drainage (only if storage is larger than storage capacity)
                D = b * (S - c)                 # in mm/day
            else:
                D = 0

            E = d * E0_int * S / c              # Calculate evaporation in mm/day
            dSdt = I - D - E                    # Now calculate the change in storage in mm/day

            return dSdt
   
        EnKF_storage.IMEquations = InterceptionModelEquations
        EnKF_storage.initialized = True         # Flag to indicate that initialization is comple

    a, b, c, d, AR_rho = par                    # Extract parameters from 'par'
    
    N = EnKF_storage.N
    m = EnKF_storage.m
    R = EnKF_storage.R
    C = EnKF_storage.C
    if np.isscalar(C):
        C = np.array([[C]])
    T = EnKF_storage.obsTime
    P = EnKF_storage.Precip
    Ep = EnKF_storage.Potevap

    S = np.full((plugin['NT'], m), np.nan)                              # Initialize matrix of storages
    S[0, :] = plugin['Y'][0]                                            # Set the first mean simulated value equal to observed value
    Meas_err = np.random.multivariate_normal(np.zeros(m), R, N).T       # Draw m x N matrix of measurement errors from N_m(0, R)
    X_a = np.add(plugin['Y'][0], Meas_err)                              # m x N matrix of initial states of ensemble members: min of zero
    X_f = np.nan * np.zeros((m, N))                                     # m x N matrix of state forecasts of ensemble members
    Mod_err = np.random.multivariate_normal(np.zeros(m), C, N).T        # Draw m x N matrix with model errors
    for t in range(plugin['NT'] - 1):                                   # Now run the model one-observation ahead
        for ii in range(N):                                             # Generate state forecast for each ensemble member
            # Solve the ODE using solve_ivp
            result = solve_ivp(EnKF_storage.IMEquations, [T[t], \
                T[t + 1]], X_a[0:m, ii], \
                args = (P[t],Ep[t],a,b,c,d), \
                method = 'RK45', t_eval = [T[t], T[t+1]])               # Run the model from one time step to next for the nth state vector
            state_forecast = result.y[:, -1]                            # Get the state forecast at the final time
            X_f[0:m, ii] = state_forecast                               # Store the simulated state forecast at final time
            # X_f[:m, ii] = state_forecast[-1, :m]
            Mod_err[:m, ii] = AR_rho * Mod_err[:m, ii] + \
                np.sqrt(1 - AR_rho**2) * np.random.multivariate_normal(np.zeros(m), C)   # Compute m x N matrix of model errors at time t using AR(1) scheme

        X_f += Mod_err                                                  # m x N matrix of state forecasts
        s = np.mean(X_f, axis = 1)                                      # Compute m x 1 vector of mean state forecast
        S[t + 1, :] = s                                                 # Store mean state forecast
        Meas_err = np.random.multivariate_normal(np.zeros(m), R, N).T   # Draw m x N matrix of measurement errors at time t
        D = np.add(plugin['Y'][t+1], Meas_err)                          # Generate m x N matrix of measured data replicates at time t + 1
        A = X_f - s[:, np.newaxis]                                      # Step 3: Subtract mean state from each vector (ensemble member) of Xf
        C = 1 / (N - 1) * np.dot(A, A.T)                                # Step 4: Compute m x m sample covariance matrix C of state forecasts
        K = np.linalg.solve(C + R, C)                                   # Step 5: Compute the m x m Kalman gain matrix
        X_a = X_f + np.dot(K, (D - X_f))                                # Step 6: Compute m x N matrix of analysis states
    
    s = S.flatten()                                                     # Return mean state forecast at each time as a single vector

    return s
