import numpy as np
from scipy.integrate import solve_ivp


############################
### Case study 26
############################
def EnKF_Lorenz63(par, plugin):
    # Lorenz, E.N. (1963), Deterministic non-periodic flow, Journal of
    #   Atmospheric Sciences, 20, 130–141, 1963
    # Vrugt, J.A., C.G.H. Diks, H.V. Gupta, W. Bouten, and J.M. Verstraten
    #   (2005), Improved treatment of uncertainty in hydrologic modeling: 
    #   Combining the strengths of global optimization and data assimilation, 
    #   Water Resources Research, 41, W01017, doi:10.1029/2004WR003059
    if not hasattr(EnKF_Lorenz63, "initialized"):   # Store local variables in memory
        EnKF_Lorenz63.N = plugin['N']               # ensemble size EnKF
        EnKF_Lorenz63.m = plugin['m']               # number of measurements at given time
        EnKF_Lorenz63.R = plugin['R']               # m x m measurement error covariance matrix
        EnKF_Lorenz63.C = plugin['C']               # m x m model error covariance matrix at given time
        def Lorenz63(t, u, sigma, rho, beta):       # Lorenz 1963 system of equations
            return np.array([ -sigma*u[0] + sigma*u[1], 
                            rho*u[0] - u[1] - u[0]*u[2], 
                            -beta*u[2] + u[0]*u[1] ])
        EnKF_Lorenz63.Lorenz63 = Lorenz63
        EnKF_Lorenz63.initialized = True            # Flag to indicate that initialization is comple

    sigma, rho, beta, AR_rho = par                  # Extract parameters from 'par'
    
    N = EnKF_Lorenz63.N
    m = EnKF_Lorenz63.m
    R = EnKF_Lorenz63.R
    C = EnKF_Lorenz63.C

    E = np.full((plugin['NT'], m), np.nan)                              # Initialize return matrix of forecast states                
    E[0, :] = plugin['Y'][:m, 0]                                        # Set the first mean simulated value equal to observed value
    Meas_err = np.random.multivariate_normal(np.zeros(m), R, N).T       # Draw m x N matrix of measurement errors from N_m(0, R)
    X_a = np.add(plugin['Y'][:m, 0][:, np.newaxis], Meas_err)           # m x N matrix of initial states of ensemble members: min of zero
    X_f = np.nan * np.zeros((m, N))                                     # m x N matrix of state forecasts of ensemble members
    Mod_err = np.random.multivariate_normal(np.zeros(m), C, N).T        # Draw m x N matrix with model errors
    for t in range(plugin['NT'] - 1):                                   # Now run the model one-observation ahead
        for ii in range(N):                                             # Generate state forecast for each ensemble member
            # Solve the ODE using solve_ivp
            result = solve_ivp(EnKF_Lorenz63.Lorenz63, [t * plugin['dt'], \
                (t + 1) * plugin['dt']], X_a[0:m, ii], \
                args = (sigma, rho, beta), method = 'RK45', \
                t_eval = [t * plugin['dt'], (t + 1) * plugin['dt']])    # Run the model from one time step to next for the nth state vector
            # state_forecast = ode45(Lorenz_1963, [(t-1) * plugin['dt'], t * plugin['dt']], X_a[:m, ii], plugin['options'])
            state_forecast = result.y[:, -1]                            # Get the state forecast at the final time
            X_f[0:m, ii] = state_forecast                               # Store the simulated state forecast at final time
#            X_f[:m, ii] = state_forecast[-1, :m]
            Mod_err[:m, ii] = AR_rho * Mod_err[:m, ii] + \
                np.sqrt(1 - AR_rho**2) * np.random.multivariate_normal(np.zeros(m), C)   # Compute m x N matrix of model errors at time t using AR(1) scheme

        X_f += Mod_err                                                  # m x N matrix of state forecasts
        e = np.mean(X_f, axis = 1)                                      # Compute m x 1 vector of mean state forecast
        E[t + 1, :] = e                                                 # Store mean state forecast
        Meas_err = np.random.multivariate_normal(np.zeros(m), R, N).T   # Draw m x N matrix of measurement errors at time t
        D = np.add(plugin['Y'][:m, t+1][:, np.newaxis], Meas_err)       # Generate m x N matrix of measured data replicates at time t + 1
        A = X_f - e[:, np.newaxis]                                      # Step 3: Subtract mean state from each vector (ensemble member) of Xf
        C = 1 / (N - 1) * np.dot(A, A.T)                                # Step 4: Compute m x m sample covariance matrix C of state forecasts
        K = np.linalg.solve(C + R, C)                                   # Step 5: Compute the m x m Kalman gain matrix
        X_a = X_f + np.dot(K, (D - X_f))                                # Step 6: Compute m x N matrix of analysis states
    
    e = E.flatten()                                                     # Return mean state forecast at each time as a single vector
    
    return e

