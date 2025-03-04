import numpy as np

############################
### Case study 20
############################
def heatflow(x):
    # ####################################################################### #
    # HEATFLOW Analytic solution of the soil temperature at depth z [cm] and  #
    # time t [hr]. Adapted to hourly time scale [= not realistic in practical #
    # application]                                                            #
    # ####################################################################### #
    
    # Persistent variables z, w, t to retain in memory after the first call
    if not hasattr(heatflow, "z"):                  # Initialize variables only once
        heatflow.z = np.array([-5, -10, -15])       # Soil depth [cm]
        heatflow.w = 2 * np.pi / 24                 # Angular frequency [hr⁻¹]
        heatflow.t = np.arange(1, 49)               # Time [hr]
        # heatflow.z = np.tile(heatflow.z, (len(heatflow.t), 1))  # Repeat z for each time step
    
    T_a = x[0]  # Annual average soil temperature [°C]
    A_0 = x[1]  # Temperature amplitude surface [°C]
    fi = x[2]   # Phase constant [hr]
    d = x[3]    # Characteristic damping depth [cm]

    # Analytic solution for soil temperature
    T_s = np.full((len(heatflow.t),3), np.nan)
    for i in range(len(heatflow.z)):
        T_s[:,i] = T_a + A_0 * np.exp(heatflow.z[i] / d) * np.sin(heatflow.w * (heatflow.t - fi) + heatflow.z[i] / d)
        
    T_s_vector = T_s.T.reshape(-1, 1)    

    # Return soil temperature as a single vector
    return T_s_vector.flatten()
