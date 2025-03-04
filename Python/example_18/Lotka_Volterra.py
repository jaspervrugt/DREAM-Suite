import numpy as np
from scipy.integrate import odeint


############################
### Case study 18
############################
def Lotka_Volterra(x):
    # ####################################################################### #
    # LOTKA_VOLTERRA Analytic solution predator-prey populations in food web  #
    # Please check the work of Massoud et al. in Ecology Letters, 2017        #
    # ####################################################################### #

    if not hasattr(Lotka_Volterra, "initialized"):                                  # Store local variables in memory
        Lotka_Volterra.y0 = [30, 4]                                                 # Initial conditions: initial population of prey and predator
        Lotka_Volterra.tout = np.arange(0, 20 + 1/12, 1/12)                         # Time vector: from 0 to 20 years with monthly intervals
        def dydt_func(y, t, alpha, beta, gamma_, delta):                            # Define the Lotka-Volterra equations as a function
            dydt1 = alpha * y[0] - beta * y[0] * y[1]                               # y[0] = prey population, y[1] = predator population
            dydt2 = -gamma_ * y[1] + delta * y[0] * y[1]
            return [dydt1, dydt2]                                                   # Return dydt of the two states

        Lotka_Volterra.dydt = dydt_func                                             # Store the dydt function
        Lotka_Volterra.initialized = True                                           # Flag to indicate that initialization is complete        

    y = Lotka_Volterra.y0                                                           # Unpack the parameters (alpha, beta, gamma_, delta) and add initial states
    result = odeint(Lotka_Volterra.dydt, y, Lotka_Volterra.tout, args = tuple(x))   # Solve the system of differential equations using odeint
    Y = result[1:, :2]                                                              # Return only the prey and predator populations
    Y_vector = Y.T.reshape(-1,1)                                                    # Flatten to return as a one-dimensional vector

    return Y_vector.flatten(), Y_vector.flatten()
