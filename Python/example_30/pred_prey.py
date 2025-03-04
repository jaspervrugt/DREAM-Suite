import numpy as np
from scipy.integrate import odeint


############################
### Case study 30
############################
def pred_prey(x, plugin):
    """
    This function runs the 1-predator-1-prey model.
    """
    # Solve the system of differential equations using odeint
    u = odeint(PP_ode_equations, plugin['u0'], plugin['t'], args=(x,))
    # Two columns of predator and prey abundances
    # Need to re-order the nx2 matrix to yield two columns underneath eachother
    u = u.T 
    # Return the simulated abundances of both species as one vector
    return u.flatten()


def PP_ode_equations(u, t, pars):
    """
    Solves the coupled ordinary differential equations of the 1-predator-1-prey population model.
    """
    # Initial states
    x = u[0]            # Prey population
    y = u[1]            # Predator population
    
    # Parameters
    r = pars[0]         # Intrinsic rate of prey natural increase
    alfa = pars[1]      # Proportionality constant linking prey mortality to prey and predator numbers
    m = pars[2]         # Mortality rate of predators
    theta = pars[3]     # Proportionality constant linking predator increase to prey and predator numbers
    K = 50              # Maximum number of preys (carrying capacity)
    
    # Ordinary differential equations
    dxdt = r * x * (1 - x / K) - alfa * x * y
    dydt = -m * y + theta * x * y
    
    # Return the derivatives as a vector
    return [dxdt, dydt]
