import numpy as np

############################
### Case study 8
############################
def ABC_func(x):
    """
    Simple one-dimensional example to illustrate ABC.

    Parameters:
        x (list or np.array): The input parameter vector (1D array).

    Returns:
        rho (float): The computed distance value.
    """
    # Draw 100 numbers from standard normal distribution
    Y = np.random.normal(x[0], 1, 100)
    # Define distance function
    if np.random.rand() < 0.5:
        # Take the mean of the absolute values of Y
        rho = np.abs(np.mean(Y))
    else:
        rho = np.abs(Y[0])
    
    return rho
