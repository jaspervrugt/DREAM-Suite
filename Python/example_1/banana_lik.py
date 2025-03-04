import numpy as np
from scipy.stats import multivariate_normal

############################
### Case study 1
############################
def banana_lik(x, plugin = None):
    # ####################################################################### #
    # BANANA_LIK d-variate banana-shaped multivariate log-likelihood.         #
    # This is equivalent to a twisted normal distribution                     # 
    #   Haario, H., E. Saksman, and J. Tamminen (1999), Adaptive proposal     #
    #       distribution for random walk Metropolis algorithm. Computational  # 
    #       Statistics 14, 375–395, https://doi.org/10.1007/s001800050022     #
    # The d-vector x enters as a horizontal vector: x = (x_1,...,x_d)         #
    # Function also works if X is a Nxd matrix of N proposals of x            #
    #                                                                         #
    # Let f_N(µ,Σ) be the density of the multivariate normal distribution     #
    # N(µ,Σ) with µ-mean and Σ covariance matrix. Then, the density of the    #
    # twisted normal distribution of Haario et al. (1999) is given by         #
    # f_{b} = f_N(µ,Σ) ◦ γ_{b} where µ = (0, 0, ..., 0_{d}) is a zero-mean    #
    # vector, the covariance matrix Σ = diag(100,1,...,1_{d}) is a dxd        #
    # diagonal matrix and the  function γ_{b} = (x_{1},x_{2} + b*x_{1}^{2}    #
    # - 100b, x_{3}, ..., x_{d})  with b = 0.1.                               #
    # ####################################################################### #
    
    global b, d, invC, log_Z

    x = np.array([x]).reshape(1,-1)                         # added for Python

    # Check if the persistent variables are initialized
    if 'invC' not in globals():                             # Local memory
        d = x.shape[1]                                      # Dimensionality of the target distribution
        b = 0.1                                             # Nonlinearity factor for the banana distribution
        C = np.eye(d)                                       # Covariance matrix (identity matrix by default)
        C[0, 0] = 100                                       # Adjust first element for the banana-shaped distribution
        invC = np.linalg.inv(C)                             # Inverse of the covariance matrix
        log_Z = np.log((2 * np.pi) ** (-d / 2) \
                       * np.linalg.det(C) ** (-1 / 2))      # Normalizing constant
        if d > 150:                                         # If d > 150, set log_Z to zero to avoid numerical issues
            log_Z = 0                       

    x[:, 1] = x[:, 1] + b * x[:, 0] ** 2 - 100 * b          # Apply banana-shaped nonlinearity: x[1] adjusted by x[0]^2

    # Calculate log-likelihood
    if d == 1:  ## Univariate normal
        loglik = log_Z - 0.5 * np.dot(x, np.dot(invC, x.T))
    else:       ## Multivariate normal
        loglik = log_Z - 0.5 * np.sum(x @ invC * x, axis = 1)
    
    return loglik

