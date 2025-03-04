import numpy as np
from scipy.stats import multivariate_normal
#from scipy.linalg import inv, det

############################
### Case study 2
############################
def normal_lik(x, plugin = None, comp_method = 3):
    # ####################################################################### #
    # NORMAL_LIK d-variate normal log-likelihood, N(µ,Σ) with µ = 0 and       #
    # Σ = dxd covariance matrix with 1:d on main diagonal and 0.5 correlation #
    # between dimensions                                                      #
    # comp_method:  1: Exact using built-in functions [slow > 10 years ago]   #
    #               2: Some protection against underflow                      #
    #               3: Some protection against underflow [faster, Note 1]     #
    # What is underflow? L(x) → zero when x is far removed from target mean   #
    #                         → log(L(x)) goes to -∞                          #
    # Note 1: comp_method = 2: N = 1e5; d = 100; X = 20*rand(N,d) - 10;       #
    # MATLAB crashes, matrices become too large                               #
    # --> use comp_method 3 (= faster for large N and d)                      #
    # ####################################################################### #
    
    # Persistent variables (treated as global here)
    global d, mu, C, log_Z, invC

    x = np.array([x]). reshape(1,-1)                                        		# added for Python
    if 'invC' not in globals():                                             		# Determine constant variables, only if not already set
        d = x.shape[1]                                                      		# Number of parameters (dimensions)
        mu = np.zeros(d)                                                    		# Mean vector (0)
        R = 0.5 * np.eye(d) + 0.5 * np.ones((d, d))                         		# Correlation matrix R
        C = np.nan * np.ones((d, d))                                        		# Covariance matrix C
        for i in range(d):
            for j in range(d):
                C[i, j] = R[i, j] * np.sqrt((i+1) * (j+1))
        invC = np.linalg.inv(C)                                             		# Inverse of covariance matrix
        log_Z = np.log(((2 * np.pi) ** (-d / 2)) * np.linalg.det(C) ** (-1 / 2))      	# Normalization constant (log(Z))

    x = x - mu                                                              		# Scale candidate point(s)

    # Calculate log-likelihood based on the chosen computation method
    if comp_method == 1:    ## Exact computation using built-in functions (slow)
        loglik = np.log(multivariate_normal.pdf(x, mean=mu, cov=C))
    elif comp_method == 2:  ## Protection against underflow (moderately fast)
        loglik = log_Z - 0.5 * np.diagonal(np.dot(np.dot(x, invC), x.T))
    elif comp_method == 3:  ## Fastest method with protection against underflow
        if d == 1:
            loglik = log_Z - 0.5 * invC * (x ** 2)
        else:
            loglik = log_Z - 0.5 * np.sum(x.T * np.dot(invC, x.T), axis = 0)
        
    return loglik