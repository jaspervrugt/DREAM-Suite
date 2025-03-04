import numpy as np
from scipy.stats import multivariate_normal

############################
### Case study 3, 10, 29
############################
def mixture_lik(x, plugin = None, comp_method = 1):
    # ####################################################################### #
    # MIXTURE_LIK d-variate normal mixture log-likelihood.                    #
    # Target distribution is w1*N(µ1,Σ1) + w2*N(µ2,Σ2) with µ1 = -5, µ2 = 5,  #
    #    Σ1/Σ2 = identity matrix and w1 = 1/3 and w2 = 2/3                    #
    # comp_method   1: Exact using built-in functions (slow)                  #
    #               2: Fast but without normalization constant [= OK]         #
    #               3: Refined/fast [= some protection against underflow]     #
    # What is underflow? L(x) → zero when x is far removed from target mean   #
    #                         → log(L(x)) goes to -∞.                         #
    # For example, if                                                         #
    # x = [11.8 17.0 19.8 18.1 24.7 23.5 11.0 4.8 29.4 0.02]                  #
    # then comp_method = 1 yield loglik = -Inf, whereas comp_method =3 yields #
    # loglik = -993.3101.                                                     #
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
    global d, k, mu, log_Z, w1, w2, C, Cc
    
    # added for Python
    # x = np.array([x]).reshape(1,-1)
    
    # Determine constant variables, only if not already set
    if 'k' not in globals():
        
        d = x.shape[0]                                                          # Dimension of the data
        k = 2                                                                   # Number of Gaussians in the mixture
        w1 = 1/3                                                                # Weights of the Gaussians
        w2 = 1 - w1                                                             # Means of the two Gaussian components
        log_prior = np.log([w1, w2])                                            # Log prior density
        mu = np.array([[-5]*d, [5]*d])                                          # Mean of components
        C = np.eye(d)                                                           # Covariance matrix (identity)
        Cc = np.linalg.cholesky(C)                                              # Cholesky factorization of C
        logDetSigma = 2 * np.sum(np.log(np.diag(Cc)))                           # Logarithm of determinant of C
        log_Z = -0.5 * logDetSigma + log_prior - d * np.log(2 * np.pi) / 2      # Normalization constants for each Gaussian
        
    ll_max = lik = 0                                                            # Initialize the log-likelihood
    # Switch between different computation methods
    if comp_method == 1:    ## Method 1: Exact using built-in functions (slow)
        lik = w1 * multivariate_normal.pdf(x, mean=mu[0], cov=C) + w2 * multivariate_normal.pdf(x, mean=mu[1], cov=C)
    elif comp_method == 2:  ## Method 2: Faster but without normalization constant (no underflow protection)
        lik = w1 * np.exp(-np.sum((x + 5)**2, axis = 1)) + w2 * np.exp(-np.sum((x - 5)**2, axis = 1))
    elif comp_method == 3:  ## Method 3: Refined/fast with protection against underflow
        ll_j = np.full(k, np.nan)
        for j in range(k):  ## Calculate log-likelihood for the j-th component
            ll_j[j] = -0.5 * np.sum(((x - mu[j]) @ np.linalg.inv(Cc)) ** 2).squeeze() + log_Z[j]

        ll_max = np.max(ll_j)                                                   # Maximum of log-likelihood
        lik = np.sum(np.exp(ll_j - ll_max), axis = 0)                           # Avoid numerical underflow    

    loglik = np.log(lik) + ll_max                                               # Final log-likelihood
    
    return loglik