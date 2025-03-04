import numpy as np
from scipy.stats import multivariate_normal
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse import issparse
from numpy.linalg import cholesky, eig
import scipy.special as sp   

############################
### Case study 11
############################
def mvt_lik(x):
    # Multivariate Student's t-distribution log-likelihood function
    
    x = x.reshape(1,-1)
    if not hasattr(mvt_lik, "initialized"):                                                         # Store local variables in memory
        mvt_lik.d = x.shape[1]                                                                      # Number of columns in x (target dimensionality)
        mvt_lik.df = 60                                                                             # Degrees of freedom
        mvt_lik.mu = np.zeros(mvt_lik.d)                                                            # Zero-mean t-distribution with df degrees of freedom    
        A = 0.5 * np.eye(mvt_lik.d) + 0.5 * np.ones((mvt_lik.d, mvt_lik.d))                         # Construct the d x d covariance matrix
        C = np.zeros_like(A)
        for i in range(mvt_lik.d):                                                                  # Rescale to variance-covariance matrix of interest
            for j in range(mvt_lik.d):
                C[i, j] = A[i, j] * np.sqrt((i + 1) * (j + 1))                                      # MATLAB indices are 1-based
        
        s = np.sqrt(np.diag(C))                                                                     # Standardize C to correlation matrix
        if np.any(s != 1):
            C = C / (s[:, None] * s[None, :])

        mvt_lik.R = cholcov(C, 0)[0].T                                                              # Make sure C is a valid covariance matrix
        mvt_lik.logNumer_A = sp.gammaln((mvt_lik.df + mvt_lik.d) / 2) - sp.gammaln(mvt_lik.df / 2)  # Define logNumer_A
        mvt_lik.logSqrtDetC = np.sum(np.log(np.diag(mvt_lik.R)))                                    # Define logSqrtDetC
        mvt_lik.logDenom = mvt_lik.logSqrtDetC + (mvt_lik.d / 2) * np.log(mvt_lik.df * np.pi)       # Define logDenom
        mvt_lik.initialized = True                                                                  # Flag to indicate that initialization is complete
        
    Z = x @ np.linalg.inv(mvt_lik.R)                                                                # Normalize x with R (i.e., standardize the data)
    logNumer_B = - ((mvt_lik.df + mvt_lik.d) / 2) * np.log(1 + np.sum(Z**2, axis = 1) / mvt_lik.df) # Define logNumer_B
    loglik = mvt_lik.logNumer_A + np.sum(logNumer_B) - mvt_lik.logDenom                             # Calculate log-density of multivariate t-distribution

    return loglik


def cholcov(Sigma, flag=1):
    """
    Cholesky-like decomposition for covariance matrix.
    
    Args:
        Sigma (ndarray): The covariance matrix.
        flag (int): Optional flag. If flag=1 (default), computes a semi-definite decomposition.
    
    Returns:
        T (ndarray): The decomposition matrix such that Sigma = T.T @ T.
        p (int or np.nan): The number of negative eigenvalues if Sigma is not positive semi-definite, else NaN.
    """
    # Test for square, symmetric
    n, m = Sigma.shape
    wassparse = issparse(Sigma)
    tol = 10 * np.finfo(float).eps * np.max(np.abs(np.diag(Sigma)))
    
    if n == m and np.all(np.abs(Sigma - Sigma.T) < n * tol):
        try:
            T = cholesky(Sigma) 	# (, lower = False)
            p = 0
        except np.linalg.LinAlgError:
            # If not positive definite, use eigenvalue decomposition
            if flag:
                # Eigenvalue decomposition
                U, D = eig((Sigma + Sigma.T) / 2)
                
                # Pick eigenvector direction so max abs coordinate is positive
                maxind = np.argmax(np.abs(U), axis=0)
                negloc = (U[maxind + np.arange(n) * n] < 0)
                U[:, negloc] = -U[:, negloc]
                
                D = np.diag(D)
                tol = np.finfo(float).eps * np.max(D) * len(D)
                t = np.abs(D) > tol
                D = D[t]
                p = np.sum(D < 0)  # number of negative eigenvalues
                
                if p == 0:
                    T = np.diag(np.sqrt(D)) @ U[:, t].T
                else:
                    T = np.zeros((0, Sigma.shape[1]))
            else:
                T = np.zeros((0, Sigma.shape[1]))
                p = np.nan
    else:
        T = np.zeros((0, Sigma.shape[1]))
        p = np.nan

    if wassparse:
        T = sparse.csr_matrix(T)  # Return sparse if input was sparse

    return T, p
