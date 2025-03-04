import numpy as np
from scipy.stats import norm, gamma, uniform, genextreme, genpareto

############################
### Case study 19
############################
def BMA_lik(par, plugin):
    """
    This function computes the BMA log likelihood for weights and sigma's.

    Parameters:
    x (numpy array): A vector containing the weights and sigma parameters.
    plugin (dict): A dictionary containing various input values for the model.
    
    Returns:
    loglik (float): The log-likelihood of the BMA model.
    """
      
    D = plugin['BMA']['D']  # Ensemble forecasts
    y = plugin['BMA']['y']  # Verifying data
    n, K = D.shape          # Number of forecasts (n) and number of ensemble members (K)
    
    w = par[:K]             # Unpack weights
    w = w / np.sum(w)       # Normalize weights along rows (axis 1)
                            # This normalization is only done to illustrate eDREAM Package
                            # To get "true" posterior, weight normalization needs to be done in eDREAM Package
                            # This is done in the MODELAVG toolbox! Please use this for BMA method
    # Variance options
    VAR = plugin['BMA']['VAR']
    
    if VAR == '1':      # Common constant variance
        sigma = par[K] * np.ones((n, K))
    elif VAR == '2':    # Individual constant variance
        sigma = np.multiply(par[K:2*K], np.ones((n, K)))
    elif VAR == '3':    # Common non-constant variance
        c = par[K]
        sigma = c * D
    elif VAR == '4':    # Individual non-constant variance
        c = par[K:2*K]
        sigma = np.multiply(c, D)
    
    sigma = np.maximum(sigma, np.finfo(float).eps)  # Ensure sigma is not zero
    
    # Conditional distribution options
    PDF = plugin['BMA']['PDF']
    
    # Calculate likelihoods
    Y = np.tile(y, (K, 1)).T  # Make K copies of verifying data
    
    if PDF == 'normal':  # Gaussian distribution
        A = D
        B = sigma
        L = norm.pdf(Y, A, B)

    elif PDF == 'gamma':  # Gamma distribution
        mu = np.abs(D)
        var = sigma ** 2
        A = mu ** 2 / var
        B = var / mu
        L = gamma.pdf(Y, A, scale=B)

    # Calculate likelihoods
    lik = np.sum(L * w, axis=1) + np.finfo(float).tiny      # BMA likelihoods
    loglik = np.sum(np.log(lik))
    
    return loglik


def setup_BMA(DREAMPar, Par_info, D, y, VAR):
    """
    Setup BMA mixture model estimation
    
    Parameters:
    - DREAMPar: Dictionary containing parameters
    - Par_info: Dictionary to store parameter information
    - D: NumPy array of forecasts (shape: n x K)
    - y: NumPy array of verifying data (shape: n x 1)
    - VAR: String indicating variance type ('1', '2', '3', or '4')

    Returns:
    - DREAMPar: Updated dictionary with new parameters
    - Par_info: Updated dictionary with parameter information
    - D_bc: Bias-corrected forecasts
    - A: Intercepts for bias correction
    - B: Slopes for bias correction
    """
    
    # Get number of forecasts (n) and number of ensemble members (K)
    n, K = D.shape
    
    adjust = True  # Linear bias correction
    
    # Perform bias correction using the ComputeAB function
    D_bc, A, B = ComputeAB(D, y, n, K, adjust)
    
    # Initialize parameter names
    par_name = [f'\\beta_{z+1}' for z in range(K)]  # Beta names for the K ensemble members
    
    # Set parameter names and adjust settings based on VAR value
    if VAR == '1':  # Common constant variance
        DREAMPar['d'] = K + 1
        Par_info['max'] = np.hstack([np.ones(K), 2 * np.std(y)])
        par_name.append('\\sigma')
    elif VAR == '2':  # Individual constant variance
        DREAMPar['d'] = 2 * K
        Par_info['max'] = np.hstack([np.ones(K), 2 * np.std(y) * np.ones(K)])
        for z in range(K):
            par_name.append(f'\\sigma_{z+1}')
    elif VAR == '3':  # Common non-constant variance
        DREAMPar['d'] = K + 1
        Par_info['max'] = np.hstack([np.ones(K), 2])
        par_name.append('c')
    elif VAR == '4':  # Individual non-constant variance
        DREAMPar['d'] = 2 * K
        Par_info['max'] = np.hstack([np.ones(K), 2 * np.ones(K)])
        for z in range(K):
            par_name.append(f'c_{z+1}')
    else:
        raise ValueError("Unknown variance option; choose between '1', '2', '3', or '4'")

    # Final adjustments to Par_info
    Par_info['names'] = par_name                    # Set parameter names
    Par_info['unit_simplex'] = 'yes'                # Weights on unit Simplex
    Par_info['min'] = np.zeros(DREAMPar['d'])       # Min. values for BMA weights/vars
    
    return DREAMPar, Par_info, D_bc, A, B


# ComputeAB function (helper for bias correction)
def ComputeAB(D, y, n, K, adjust = True):
    """
    Performs linear bias correction if indicated
    
    Parameters:
    - D: NumPy array of ensemble forecasts (shape: n x K)
    - y: NumPy array of verifying data (shape: n x 1)
    - n: Number of forecasts
    - K: Number of ensemble members
    - adjust: Boolean flag to indicate if bias correction should be applied

    Returns:
    - D_bc: Bias-corrected ensemble forecasts
    - A: Intercepts of the linear bias correction
    - B: Slopes of the linear bias correction
    """
    
    if adjust:
        B = np.zeros(K)
        A = np.zeros(K)
        for k in range(K):
            T = np.cov(D[:, k], y)[0, 1]  # Covariance
            B[k] = T / np.var(D[:, k])  # Slope of the linear regression
            A[k] = np.mean(y) - B[k] * np.mean(D[:, k])  # Intercept of the regression
    else:
        B = np.ones(K)
        A = np.zeros(K)
    
    # Compute the bias-corrected forecasts
    D_bc = np.full_like(D, np.nan)
    for k in range(K):
        D_bc[:, k] = B[k] * D[:, k] + A[k]
    
    return D_bc, A, B
