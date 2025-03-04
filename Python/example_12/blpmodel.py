import numpy as np
from scipy.stats import multivariate_normal
from scipy.sparse import lil_matrix, csr_matrix
from numpy.linalg import cholesky, eig
import scipy.special as sp   
import matplotlib.pyplot as plt

############################
### Case study 12
############################
def blpmodel(par):

    if not hasattr(blpmodel, "initialized"):                                        # Store local variables in memory
        data = np.loadtxt('forest.txt')                                             # Forest data
        x = data[:, 0]                                                              # x coordinates
        y = data[:, 1]                                                              # y coordinates
        blpmodel.z = data[:, 2]                                                     # z values (the observed data)
        X = np.column_stack((x, y))                                                 # Create matrix of coordinates
        blpmodel.M = np.ones((len(blpmodel.z), 1))                                  # Create trend vector (constant term)
        blpmodel.H = distmat(X)                                                     # Compute the distance matrix between all observations
        blpmodel.np = 1                                                             # Number of parameters in the linear model
        blpmodel.nt = 5                                                             # Number of parameters in the variogram function
        blpmodel.initialized = True                                                 # Flag to indicate that initialization is complete
    
    beta = par[:blpmodel.np]                                                        # Parameters of the linear model
    S = blpmodel.M @ beta                                                           # The linear model
    th = par[blpmodel.np:blpmodel.nt]                                               # Parameters of the variogram
    C = cova_matern(th[0], th[1], th[2], th[3], blpmodel.H)                         # Compute covariance structure of the data using Matern covariance
    dC = np.linalg.det(C) + np.finfo(np.float64).tiny                               # Use realmin for numerical stability
    L1 = np.dot((blpmodel.z - S).T, np.dot(np.linalg.inv(C), (blpmodel.z - S)))     # Calculate the L1 norm
    loglik = -0.5 * (np.log(dC) + L1)                                               # Calculate the log-likelihood
    
    return loglik


def calc_vario(xd, yd, z, nlag = 30, hmax = 0, tolag = 0):
    """
    Calculate distance matrix & empirical variogram.
    
    Parameters:
    xd     : numpy array - x coordinates
    yd     : numpy array - y coordinates
    z      : numpy array - data values
    nlag   : int (optional) - number of lags (default is 30)
    hmax   : float (optional) - maximum distance for variogram calculation (default is 0)
    tolag  : float (optional) - percent of lag tolerance (default is 0)
    
    Returns:
    H      : numpy array - distance matrix for all points
    G      : numpy array - semivariance matrix for all points
    semv   : numpy array - empirical semivariogram [distance, semivariance, no. points, std.dev of semivar]
    vcloud : numpy array - variogram cloud [distance, semivariance]
    """
    
    X = np.column_stack((xd, yd))
    n = len(X)
    # Calculate the distance matrix (H)
    x = np.reshape(X, (n, 1, 2))
    y = np.reshape(X, (1, n, 2))
    # Use np.tile to repeat arrays along the necessary dimensions
    x_tiled = np.tile(x, (1, n, 1))  # Shape becomes (108, 108, 2)
    y_tiled = np.tile(y, (n, 1, 1))  # Shape becomes (108, 108, 2)
    # Calculate the squared differences, sum along the third dimension (axis=2), then take the square root
    H = np.sqrt(np.sum((x_tiled - y_tiled)**2, axis=2))
    # Calculate the semivariance matrix (G)
    z1 = np.reshape(z, (n, 1))  # Shape: (n, 1) 
    z2 = np.reshape(z, (1, n))  # Shape: (1, n)
    # Subtract z1 from z2 using broadcasting
    # result = z1 - z2  # Resulting shape: (n, n)
    G = 0.5 * (((z1 - z2)**2))

    # Extract distance & gamma values to create the variogram cloud
    row, col = np.where((np.tril(G) > 0).T == True) 
    ij_linear = np.ravel_multi_index((row, col), G.shape)
    hd = H.ravel()[ij_linear]
    gd = G.ravel()[ij_linear]
    vcloud = np.column_stack((hd, gd))
    # Calculate empirical variogram using method of moments
    if hmax == 0:
        hmax = 0.8 * np.max(hd)
    
    step = hmax / (nlag + 1)
    xtol = step * tolag / 100               # Lag tolerance
    h = np.arange(0, hmax + 1, step)
    nlag = len(h)
    semv = np.zeros((nlag, 4))
    
    for i in range(nlag - 1):
        ll = h[i] - xtol
        ul = h[i] + step + xtol
        ij = np.where((hd >= ll) & (hd <= ul))[0]
        nij = len(ij)
        
        if nij > 0:
            semv[i, 0] = np.mean(hd[ij])        # Mean distance
            semv[i, 1] = np.mean(gd[ij])        # Mean semivariance
            semv[i, 2] = np.std(gd[ij],ddof=1)  # Standard deviation of the semivariance
            semv[i, 3] = nij                    # Number of pairs
    
    semv = semv[semv[:, 3] > 1]             # Filter semv to exclude rows with less than 2 pairs

    return H, G, semv, vcloud


def vmatern(c0, c1, range_, nu, H, calc_method = 1):
    """
    Matern variogram function.
    
    Parameters:
    c0    : float - sill (the range value at which the variogram flattens)
    c1    : float - nugget (small value added to prevent divergence)
    range_ : float - the range parameter (distance at which the variogram reaches the sill)
    nu    : float - smoothness parameter
    H     : numpy array - separation distances
    
    Returns:
    V     : numpy array - variogram values
    """
    Hr = H / range_                     # Normalize the distance by range
    r1 = 2**(nu - 1) * sp.gamma(nu)     # Scale factor (Gamma function)

    if calc_method == 1:    ## MATLAB implementation
        bes = sp.kv(nu, Hr)                 # creates warning in kv as H = 0 on diagonal elements
        C = (1 / r1) * (Hr ** nu) * bes
    elif calc_method == 2:  ## Python implementation [avoid warning]
        Hr_copy = Hr.copy()
        np.fill_diagonal(Hr_copy, 0.001)    # avoids warning by setting diagonal element (H == 0) to small value 
        bes_copy = sp.kv(nu, Hr_copy)                  
        C = (1 / r1) * (Hr_copy ** nu) * bes_copy

    V = c0 + c1 * (1 - C)               # Variogram
    V[H == 0] = c0                      # Variogram value should be c0 when the distance is zero
    
    return V


def cova_matern(c0, c1, range_, nu, H, calc_method = 2):
    """
    Covariance function of the Matern model.
    
    Parameters:
    c0    : float - nugget effect (value at zero distance)
    c1    : float - sill (value at large distances)
    range : float - range parameter
    nu    : float - smoothness parameter
    H     : numpy array - distances
    
    Returns:
    S     : numpy array - covariance values corresponding to the distances in H
    """
    Hr = H / range_
    
    if nu == 0.5:                               # Use exponential model for nu = 0.5
        S = c1 * np.exp(-Hr)
    else:
        r1 = (2 ** (nu - 1)) * sp.gamma(nu)     # scalar
        if calc_method == 1:    ## MATLAB implementation
            bes = sp.kv(nu, Hr)                 # creates warning in kv as H = 0 on diagonal elements
            F = (1 / r1) * (Hr ** nu) * bes
        elif calc_method == 2:  ## Python implementation [avoid warning]
            Hr_copy = Hr.copy()
            np.fill_diagonal(Hr_copy, 0.001)    # avoids warning by setting diagonal element (H == 0) to small value 
            bes_copy = sp.kv(nu, Hr_copy)                  
            F = (1 / r1) * (Hr_copy ** nu) * bes_copy

        S = c1 * F
            
    S[H == 0] = c0 + c1                         # Apply nugget effect for zero distance [resolves issue with Inf of bes]
    
    return S


def distmat(x, y = None):
    """
    Compute the distance matrix (Euclidean distance) between two sets of points.
    
    Parameters:
        x (ndarray): A 2D array of size (m, p), where m is the number of points and p is the dimension of each point.
        y (ndarray, optional): A 2D array of size (n, p), where n is the number of points in the second set. If not provided, y is assumed to be the same as x.
        
    Returns:
        ndarray: A 2D distance matrix of size (m, n) or (m, m) if y is not provided.
    """
    
    if y is None:  # If only one matrix is provided, compute the distance matrix of X with itself
        
        if x.ndim != 2:
            raise ValueError("Input must be a 2D matrix.")
        
        m = x.shape[0]
        d = np.zeros((m, m))  # initialize the distance matrix
        # Use broadcasting to calculate the pairwise Euclidean distances
        for i in range(m):
            for j in range(i+1, m):
                d[i, j] = np.linalg.norm(x[i, :] - x[j, :]) # Euclidean distance
                d[j, i] = d[i, j]                           # The distance matrix is symmetric
                
    else:  # If both X and Y matrices are provided
        if x.ndim != 2 or y.ndim != 2:
            raise ValueError("Both inputs must be 2D matrices.")
        
        mx, nx = x.shape
        my, ny = y.shape
        if nx != ny:
            raise ValueError("Both matrices must have the same number of columns.")
        
        # Compute the distance matrix between X and Y
        d = np.zeros((mx, my))  # Distance matrix with shape (m, n)
        for i in range(mx):
            for j in range(my):
                d[i, j] = np.linalg.norm(x[i, :] - y[j, :])  # Euclidean distance
    
    return d


def postproc_variogram(DREAMPar, P):

    idx = np.argmax(P[:, -1])                               # Find index of max log-density
    idx = idx[0] if isinstance(idx, np.ndarray) else idx    # Handle cases where idx is an array
    Pars = P[int(0.75 * len(P)):, :DREAMPar['d']]           # Take the last 25% of the posterior samples

    # Displaying message
    print('POSTPROC_VARIOGRAM: PLEASE RESORT TO TABLES AND FIGURES OF DREAM POSTPROCESSOR FOR RESULTS INFERENCE')
    print('POSTPROC_VARIOGRAM: THIS POSTPROCESSOR PLOTS THE FINAL VARIOGRAM AND ITS 95% PREDICTION INTERVALS')

    # Load the forest data (replace with actual file path)
    data = np.loadtxt('forest.txt')
    x, y, z = data[:, 0], data[:, 1], data[:, 2]        # Define coordinates
    _, _, semv, _ = calc_vario(x, y, z)                 # Calculate the empirical variogram

    N = Pars.shape[0]                                   # Number of posterior samples to use
    hd = np.arange(0, 801, 2)                           # Separation distance (hd)
    # Identify unique posterior solutions
    B, I, J = np.unique(Pars[:, :DREAMPar['d']], axis = 0, return_inverse = True, return_index = True)
    V = np.zeros((B.shape[0], len(hd)))                 # Loop over each sample and calculate variograms
    for j in range(B.shape[0]):
        # Assuming vmatern function returns the variogram for given parameters
        V[j, :] = vmatern(B[j, 1], B[j, 2], B[j, 3], B[j, 4], hd)

    V = V[J, :]                                         # Rearrange matrix sim_out
    V_low = np.percentile(V, 2.5, axis=0)               # 2.5 percentile posterior variogram
    V_high = np.percentile(V, 97.5, axis=0)             # 97.5 percentile posterior variogram
    fig, ax = plt.subplots(figsize=(12, 8))             # Plot results
    # Plot the mean variogram and percentiles
    ax.plot(hd, np.mean(V, axis = 0), '-', color = (0.5, 0.5, 0.5), linewidth = 3)
    ax.plot(hd, V_low, 'k-.', linewidth = 3)
    ax.plot(hd, V_high, 'k-.', linewidth = 3)
    # Plot empirical variogram
    ax.plot(semv[:, 0], semv[:, 1], 'r.', markersize = 20)

    # Labeling the axes
    ax.set_xlabel('Separation distance (lag) [m]', fontsize = 18, fontweight = 'bold')
    ax.set_ylabel('Semivariance', fontsize = 18, fontweight = 'bold')
    ax.set_title('Posterior variogram: prediction uncertainty and empirical variogram', fontsize = 18, fontweight = 'bold')
    # Add interval
    ax.set_xlim([0, 800])
    ax.set_ylim([0, 800])
    ax.minorticks_on()          
    # Adding legends
    ax.legend(['Mean Variogram', '2.5% Percentile', '97.5% Percentile', 'Empirical Variogram'], loc = 'upper right', fontsize = 14)
    ax.grid(True)
    plt.show()

    return 