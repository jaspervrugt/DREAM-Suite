import numpy as np
from scipy.optimize import minimize
from scipy.stats import percentileofscore, cumfreq
from hmodel import *

############################
### Case study 14
############################
def hmodel_summary_metrics(par, plugin):
    """
    hmodel simulation of discharge according to Schoups et al. 2010
    Returns summary statistics of simulated discharge record

    Parameters:
        par (array): Parameter values for the model.
        plugin (dict): Dictionary containing plugin information, including time vector, model options, observed data, etc.
    
    Returns:
        F (array): RMSE of driven and non-driven part hydrograph.
        Y_sim (array): Simulated discharge after burn-in period.
    """
    
    # Run model function (assumed to be implemented in Python or using C extension)
    y_sim = hmodel(par, plugin)

    # Calculate summary statistics
    S_sim = calc_signatures(y_sim).T 
    
    return S_sim, y_sim


def calc_signatures(y, method = 1):
    """
    Calculate the summary metrics of measured/simulated discharge record.

    Parameters:
        y (numpy array): nx1 vector of discharge values
        method (int): Computational method for Flow Duration Curve (= FDC)
                      1 -> Weibull plotting position
                      other -> Built-in plotting position

    Returns:
        S (numpy array): 1x4 vector with summary metrics of streamflow record
                         [S(1), S(2), S(3), S(4)]
                         S(1) -> Annual runoff index
                         S(2) -> Annual baseflow index
                         S(3) -> Air-entry (= parameter 1) of 2-parameter VG-inspired FDC function
                         S(4) -> Slope (= parameter 2) of 2-parameter VG-inspired FDC function
    """
    # Load local data
    try:
        daily_data = np.loadtxt('03451500.dly')  # Adjust path to your data
        P = daily_data[730:, 3]  # Measured precip (2y spin-up)
    except FileNotFoundError:
        raise FileNotFoundError("File '03451500.dly' not found. Ensure the path is correct.")
    
    ny = len(y)  # Number of streamflow data points
    fi = 0.925  # Low pass filter (baseflow index)
    
    # Make sure y is a column vector (1D array)
    y = np.array(y).reshape(-1)

    # Initialize output array for summary metrics
    S = np.nan * np.ones(4)
    S[0] = np.sum(y) / np.sum(P)  # Annual runoff ratio (S1)

    # Initialize baseflow contribution using low pass filter
    yb = np.full(ny, np.nan)
    yb[0] = y[0] / 4

    # Apply low-pass filter for baseflow
    for j in range(1, ny):
        yb[j] = min(fi * yb[j - 1] + 0.5 * (1 - fi) * (y[j] + y[j - 1]), y[j])

    S[1] = np.sum(yb) / np.sum(y)  # Annual baseflow index (S2)

    # Calculate empirical FDC and exceedance probabilities
    Y_s, E_p = calc_efdc(y, ny, method)

    # Initial parameter values for Nelder-Mead optimization
    x0 = [1, 1]
    
    # Optimize parameters (air-entry and slope) for VG-inspired FDC fitting
    result = minimize(lambda x: VG_Ret(x, Y_s, E_p, ny), x0, method='Nelder-Mead')
    S[2:4] = result.x  # Extract fitted parameters (S3, S4)

    return S


def calc_efdc(Y, nY, method):
    """
    Calculate exceedance probability of each flow level.
    
    Parameters:
        Y (numpy array): The discharge data
        nY (int): The number of data points
        method (int): The method to use for FDC calculation
                      1 -> Weibull plotting position
                      other -> Built-in empirical CDF

    Returns:
        Y_s (numpy array): Sorted discharge values
        E_p (numpy array): Exceedance probabilities
    """
    if method == 1:
        Y_s = np.sort(Y)  # Sort discharge in ascending order
        E_p = 1 - (np.arange(1, nY + 1) - 0.5) / nY  # Weibull exceedance probabilities
    else:
        # Built-in empirical CDF
        cdf = cumfreq(Y, numbins=nY)
        E_p = 1 - cdf.cumcount / nY  # Exceedance probabilities

    return Y_s, E_p


def VG_Ret(pars, Y_s, E_ptrue, nY):
    """
    Compute the RMSE between actual exceedance probabilities (E_ptrue)
    and the VG model's predicted exceedance probabilities (E_ppred).

    Parameters:
        pars (list or array): The parameters (alpha, n)
        Y_s (numpy array): Sorted discharge values
        E_ptrue (numpy array): True exceedance probabilities
        nY (int): The number of data points

    Returns:
        RMSE (float): The root mean square error between true and predicted E_p values
    """
    alpha, n = pars  					# Extract parameter values
    E_p = (1 + (alpha * Y_s)**n)**(1/n - 1)  		# Exceedance probability according to VG WRF
    RMSE = np.sqrt(np.sum((E_p - E_ptrue)**2) / nY)  	# Calculate RMSE

    return RMSE
