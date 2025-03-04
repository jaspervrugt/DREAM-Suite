import numpy as np

############################
### Case study 35
############################
def preprocess_data(data, Parameters, n_soil, n_data, interpolation_method, interpolation_scheme, I_max=None):
    """
    This function preprocesses the data and returns a list with final data sets and corresponding true values
    of Ks [cm/h] and S [cm/h^1/2].
    
    Parameters:
        data (list of np.ndarray): List of soil data, each item is a 2-column array of time (hrs) and cumulative infiltration (cm).
        Parameters (list): A list containing the true parameters for each soil type.
        n_soil (int): Number of soil types to process.
        n_data (int): Number of interpolated data points.
        interpolation_method (int): Defines the interpolation method.
        interpolation_scheme (str): The interpolation scheme ('linear', 'logarithmic', or 'square_root').
        I_max (float, optional): Maximum infiltration for normalization (optional).
    
    Returns:
        data_new (list of np.ndarray): Processed data with interpolated time and infiltration values.
        true_pars (np.ndarray): Array with true values of Ks and S for each soil type.
    """
    # Initialize data_new and true_pars
    data_new = [None] * 12
    true_pars = np.full((12, 2), np.nan)

    # Process each soil type
    for soil_type in range(n_soil):
        # Extract time (hours) and cumulative infiltration (cm)
        dat = data[soil_type]
        # Python
        dat = dat[0]
        # Remove NaN values (appear at the end)
        dat = dat[~np.isnan(dat[:, 1]), :2]
        # Remove duplicate values of time
        id = [np.True_] + list(np.diff(dat[:, 0]) > 0)
        # dat = np.vstack([dat[0, :], dat[np.diff(dat[:, 0]) > 0, :]])
        dat = dat[id,:]
        # Remove duplicate values of infiltration
        id = [np.True_] + list(np.diff(dat[:, 1]) > 0)
        # dat = np.vstack([dat[0, :], dat[np.diff(dat[:, 1]) > 0, :]])
        dat = dat[id,:]
        # Extract time and infiltration values
        t_meas = dat[:, 0]
        I_meas = dat[:, 1]
        t_end = t_meas[-1]

        # Interpolation based on method
        if interpolation_method == 1:
            # Divide time interval into n_data + 1 values
            if interpolation_scheme.lower() == 'linear':
                t_int = np.linspace(0, t_end, n_data + 1)
            elif interpolation_scheme.lower() == 'logarithmic':
                t_int = np.logspace(np.log10(1), np.log10(t_end + 1), n_data + 1) - 1
            elif interpolation_scheme.lower() == 'square_root':
                dT = np.sqrt(t_end) / n_data
                t_sqrt_t_int = np.arange(0, np.sqrt(t_end), dT)
                t_int = t_sqrt_t_int**2
                t_int[-1] = t_end
            else:
                raise ValueError("preprocess_data: Unknown interpolation scheme")

            # Interpolate I_meas to t_int
            I_int = np.interp(t_int, t_meas, I_meas)
        elif interpolation_method == 2:
            # Normalize measured infiltration with I_max
            I_norm = I_meas / I_max
            # Interpolation of normalized infiltration
            if interpolation_scheme.lower() == 'linear':
                I_norm_int = np.linspace(0, 1, n_data + 1)
            elif interpolation_scheme.lower() == 'logarithmic':
                I_norm_int = np.logspace(np.log10(1), np.log10(2), n_data + 1) - 1
            else:
                raise ValueError("preprocess_data: Unknown interpolation scheme")

            # Determine corresponding measurement times
            t_int = np.interp(I_norm_int, I_norm, t_meas)
            # Denormalize cumulative infiltration values
            I_int = I_norm_int * I_max
        else:
            raise ValueError("preprocess_data: Unknown interpolation method")

        # Combine (t_int, I_int) and remove time zero (I_int = 0)
        data_new[soil_type] = np.column_stack((t_int[1:], I_int[1:]))

        # Extract true values of Ks [cm/h] and S [cm/h^1/2] from Parameters
        true_pars[soil_type, 0] = Parameters[soil_type + 2][8]  # Ks
        true_pars[soil_type, 1] = Parameters[soil_type + 2][9]  # S

    return data_new, true_pars
