import numpy as np
from hmodel import hmodel

############################
### Case study 28
############################
def rainfall_runoff(par, plugin):
    # Run hmodel using C++/python source code. The vector par stores the model
    # parameters, and possibly also rainfall multipliers. 
    # See Vrugt et al., 2008. DREAM paper in Water Resources Research

    # Extract values of multipliers
    mult = np.ones(plugin['Tmax'])
    for j in range(plugin['id'].shape[0]):
        mult[plugin['id'][j, 0]:plugin['id'][j, 1] + 1] = par[plugin['nmod'] + j]
    
    # Compute new hyetograph
    plugin['data']['P'] = mult * plugin['P'].copy()

    # Now execute hmodel with hyetograph      
    y = hmodel(par, plugin, 1)

    return y.flatten()


def check_rainfall(plugin):
    # Locates individual storm events in rainfall-discharge record

    # Analyze original rainfall and discharge record
    plugin['id'], P_new = analyze_rainfall(plugin['P'])
    
    # How many multipliers?
    plugin['n_mult'] = plugin['id'].shape[0]
    
    # Now assign the ranges of the multiplier
    mult = {'min': 0.05 * np.ones(plugin['n_mult']), 'max': 5 * np.ones(plugin['n_mult'])}
    
    # Now determine names of multipliers
    mult_names = []
    for ii in range(plugin['n_mult']):
        mult_names.append(f'\\beta_{{{ii + 1}}}')  # 1-based index for MATLAB-style
        
    fpar_mult = np.ones(plugin['n_mult'])
    
    return plugin, mult, mult_names, fpar_mult


def analyze_rainfall(P):
    # Identifies storm events from precipitation data record only
    # Option (commented out) also allows to jointly use discharge as well
    # Indeed: discharge can go up while rainfall is zero --> suspicious

    # Length of P?
    n = P.shape[0]
    
    # Initialize variables
    ct = 0  # Counter for storm events
    flag = 0  # Flag to mark the start of a storm
    id_s = []  # List to store the start indices of storms
    id_e = []  # List to store the end indices of storms
    
    # Loop over each day of precipitation data set
    for j in range(n):
        # Check whether rainfall is larger than zero
        if P[j] > 0:
            if flag == 0:
                # A new storm just started
                flag = 1
                id_s.append(j)      # Record the start index
        elif P[j] == 0 and flag == 1:
            # End of storm event
            id_e.append(j - 1)      # Record the end index
            ct += 1                 # Increment storm counter
            flag = 0                # Reset flag
    
    # Make sure the final idx_end is correct if the last storm didn't end
    if P[n - 1] > 0:
        id_e.append(n - 1)
    
    # Return both start and end indices as a combined matrix
    id = np.array([id_s, id_e]).T  # Convert to numpy array and transpose for consistency
    P_new = P  # Return the original precipitation data (as it's not modified)
    
    return id, P_new


def load_data_dly(ID_watershed):
    """
    This function loads the data of a watershed in dly format.
    
    Args:
    - ID_watershed: str, watershed identifier (used to open the corresponding .dly file)
    
    Returns:
    - data: numpy array, containing the loaded data
    """
    # Construct the file name from the ID_watershed
    filename = f"{ID_watershed}.dly"
    # Try to open the file
    try:
        with open(filename, 'r') as file:
            # Initialize data array (assuming a max of 1e5 rows and 8 columns)
            data = np.full((int(1e5), 8), np.nan)
            # Initialize counter for rows
            counter = 0
            # Read the file line by line
            for line in file:
                try:
                    # Extract year, month, and day from the first part of the line
                    year = int(line[0:4])
                    month = int(line[4:6])
                    day = int(line[6:8])
                    # Extract the remaining data values
                    dat = list(map(float, line[8:].split()))
                    # Store the data (streamflow is assumed to be in the third column)
                    data[counter, :] = [year, month, day] + dat
                    # Increment the counter
                    counter += 1

                except Exception as e:
                    # If there’s an error (e.g., bad formatting), break the loop
                    print(f"Error processing line: {line}")
                    break

            # Trim data to remove excess NaNs (the first NaN row marks the end of valid data)
            data = data[:counter, :]

            return data

    except FileNotFoundError:
        # Handle case when file doesn't exist
        raise ValueError(f"Wrong ID_watershed: {ID_watershed} not found.")
