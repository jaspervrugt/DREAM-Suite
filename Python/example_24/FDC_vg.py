import numpy as np

############################
### Case study 24
############################
def FDC_vg(x, plugin):
    # ####################################################################### #
    # FDC_VG flow duration curve according t van Genuchten soil water         #
    # characteristic. Exceedance probability is independent variable          #
    # Please check the work of Sadegh et al., Journal of Hydrology, 2016      #
    # ####################################################################### #

    # Define parameter values
    a = x[0]
    b = x[1]
    c = x[2]
    # Exceedance probability in --> discharge out
    q = (1 / a) * (plugin['E']**(-1 / c) - 1)**(1 / b)

    return q


def calc_FDC(ID, type, col):
    """
    This function calculates the flow duration curve (FDC) for a given watershed.
    
    Args:
    - ID: str or numpy array, watershed ID or already loaded data
    - type: str, type of FDC: 'day', 'week', 'month', 'annual'
    - col: int, column number where discharge data is located (1-indexed in MATLAB, 0-indexed in Python)
    
    Returns:
    - e_n: numpy array, normalized exceedance probabilities for FDC type
    - y_s: numpy array, sorted discharge data for FDC type
    - p_0: float, probability of zero flows for the given FDC type
    """
    if isinstance(ID, (np.ndarray, list)):
        raise ValueError("calc_FDC: First input argument must be a string (file path)")
    elif isinstance(ID, str):
        # Load the data if ID is a string (file path)
        data = load_data_dly(ID)
        # Take one column for Python, JAV
        col = col - 1
        # Take only zero or positive streamflow values
        data = data[data[:, col] >= 0,:]
        # Number of elements in data
        N = data.shape[0]
    else:
        raise ValueError("ID must be a numeric array or a string representing a file path.")
    
    # Initialize probability of zero flows
    p_0 = 0
    Y = []
    # Switch over FDC type
    if type == 'day':  # Daily FDC
        # Probability of zero flows
        p_0 = np.sum(data[:, col] == 0) / N
        # Remove zero values
        Y = data[data[:, col] != 0, col]

    elif type == 'week':  # Weekly FDC
        ct = 0
        p = 0
        for i in range(0, N, 7):  # Iterate over data in blocks of 7 days
            week_data = data[i:min(i + 7, N), col]
            week_data = week_data[week_data > 0]  # Remove non-positive values
            if len(week_data) > 0:
                Y.append(np.mean(week_data))
                ct += 1
            else:
                p += 1
        p_0 = p / (ct + p) if ct + p > 0 else 0.0

    elif type == 'month':  # Monthly FDC
        ct = 0
        p = 0
        id_month = np.concatenate(([0], np.where(np.diff(data[:, 1]) != 0)[0] + 1))  # Find start of each month
        N_month = len(id_month)
        for i in range(1, N_month):
            month_data = data[id_month[i-1]:id_month[i], col]
            month_data = month_data[month_data > 0]
            if len(month_data) > 0:
                Y.append(np.mean(month_data))
                ct += 1
            else:
                p += 1
        p_0 = p / (ct + p) if ct + p > 0 else 0.0

    elif type == 'annual':  # Annual FDC
        ct = 0
        p = 0
        id_year = np.concatenate(([0], np.where(np.diff(data[:, 0]) != 0)[0] + 1))  # Find start of each year
        N_year = len(id_year)
        for i in range(1, N_year):
            year_data = data[id_year[i-1]:id_year[i], col]
            year_data = year_data[year_data > 0]
            if len(year_data) > 0:
                Y.append(np.mean(year_data))
                ct += 1
            else:
                p += 1
        p_0 = p / (ct + p) if ct + p > 0 else 0.0

    else:
        raise ValueError(f"Invalid type '{type}' -- must be 'day', 'week', 'month', or 'annual'")

    # Remove NaNs (if any) due to initialization
    Y = np.array(Y)
    Y = Y[~np.isnan(Y)]
    # Sort discharge data in ascending order
    y_s = np.sort(Y)
    # Exceedance probabilities (based on N only)
    N = len(y_s)
    e_n = 1 - (np.arange(1, N+1) - 0.5) / N

    return e_n, y_s, p_0


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
