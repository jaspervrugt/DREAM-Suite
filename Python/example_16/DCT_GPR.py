import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.fft import dct, idct
import matplotlib.pyplot as plt

############################
### Case study 16
############################
def DCT_GPR(x):
    # DCT inversion of geophysical travel time data measured with GPR
    
    # Store local variables in memory
    global func, dummy
    if 'dummy' not in globals():    # Initialize local variables only once
        func = {'error': 0.5,       # Standard deviation of Gaussian traveltime error
                'parx': 8,          # Inversion parameters in x-direction (DCT order)
                'parz': 8}          # Inversion parameters in z-direction (DCT order)
        func = setup_GPR(func)      # Setup forward problem (assuming this function is implemented elsewhere)
        dummy = 1                   # Update dummy

    # -------------------------------------------------------------------------
    #                               Model script
    # -------------------------------------------------------------------------   
    model_DCT = np.zeros((func['dimver'], func['dimhor']))  # Initialize DCT model
    for k in range(func['parz']):                           # Assign proposed model
        for i in range(func['parx']):
            model_DCT[k, i] = x[(k * func['parx']) + i]

    model = idct2(model_DCT)                                # Do inverse DCT
    model = 10 ** model                                     # Transform from logarithmic scale
    model = model.T                                         # Transpose
    minvalue = np.min(model)                                # Find minimum values
    maxvalue = np.max(model)                                # Find maximum values
    data = func['J'] @ model.flatten() - func['datasim']    # Calculate residual
    data = data / func['error']                             # Weighted residual
    data = data.squeeze().reshape(-1,1)
    loglik = data.T @ data                                  # = np.sum(np.power(data, 2))

    if minvalue < (1 / 0.17):                               # Penalty for models outside range
        loglik *= (1 + (1 / 0.17) / minvalue)
    
    if maxvalue > (1 / 0.05):
        loglik *= (1 + maxvalue / (1 / 0.05))
    
    loglik = -0.5 * loglik                                  # Return log-likelihood
    
    return loglik


def setup_GPR(func):
    # Set-up forward kernel and create synthetic data

    x = np.arange(0, 3.1, 0.1)                          # Boundaries of uniform x-grid (3 by 3 m grid); 0.1 m discretization
    z = np.arange(0, 3.1, 0.1)                          # Boundaries of uniform z-grid

    sourcex = 0.01                                      # x-position of source
    sourcez = np.arange(0.05, 3, 0.1)                   # z-positions of sources
    receiverx = 2.99                                    # x-position of receivers
    receiverz = np.arange(0.05, 3, 0.1)                 # z-positions of receivers
    nsource = len(sourcez)
    nreceiver = len(receiverz)

    # Calculate acquisition geometry (multiple-offset gather)
    data = []
    for j in range(nsource):
        for i in range(nreceiver):
            data.append([sourcex, sourcez[j], receiverx, receiverz[i]])

    data = np.array(data)                               # Convert data to numpy array

    # Calculate forward modeling kernel
    func['J'] = tomokernel_straight(data, x, z)         # Distance of ray-segment in each cell for each ray

    # Grid-cells in horizontal and vertical direction
    func['dimhor'] = len(x) - 1
    func['dimver'] = len(z) - 1

    func['por'] = 0.36                                  # Porosity field
    nz = 30                                             # Original discretization of water saturation model
    nx = 40
    wcon = np.loadtxt('Sw.dat')                         # Load water saturation model

    # Create wcont by slicing and reversing the rows of wcon
    wcont = np.zeros((nz, nx - 10))
    for k in range(nz):
        for i in range(6, nx - 5):                      # Make model 3 by 3 m
            wcont[k, i - 5] = wcon[nz - k - 1, i]

    wcont *= func['por']                                # Transform saturation into water content
    wcont = dct(wcont, axis=0, norm='ortho')     # DCT transform along vertical axis (column-wise)

    wtrunc = np.zeros((func['dimver'], func['dimhor']))
    for j in range(func['parz']):
        for i in range(func['parx']):
            wtrunc[j, i] = wcont[j, i]                  # Truncate at the same order as inverse parameterization

    wtrunc = idct2(wtrunc)                              # Inverse DCT along vertical axis
    wtrunc = wtrunc.T                                   # Transpose to match model orientation

    # Translate water content into slowness using Refractive Index/CRIM model
    func['pw'] = 81                                     # Permittivity of water
    func['pa'] = 1                                      # Permittivity of air
    func['ps'] = 5                                      # Permittivity of mineral grains
    slowtrue = (wtrunc.flatten() * np.sqrt(func['pw']) +
                (func['por'] - wtrunc.flatten()) * np.sqrt(func['pa']) +
                (1 - func['por']) * np.sqrt(func['ps']))
    slowtrue /= 0.3                                     # True slowness field

    # Add Gaussian uncorrelated noise with a standard deviation of 1 ns.
    # func['datasim'] = np.dot(func['J'], slowtrue) + func['error'] * np.random.randn(nsource * nreceiver)
    
    # Create the random noise as a sparse matrix (same shape as the result of np.dot)
    func['datasim'] = np.dot(func['J'].todense(), slowtrue) + func['error'] * np.random.randn(nsource * nreceiver)

    return func


def tomokernel_straight(data, x, z):
    """
    Computes the kernel matrix for a straight ray tomographic inversion.
    
    Args:
    - data: A 2D numpy array of shape (nrays, 4), where each row contains:
      [source_x, source_z, receiver_x, receiver_z] of a ray.
    - x: A numpy array representing the x-positions of cell boundaries.
    - z: A numpy array representing the z-positions of cell boundaries.

    Returns:
    - J: The Jacobian matrix in sparse format.
    """
    # Check if data are within bounds set by x and z
    xmin, xmax = x[0], x[-1]
    zmin, zmax = z[0], z[-1]
    
    if xmin > np.min([data[:, 0], data[:, 2]]) or xmax < np.max([data[:, 0], data[:, 2]]) or \
       zmin > np.min([data[:, 1], data[:, 3]]) or zmax < np.max([data[:, 1], data[:, 3]]):
        print("Tomokernel_Straight Error: Data outside of range of min and max values")
        return None
    
    # Determine some initial parameters
    dx = x[1] - x[0]                                        # Horizontal discretization
    dz = z[1] - z[0]                                        # Vertical discretization
    xmid = np.arange(xmin + dx / 2, xmax - dx / 2, dx)      # x-coordinates of cell midpoints
    zmid = np.arange(zmin + dz / 2, zmax - dz / 2, dz)      # z-coordinates of cell midpoints
    nrays = data.shape[0]                                   # Number of rays to consider
    nx = len(x) - 1                                         # Number of cells in x-direction
    nz = len(z) - 1                                         # Number of cells in z-direction
    
    # Initialize sparse matrix storage
    maxelem = round(nrays * np.sqrt(nx**2 + nz**2))
    irow = np.zeros(maxelem, dtype=int)
    icol = np.zeros(maxelem, dtype=int)
    jaco = np.zeros(maxelem)

    # Determine elements of Jacobian matrix
    count = 0
    for i in range(0,1): #nrays):
        xs, zs = data[i, 0], data[i, 1]                     # x-position, z-position of source
        xr, zr = data[i, 2], data[i, 3]                     # x-position, z-position of receiver
        
        if xs == xr:
            xr += 1e-10  # if ray is vertical, add small value for stability
        
        slope = (zr - zs) / (xr - xs)                       # Slope of raypath

        # Determine x-positions of vertical cell boundaries hit by the ray, including ray endpoints
        xcellb = x[(x > min(xs, xr)) & (x < max(xs, xr))]
        xcellb = np.concatenate([xcellb, [xs, xr]])

        # Determine z-positions of horizontal cell boundaries, including ray endpoints
        zcellb = z[(z > min(zs, zr)) & (z < max(zs, zr))]
        zcellb = np.concatenate([zcellb, [zs, zr]])

        # Form matrix containing all intersection points of ray with cell boundaries
        ipoint = np.zeros((len(xcellb)+len(zcellb), 2))
        ipoint[:, 0] = np.concatenate([xcellb, xs + (zcellb - zs) / (slope + 1e-20)])   # x-coordinates
        ipoint[:, 1] = np.concatenate([zs + (xcellb - xs) * slope, zcellb])             # z-coordinates
        # Sort intersection points by x-coordinate
        ipoint = ipoint[np.argsort(ipoint[:, 0])]
        
        # Calculate lengths and midpoints of the ray segments
        xlength = np.abs(ipoint[1:, 0] - ipoint[:-1, 0])    # x-component of length
        zlength = np.abs(ipoint[1:, 1] - ipoint[:-1, 1])    # z-component of length
        clength = np.sqrt(xlength**2 + zlength**2)          # Total length of the ray segment
        cmidpt = 0.5 * (ipoint[:-1, :] + ipoint[1:, :])     # Midpoints of ray segments
        
        # Calculate which slowness cell each ray segment belongs to
        srow = np.ceil((cmidpt[:, 0] - xmin) / dx).astype(int)
        scol = np.ceil((cmidpt[:, 1] - zmin) / dz).astype(int)
        
        srow[srow<1] = 1;  srow[srow>nx] = nx
        scol[scol<1] = 1;  scol[scol>nz] = nz
#        srow = np.clip(srow, 1, nx)                         # Ensure the indices are within bounds
#        scol = np.clip(scol, 1, nz)
        njaco = len(srow)

        # Store values in the sparse matrix arrays
        irow[count:(count + njaco)] = (i+1) * np.ones(njaco, dtype=int)
        icol[count:(count + njaco)] = (scol - 1) * nx + srow
        jaco[count:(count + njaco)] = clength
        
        count += njaco

    # Convert sparse storage arrays to sparse matrix
    index = np.where(jaco > 0)[0]
    irow = irow[index]
    icol = icol[index]
    jaco = jaco[index]
    # Construct sparse matrix J
    J = lil_matrix((nrays, nx * nz))
    J[irow, icol] = jaco

    return J


def GPR_par_ranges(func):
    # Initialize the dictionary to store parameter ranges
    Par_info = {'min': [], 'max': []}

    # Give the parameter ranges (minimum and maximum values)
    Par_info['min'].append(30 * 0.7696)
    Par_info['max'].append(30 * 1.301)      # Corresponds to 1/0.05 to 1/0.17 ns/m in logarithmic units

    # Define grid boundaries
    x = np.arange(0, 3.1, 0.1)              # Boundaries of uniform x-grid (3 by 3 m grid); 0.1 m discretization
    z = np.arange(0, 3.1, 0.1)              # Boundaries of uniform z-grid

    # Grid-cells in horizontal and vertical direction
    func['dimhor'] = len(x) - 1
    func['dimver'] = len(z) - 1

    # Scale DCT coefficients such that all models are possible
    dum = np.zeros((func['dimver'], func['dimhor']))
    count = 1
    for i in range(func['parz']):
        for j in range(func['parx']):
            if i > 0 or j > 0:
                count += 1
                dum[i, j] = 1
                dummy = idct2(dum)
                Par_info['min'].append(-1.7 / np.max(np.abs(dummy)))    # 0.2657
                Par_info['max'].append(1.7 / np.max(np.abs(dummy)))     # 0.2657
                dum[i, j] = 0

    return Par_info


def dct2(arg1, mrows = None, ncols = None):
    """
    2-D Discrete Cosine Transform (DCT).
    
    B = dct2(A) returns the discrete cosine transform of A.
    The matrix B is the same size as A and contains the
    discrete cosine transform coefficients.
    
    B = dct2(A, [M, N]) or B = dct2(A, M, N) pads the matrix A with
    zeros to size M-by-N before transforming. If M or N is smaller than
    the corresponding dimension of A, dct2 truncates A.
    
    This transform can be inverted using idct2.
    
    Parameters:
    -----------
    arg1 : ndarray
        Input matrix A.
    mrows : int, optional
        Number of rows for the transformed matrix (if specified).
    ncols : int, optional
        Number of columns for the transformed matrix (if specified).
    
    Returns:
    --------
    b : ndarray
        The 2-D discrete cosine transformed matrix.
    """
    a = np.array(arg1)
    m, n = a.shape

    if mrows is None:
        mrows = m
    if ncols is None:
        ncols = n

    # Padding for vector input
    mpad, npad = mrows, ncols

    if m == 1 and mpad > m:
        a = np.vstack([a, np.zeros((1, n))])
        m = 2
    if n == 1 and npad > n:
        a = np.hstack([a, np.zeros((m, 1))])
        n = 2
    if m == 1:
        mpad = npad
        npad = 1  # For row vector

    # Perform DCT transformation
    b = dct(a, axis=0, norm='ortho')  # DCT along rows
    b = dct(b, axis=1, norm='ortho')  # DCT along columns

    return b


def idct2(arg1, mrows = None, ncols = None):
    """
    Perform the 2-D Inverse Discrete Cosine Transform.
    
    Parameters:
    arg1 : ndarray
        The input 2D array to be transformed.
    mrows : int, optional
        The number of rows to transform. If not provided, it uses the input shape.
    ncols : int, optional
        The number of columns to transform. If not provided, it uses the input shape.
    
    Returns:
    a : ndarray
        The 2D inverse DCT of the input array.
    """
    m, n = arg1.shape
    
    # If no mrows or ncols are provided, use the original shape of arg1
    if mrows is None:
        mrows = m
    if ncols is None:
        ncols = n

    b = arg1.copy()

    # Padding for vector input
    if m == 1 and mrows > m:
        b = np.pad(b, ((0, 1), (0, 0)), mode = 'constant')
        m = 2
    if n == 1 and ncols > n:
        b = np.pad(b, ((0, 0), (0, 1)), mode = 'constant')
        n = 2
    if m == 1:
        mrows = ncols
        ncols = 1

    # Perform the 2D inverse DCT
    a = idct(b, axis = 0, norm = 'ortho')
    if m > 1 and n > 1:
        a = idct(a, axis = 1, norm = 'ortho')
    
    return a


def GPR_postprocessor(slowtrue, x, z, nx, nz, plugin, chain, output):
# Assuming dd, plugin, and output are defined as in the MATLAB code

    dd = x[1] - x[0]
    nn = len(plugin['datasim'])

    # Plot true model
    plt.figure(3)
    model = np.zeros((30, 30))
    for j in range(nz):
        for i in range(nz):
            model[j, i] = slowtrue[(j - 1) * nz + i]

    model = 1 / model
    for k in range(plugin['dimver']):
        zpatch = np.array([[z[k], z[k] + dd, z[k] + dd, z[k]]]).T
        for i in range(plugin['dimhor']):
            xpatch = np.array([[x[i], x[i], x[i] + dd, x[i] + dd]]).T
            plt.fill(xpatch[:, 0], -zpatch[:, 0], 1000 * model[k, i], linestyle='none')

    plt.xlabel('Distance (m)')
    plt.ylabel('Depth (m)')
    plt.ylim([-3, 0])
    plt.axis('equal')
    plt.title('True model (m/micro.s)')
    plt.clim(50, 170)
    plt.colorbar()

    numb = np.count_nonzero(output['AR'][:,1])
    # Additional figures for realization orders
    for ii in range(4, 9, 2):
        titles = f'Realization order: {ii}'
        plt.figure(ii)
        for kk in range(3):
            for jj in range(3):
                # Plot example model and true model
                model = np.zeros((30, 30))
                for j in range(ii):
                    for i in range(ii):
                        model[j, i] = chain[int(np.round(numb * (0.6 + 0.2 * (jj - 1)))), (j) * plugin['parx'] + i, kk]
                model = np.fft.idct2(model)     # Inverse DCT
                model = 10 ** model             # Slowness field
                model = 1 / model
                # model = (0.3 * model - plugin['por'] * np.sqrt(plugin['pa']) - (1 - plugin['por']) * np.sqrt(plugin['ps'])) / \
                #     (np.sqrt(plugin['pw']) - np.sqrt(plugin['pa']))  # Water content
                
                plt.subplot(3, 3, (kk * 3) + jj + 1)
                for k in range(plugin['dimver']):
                    zpatch = np.array([[z[k], z[k] + dd, z[k] + dd, z[k]]]).T
                    for i in range(plugin['dimhor']):
                        xpatch = np.array([[x[i], x[i], x[i] + dd, x[i] + dd]]).T
                        plt.fill(xpatch[:, 0], -zpatch[:, 0], 1000 * model[k, i], linestyle='none')

                plt.xlabel('Distance (m)')
                plt.ylabel('Depth (m)')
                plt.ylim([-3, 0])
                plt.axis('equal')
                plt.title(titles)
                plt.clim(50, 170)
                plt.colorbar()
