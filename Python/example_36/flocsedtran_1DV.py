import numpy as np

############################
### Case study 36
############################
def flocsedtran_1DV(par, plugin):
    """
    2-parameter flocculation model
    
    Args:
        par (array): Array of parameters.
            par[0] = M (Empirical erosion parameter) (kg/m2/s)
            par[1] = τ_c (Critical shear stress for erosion) (Pa)
        plugin (dict): Dictionary of model inputs.
        
    Returns:
        np.ndarray: Simulated temperature data as a single vector.
    """
    # Unpack parameter values
    prFloc = plugin['prFloc']
    prFloc['alpha_pbe1'] = par[0]
    prFloc['alpha_pbe2'] = par[0]
    prFloc['alpha_pbe3'] = par[0]
    prFloc['frac_df'] = par[1]
    np_grid = plugin['prGrid']['ngp']
    
    # Initial states
    ts_onemab = plugin['ts_onemab']
    sig_b = plugin['sig_b']
    time = plugin['time']
    
    # Execute Initial Conditions
    xn1, xn2, xn3, xnt, conc, fdia, xnc, floc_den, wsf, wsp, pri_tmass, tde1, te, shear = InitialConditions(time, plugin['prIni'], prFloc, plugin['prGrid'], plugin['prFlow'], 0)

    for itm in range(32400):  # Time loop
        time += plugin['prTime']['dt']
        psig_b = sig_b
        
        # Update Flow Fields
        ustar, tau_b, tde1, te, shear = updataFlowfield(time, prFloc, plugin['prGrid'], plugin['prFlow'], tde1, te, shear)
        
        # Explicit solver for tcpbe
        xn1, xn2, xn3 = explicitTCPBEcaculation(plugin['beta'], np_grid, prFloc, plugin['prTime'], xnc, fdia, xn1, xn2, xn3, shear, wsp, wsf)
        
        # GS iteration for solving implicit 1st order sediment transport equation
        if np_grid == 1:
            xn1, xn2, xn3, xnt, conc, fdia, xnc, floc_den, wsp, wsf, ers_np, dep_nt = depthaveSedtran(np_grid, tau_b, plugin['delz'], prFloc, plugin['prTime'], plugin['prBot'], xn1, xn2, xn3, xnt, conc, fdia, xnc, floc_den, tde1, wsp, wsf)
        else:
            xn1, xn2, xn3, xnt, conc, fdia, xnc, floc_den, wsp, wsf, ers_np, dep_nt = implicitSedTransCalculation(plugin['ae2'], plugin['ae3'], np_grid, tau_b, plugin['delz'], prFloc, plugin['prTime'], plugin['prBot'], xn1, xn2, xn3, xnt, conc, fdia, xnc, floc_den, tde1, wsp, wsf)
        
        # Update thickness of bottom layer (bl)
        sig_b = psig_b + plugin['prTime']['dt'] * plugin['prBot']['fntc'] / prFloc['cgel'] * (-ers_np - dep_nt)
        
        # Post processing at each time
        if np.any((plugin['ts_obs'] >= time - plugin['prTime']['dt']) & (plugin['ts_obs'] < time)):
            onemab_pfrac = xn1[0] / (xn1[0] + xn3[0])
            onemab_fdia = fdia[0] * 1.0e3       # millimeter
            onemab_conc = conc[0]
            timeHour = time / 3600
            onemab_combined = np.column_stack((timeHour, ustar, onemab_pfrac, onemab_fdia, onemab_conc))
            ts_onemab = np.vstack((ts_onemab, onemab_combined))

    # Post processing at the end of the time loop
    ts_onemab = ts_onemab[1:, 2:]   # Remove the first row
    # floc_sim = ts_onemab.flatten()  # Return as a single vector
    floc_sim = np.concatenate((ts_onemab[:,0], ts_onemab[:,1], ts_onemab[:,2])) 

    return floc_sim


def InitialConditions(time, prIni, prFloc, prGrid, prFlow, pri_tmass):
    """
    Function defines the initial conditions for sediment transport.

    Args:
        time (float): Current time (s).
        prIni (dict): Initial conditions parameters.
        prFloc (dict): Flocculation parameters.
        prGrid (dict): Grid parameters.
        prFlow (dict): Flow parameters.
        pri_tmass (float): Previous total mass.

    Returns:
        tuple: Arrays of initial conditions including xn1, xn2, xn3, etc.
    """
    np_grid = prGrid['ngp']

    # Preallocate arrays
    xn1 = np.zeros(np_grid)
    xn2 = np.zeros(np_grid)
    xn3 = np.zeros(np_grid)
    xnt = np.zeros(np_grid)
    conc = np.zeros(np_grid)
    fdia = np.zeros(np_grid)
    xnc = np.zeros(np_grid)
    floc_den = np.zeros(np_grid)
    wsp = np.zeros(np_grid)
    wsf = np.zeros(np_grid)
    tde1 = np.zeros(np_grid)
    te = np.zeros(np_grid)
    shear = np.zeros(np_grid)

    # ustar0 computation
    ustar0 = abs(0.036 + 0.000 * np.sin(2 * np.pi / 12.4 * (time / 3600) + 0.5) - 
                 0.020 * np.cos(2 * (2 * np.pi / 12.4 * (time / 3600) + 0.5)))

    # Loop for each grid point
    for k in range(0, np_grid):
        # Sediment and floc profile
        xn1[k] = prIni['pnconc'] * 0.90
        xn3[k] = prIni['pnconc'] * 0.10
        fdia[k] = prFloc['ffdia']
        xnc[k] = (fdia[k] / prFloc['pdia']) ** prFloc['frac_df']
        xn2[k] = xn3[k] / xnc[k]
        xnt[k] = xn1[k] + xn3[k]
        conc[k] = xnt[k] * prFloc['par_den'] * prFloc['pri_vol']
        floc_den[k] = prFloc['wat_den'] + (prFloc['par_den'] - prFloc['wat_den']) * \
                      (prFloc['pdia'] / fdia[k]) ** (3 - prFloc['frac_df'])
        re_p = fdia[k] * -wsf[k] * prFloc['wat_den'] / prFloc['vmu']
        wsf[k] = -(prFloc['par_den'] - prFloc['wat_den']) * prFloc['g'] * \
                 prFloc['pdia'] ** (3 - prFloc['frac_df']) * fdia[k] ** (prFloc['frac_df'] - 1) / \
                 (18.0 * prFloc['vmu']) / (1 + 0.15 * re_p ** 0.687)
        wsp[k] = -(prFloc['pdia'] ** 2) * prFloc['g'] * (prFloc['par_den'] - prFloc['wat_den']) / \
                 (18.0 * prFloc['vmu'])
        pri_tmass += xnt[k] * prGrid['dz']

        # Flow and turbulence profile
        tde1[k] = prFlow['v_karm'] * ustar0 * (k + 1) * prGrid['dz'] * (1.0 - (k + 1) * prGrid['dz'] / prGrid['t_height'])
        te[k] = ustar0 ** 3 / (prFlow['v_karm'] * prGrid['dz'] * (k + 1))
        shear[k] = np.sqrt(te[k] / (prFloc['vmu'] * 1.0e-3))

    return xn1, xn2, xn3, xnt, conc, fdia, xnc, floc_den, wsf, wsp, pri_tmass, tde1, te, shear


def updataFlowfield(time, prFloc, prGrid, prFlow, tde1, te, shear):
    """
    Update the flow field parameters like shear, ustar, eddy viscosity, etc.

    Args:
        time (float): Current time in seconds.
        prFloc (dict): Flocculation parameters.
        prGrid (dict): Grid parameters.
        prFlow (dict): Flow parameters.
        tde1 (numpy.ndarray): Eddy viscosity related to the bottom shear.
        te (numpy.ndarray): Turbulent energy.
        shear (numpy.ndarray): Shear stress.

    Returns:
        tuple: Updated values for ustar, tau_b, tde1, te, and shear.
    """
    np_grid = prGrid['ngp']

    # Calculate bottom shear (ustar) and shear rate profile
    ustar = abs(0.036 + 0.000 * np.sin(2.0 * np.pi / 12.4 * (time / 3600) + 0.5) -
                0.020 * np.cos(2.0 * (2.0 * np.pi / 12.4 * (time / 3600) + 0.5)))
    tau_b = prFloc['wat_den'] * ustar**2

    # Update eddy viscosity and turbulence-related parameters
    if np_grid == 1:  # Depth-averaged Model
        k = 1
        tde1[k-1] = prFlow['v_karm'] * ustar * float(k) * (prGrid['dz'] / 2.0) * \
                    (1.0 - float(k) * (prGrid['dz'] / 2) / prGrid['t_height'])
        te[k-1] = ustar**3 / (prFlow['v_karm'] * (prGrid['dz'] / 2) * float(k))
        shear[k-1] = np.sqrt(te[k-1] / (prFloc['vmu'] * 1.0e-3))
    else:
        for k in range(np_grid):
            tde1[k] = prFlow['v_karm'] * ustar * float(k+1) * prGrid['dz'] * \
                      (1.0 - float(k+1) * prGrid['dz'] / prGrid['t_height'])
            te[k] = ustar**3 / (prFlow['v_karm'] * prGrid['dz'] * float(k+1))
            shear[k] = np.sqrt(te[k] / (prFloc['vmu'] * 1.0e-3))

    return ustar, tau_b, tde1, te, shear


def explicitTCPBEcaculation(beta, ngp, prFloc, prTime, xnc, fdia, xn1, xn2, xn3, shear, wsp, wsf):
    """
    Perform the explicit TCPBE calculation.

    Args:
        beta (dict): Dictionary containing `bm1`, `sh1`, `pbe1`, etc.
        ngp (int): Number of grid points.
        prFloc (dict): Flocculation parameters.
        prTime (dict): Time-related parameters, including `dt`.
        xnc (numpy.ndarray): Array for particle concentration.
        fdia (numpy.ndarray): Array for particle diameters.
        xn1 (numpy.ndarray): Array for first scalar concentration.
        xn2 (numpy.ndarray): Array for second scalar concentration.
        xn3 (numpy.ndarray): Array for third scalar concentration.
        shear (numpy.ndarray): Array for shear stress values.
        wsp (numpy.ndarray): Array for settling velocity in pure water.
        wsf (numpy.ndarray): Array for sediment settling velocity.

    Returns:
        tuple: Updated arrays for `xn1`, `xn2`, `xn3`.
    """
    # Initialize an array for brk_s
    brk_s = np.zeros(ngp)

    # Copy the values of xn1, xn2, xn3 for previous scalar concentrations
    ppxn1 = xn1.copy()
    ppxn2 = xn2.copy()
    ppxn3 = xn3.copy()

    # Loop through each particle in the grid
    for k in range(ngp):
        xnc[k] = xn3[k] / xn2[k]  # Update xnc
        fdia[k] = xnc[k] ** (1.0 / prFloc['frac_df']) * prFloc['pdia']  # Update fdia

        # Compute the beta values for different parameters
        beta['bm1'][k] = 2.0 * prFloc['bolz'] * prFloc['temp'] / (3.0 * prFloc['vmu']) * 4.0
        beta['sh1'][k] = (1.0 / 6.0) * shear[k] * (2.0 * prFloc['pdia']) ** 3.0
        beta['pbe1'][k] = prFloc['alpha_pbe1'] * (beta['bm1'][k] + beta['sh1'][k])

        beta['bm2'][k] = 2.0 * prFloc['bolz'] * prFloc['temp'] / (3.0 * prFloc['vmu']) * \
                         (1.0 / prFloc['pdia'] + 1.0 / fdia[k]) * (prFloc['pdia'] + fdia[k])
        beta['ds2'][k] = np.pi / 4.0 * (fdia[k] + prFloc['pdia']) ** 2.0 * abs(wsf[k] - wsp[k])
        beta['sh2'][k] = (1.0 / 6.0) * shear[k] * (prFloc['pdia'] + fdia[k]) ** 3.0
        beta['pbe2'][k] = prFloc['alpha_pbe2'] * (beta['bm2'][k] + beta['ds2'][k] + beta['sh2'][k])

        beta['bm3'][k] = 2.0 * prFloc['bolz'] * prFloc['temp'] / (3.0 * prFloc['vmu']) * 4.0
        beta['sh3'][k] = (1.0 / 6.0) * shear[k] * (fdia[k] + fdia[k]) ** 3.0
        beta['pbe3'][k] = prFloc['alpha_pbe3'] * (beta['bm3'][k] + beta['ds3'][k] + beta['sh3'][k])

        # Calculate breakage term brk_s
        brk_s[k] = prFloc['es'] * shear[k] * (fdia[k] / prFloc['pdia'] - 1.0) ** (3.0 - prFloc['frac_df']) * \
                   (prFloc['vmu'] * shear[k] * fdia[k] ** 2.0 * prFloc['ffy']) ** 1.00

        # Update xn1, xn2, xn3 based on the explicit formula
        xn1[k] = ppxn1[k] + prTime['dt'] * (-0.5 * beta['pbe1'][k] * ppxn1[k] ** 2 * (xnc[k] / (xnc[k] - 1.0)) -
                                            beta['pbe2'][k] * ppxn1[k] * ppxn2[k] +
                                            prFloc['brk_f'] * xnc[k] * brk_s[k] * ppxn2[k])

        xn2[k] = ppxn2[k] + prTime['dt'] * (0.5 * beta['pbe1'][k] * ppxn1[k] ** 2 * (1.0 / (xnc[k] - 1.0)) -
                                            0.5 * beta['pbe3'][k] * ppxn2[k] ** 2 +
                                            brk_s[k] * ppxn2[k])

        xn3[k] = -xn1[k] + ppxn1[k] + ppxn3[k]

    return xn1, xn2, xn3


def implicitSedTransCalculation(ae2, ae3, ngp, tau_b, delz, prFloc, prTime, prBot, 
                                xn1, xn2, xn3, xnt, conc, fdia, xnc, floc_den, tde1, wsp, wsf):
    """
    Implicit Sediment Transport Calculation.

    Args:
        ae2 (numpy.ndarray): Array for second diffusion coefficient.
        ae3 (numpy.ndarray): Array for third diffusion coefficient.
        ngp (int): Number of grid points.
        tau_b (float): Bottom shear stress.
        delz (numpy.ndarray): Array of depth intervals.
        prFloc (dict): Dictionary with flocculation-related parameters.
        prTime (dict): Dictionary with time-related parameters.
        prBot (dict): Dictionary with bottom boundary parameters.
        xn1 (numpy.ndarray): Array for first scalar concentration.
        xn2 (numpy.ndarray): Array for second scalar concentration.
        xn3 (numpy.ndarray): Array for third scalar concentration.
        xnt (numpy.ndarray): Array for total scalar concentration.
        conc (numpy.ndarray): Array for concentration.
        fdia (numpy.ndarray): Array for particle diameters.
        xnc (numpy.ndarray): Array for particle concentration.
        floc_den (numpy.ndarray): Array for floc density.
        tde1 (numpy.ndarray): Array for turbulent diffusion coefficients.
        wsp (numpy.ndarray): Array for settling velocity in pure water.
        wsf (numpy.ndarray): Array for sediment settling velocity.

    Returns:
        tuple: Updated arrays (xn1, xn2, xn3, xnt, conc, fdia, xnc, floc_den, wsp, wsf, ers_np, dep_nt).
    """
    # Initialization for sediment transport
    change_max = 1.0
    j = 0

    pxn1 = xn1.copy()
    pxn2 = xn2.copy()
    pxn3 = xn3.copy()

    # Deposition/Erosion (sink/source) calculation - Le Hir (2011, csr)
    if tau_b > prBot['tau_c']:
        ers_np = prBot['ers_m'] / (1.0 / 6.0 * np.pi * (prFloc['pdia']**3) * prFloc['par_den']) * \
                 (tau_b / prBot['tau_c'] - 1.0)
    else:
        ers_np = 0.0

    dep_np = wsp[0] * xn1[0]
    dep_nf = wsf[0] * xn2[0]
    dep_nt = wsf[0] * xn3[0]
    ers_np += dep_np

    while change_max > 1.0e-10:
        j += 1
        cxn3 = xn3.copy()  # Just for checking numerical solution

        # Update mass and volume concentration
        for k in range(ngp):
            xnc[k] = xn3[k] / xn2[k]
            fdia[k] = xnc[k] ** (1.0 / prFloc['frac_df']) * prFloc['pdia']
            xnt[k] = xn1[k] + xn3[k]
            conc[k] = xnt[k] * prFloc['par_den'] * prFloc['pri_vol']
            floc_den[k] = prFloc['wat_den'] + (prFloc['par_den'] - prFloc['wat_den']) * \
                          (prFloc['pdia'] / fdia[k]) ** (3.0 - prFloc['frac_df'])
            wsf[k] = -0.055556 * prFloc['g'] * (floc_den[k] / prFloc['wat_den'] - 1.0) * \
                     fdia[k]**2 / (prFloc['vmu'] * 1.0e-3)
            wsp[k] = -0.055556 * prFloc['g'] * (prFloc['par_den'] / prFloc['wat_den'] - 1.0) * \
                     prFloc['pdia']**2 / (prFloc['vmu'] * 1.0e-3)

        # Diffusion coefficient for 2nd order approximation w.r.t. time
        for k in range(ngp):
            if k == 0:
                ae2[k] = prTime['dt'] * (tde1[k + 1] + tde1[k]) / (2.0 * delz[k]**2)
            elif k == ngp - 1:
                ae3[k] = prTime['dt'] * (tde1[k] + tde1[k - 1]) / (2.0 * delz[k - 1]**2)
            else:
                ae2[k] = prTime['dt'] * (tde1[k + 1] + tde1[k]) / (2.0 * delz[k]**2)
                ae3[k] = prTime['dt'] * (tde1[k] + tde1[k - 1]) / (2.0 * delz[k - 1]**2)

        # Scalar (sediment concentration) update inside the nodal system
        for k in range(ngp):
            if k == 0:  # Node calculation at the bottom boundary (open)
                xn1[k] = (pxn1[k] - prTime['dt'] / delz[k] * (wsp[k + 1] * xn1[k + 1] - ers_np) + ae2[k] * xn1[k + 1]) / \
                          (1.0 + ae2[k])
                xn2[k] = (pxn2[k] - prTime['dt'] / delz[k] * (wsf[k + 1] * xn2[k + 1] - dep_nf) + ae2[k] * xn2[k + 1]) / \
                          (1.0 + ae2[k])
                xn3[k] = (pxn3[k] - prTime['dt'] / delz[k] * (wsf[k + 1] * xn3[k + 1] - dep_nt) + ae2[k] * xn3[k + 1]) / \
                          (1.0 + ae2[k])

            elif k == ngp - 1:  # Node calculation at the top boundary (closed)
                xn1[k] = (pxn1[k] + ae3[k] * xn1[k - 1]) / \
                          (1.0 - prTime['dt'] / delz[k] * wsp[k] + ae3[k])
                xn2[k] = (pxn2[k] + ae3[k] * xn2[k - 1]) / \
                          (1.0 - prTime['dt'] / delz[k] * wsf[k] + ae3[k])
                xn3[k] = (pxn3[k] + ae3[k] * xn3[k - 1]) / \
                          (1.0 - prTime['dt'] / delz[k] * wsf[k] + ae3[k])

            else:  # Node calculation inside boundary
                xn1[k] = (pxn1[k] - prTime['dt'] / delz[k] * wsp[k + 1] * xn1[k + 1] + ae2[k] * xn1[k + 1] + ae3[k] * xn1[k - 1]) / \
                          (1.0 - prTime['dt'] / delz[k] * wsp[k] + ae2[k] + ae3[k])
                xn2[k] = (pxn2[k] - prTime['dt'] / delz[k] * wsf[k + 1] * xn2[k + 1] + ae2[k] * xn2[k + 1] + ae3[k] * xn2[k - 1]) / \
                          (1.0 - prTime['dt'] / delz[k] * wsf[k] + ae2[k] + ae3[k])
                xn3[k] = (pxn3[k] - prTime['dt'] / delz[k] * wsf[k + 1] * xn3[k + 1] + ae2[k] * xn3[k + 1] + ae3[k] * xn3[k - 1]) / \
                          (1.0 - prTime['dt'] / delz[k] * wsf[k] + ae2[k] + ae3[k])

        # Check convergence
        change = np.abs((cxn3 - xn3) / xn3)
        change_max = np.max(change)

    return xn1, xn2, xn3, xnt, conc, fdia, xnc, floc_den, wsp, wsf, ers_np, dep_nt


def model_properties():
    """
    Initialize model properties for sediment transport and flocculation.
    Returns:
        dict: A dictionary containing various model properties.
    """
    # Initialize dictionaries for the model properties
    prFloc = {}
    prIni = {}
    prBot = {}
    prFlow = {}
    prGrid = {}
    prTime = {}
    plugin = {}

    # Aggregation/Break-up Kinetic Constants
    prFloc['alpha_pbe1'] = 0.10
    prFloc['alpha_pbe2'] = 0.10
    prFloc['alpha_pbe3'] = 0.10

    # Correction Coefficient for Diff Settling-Mediated Flocculation
    prFloc['es'] = 2.0e-5
    prFloc['ffy'] = 1.0e+10
    prFloc['brk_f'] = 1.0

    # Constant for Fractal Theory - Check Floc et al (AICHE, 1999)
    prFloc['frac_df'] = 2.3

    # Physicochemical Properties of Solid and Liquid
    prFloc['par_den'] = 1.80e3
    prFloc['wat_den'] = 1.05e3
    prFloc['g'] = 9.81
    prFloc['vmu'] = 1.002e-3
    prFloc['bolz'] = 1.38e-23
    prFloc['temp'] = 293.0

    # Constants (make 'MKS' units) - Check Spicer & Pratsinis (AICHE, 1996)
    prFloc['pdia'] = 15.0e-6
    prFloc['ffdia'] = 100.0e-6
    prFloc['cgel'] = 100.0
    prFloc['pri_vol'] = (1.0 / 6.0) * np.pi * prFloc['pdia']**3

    # Initial Concentration
    prIni['pconc'] = 0.32
    prIni['pnconc'] = prIni['pconc'] / prFloc['par_den'] / prFloc['pri_vol']

    # Erosion/Deposition Constants
    prBot['ers_m'] = 0.25e-3
    prBot['tau_c'] = 1.0
    prBot['fntc'] = (1.0 / 6.0) * np.pi * prFloc['pdia']**3 * prFloc['par_den']

    # Coefficient for boundary conditions for velocity and turbulent equations
    prFlow['v_karm'] = 0.4

    # Grid Information
    prGrid['t_height'] = 7.0
    prGrid['dz'] = 1.0
    prGrid['ngp'] = int(prGrid['t_height'] / prGrid['dz'])

    # Time step size
    prTime['dt'] = 10.0
    printHour = 0.5  # hours
    prTime['print'] = int(printHour * 3600 / prTime['dt'])

    # Save in plugin structure
    plugin['prFloc'] = prFloc
    plugin['prIni'] = prIni
    plugin['prBot'] = prBot
    plugin['prFlow'] = prFlow
    plugin['prGrid'] = prGrid
    plugin['prTime'] = prTime

    # Initialize beta for explicit solution
    plugin['beta'] = {
        'pbe1': np.zeros(prGrid['ngp']),
        'bm1': np.zeros(prGrid['ngp']),
        'sh1': np.zeros(prGrid['ngp']),
        'pbe2': np.zeros(prGrid['ngp']),
        'bm2': np.zeros(prGrid['ngp']),
        'ds2': np.zeros(prGrid['ngp']),
        'sh2': np.zeros(prGrid['ngp']),
        'pbe3': np.zeros(prGrid['ngp']),
        'bm3': np.zeros(prGrid['ngp']),
        'ds3': np.zeros(prGrid['ngp']),
        'sh3': np.zeros(prGrid['ngp'])
    }

    # Initialize ae2 and ae3 for implicit solution
    plugin['ae2'] = np.zeros(prGrid['ngp'])
    plugin['ae3'] = np.zeros(prGrid['ngp'])

    # Grid info (delz of settling column)
    plugin['delz'] = np.ones(prGrid['ngp']) * prGrid['dz']

    plugin['ngp'] = prGrid['ngp']       # Number of nodes
    plugin['ts_onemab'] = np.zeros(5)   # Initial condition (initial seed)
    plugin['sig_b'] = 0.0               # Initial sig_b
    plugin['time'] = 0.0                # Initial time

    return plugin
