import numpy as np
import subprocess

############################
### Case study 5
############################
def hymodFORTRAN(x):
    # Define local variables
    global F, MaxT

    # Calculate area factor - only once
    if 'F' not in globals():  # Equivalent to checking if persistent variable is empty in MATLAB
        MaxT = 795                                          # Define maximum time
        F = 1944 * (1000 * 1000) / (1000 * 60 * 60 * 24)    # Area factor to translate HYMOD output in mm/d to m3/s 
                                                            # Leaf River area is about 1944 km2

    np.savetxt('Param.in', x, delimiter = ' ')              # Write the parameter values to a file Param.in

    try:    ## Execute the model -- this model reads the current parameter values from Param.in
        result = subprocess.run(['HYMODsilent.exe'], capture_output = True, text = True, creationflags = subprocess.CREATE_NO_WINDOW)
        status = result.returncode
        output = result.stdout
        # Load the output of the model 
        SimRR = F * np.loadtxt('Q.out')
        # Equivalent to SimRR(65:MaxT) in MATLAB
        SimRR = SimRR[64:MaxT]              
    except Exception as e:
        # Return "bad" simulated value if HYMOD FOrtran code did not converge properly
        SimRR = np.zeros(MaxT - 64)  # Equivalent to 0 * [65:MaxT] in MATLAB

    return SimRR.flatten()