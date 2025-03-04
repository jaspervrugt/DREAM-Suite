import numpy as np
import os
import subprocess
import scipy.io


############################
### Case study 7, 21
############################
def HYDRUS(x, data_hydrus):

    # Store local variable (is initially empty after declaration)
    if not hasattr(HYDRUS, "fid"):
        # Write level_01.dir (needed to execute HYDRUS-1D) one time
        with open('level_01.dir', 'w+') as fid:
            fid.write(os.getcwd())      # Write the current directory to the file

        HYDRUS.fid = True               # Mark that the file has been written
    
    # Back-transform log10 sampled VGM parameters (alpha, n and Ks)
    x[2:5] = 10 ** x[2:5]

    try:
        # Run HYDRUS-1D
        sim_hydrus = runH1D(x, data_hydrus)
        # Filter simulated water contents to return values at measurement times
        ind = np.zeros_like(data_hydrus['hoy'], dtype = int)
        for i in range(len(data_hydrus['hoy'])):
            ind[i] = np.where(sim_hydrus['hoy'] == data_hydrus['hoy'][i])[0][0]

        # Return simulated soil moisture contents
        sim_theta = sim_hydrus['water'][ind]
  
    except Exception as e:
        # Return "bad" simulated value if HYDRUS-1D did not converge properly
        # This will result in a very low log-likelihood value
        sim_theta = 1000 * np.ones_like(data_hydrus['hoy'])

    return sim_theta.flatten()


def runH1D(x, data_hydrus):

    # Update the values of "x" in the HYDRUS input file SELECTOR.IN
    modify_selector_in(x)
    # Make a local copy so that it can be modified
    copied_dict = data_hydrus.copy()
    # Modify file PROFILE.DAT using the latest value of bottom boundary condition
    copied_dict['initial'] = np.array(copied_dict['initial'], dtype=float)
    copied_dict['initial'][:, 2] = x[6]
    modify_profile_dat(copied_dict)
    # Modify ATMOSPH.IN with the new boundary condition value
    copied_dict['boundcon'] = np.array(copied_dict['boundcon'], dtype=float)
    copied_dict['boundcon'][:, 6] = x[6]
    modify_atmosph_in(copied_dict)

    # Run HYDRUS-1D using a DOS (or system) command
    try:
        # subprocess.run runs the command and waits for it to finish
        result = subprocess.run(['H1D_CALC.EXE'], capture_output = True, text = True, creationflags = subprocess.CREATE_NO_WINDOW)
        status = result.returncode
        output = result.stdout
    except Exception as e:
        status = -1
        output = str(e)

    # Read the output from OBS_NODE.OUT
    sim_theta = readobsnodeout()

    return sim_theta


def modify_selector_in(x):
    # Modify the SELECTOR.IN file specifying the general setup of HYDRUS

    # Open 'SELECTOR.IN' file in read-write mode
    with open('SELECTOR.IN', 'r+') as file:
        # Read all lines
        lines = file.readlines()
        # Find the line that contains the values to replace
        target_line = '   thr     ths    Alfa      n         Ks       l'
        replacement_line = '   {:.3e}   {:.3e}   {:.3e}   {:.3e}   {:.3e}   {:.3e}\n'.format(*x)
        # Replace the specific line with the new values
        for i, line in enumerate(lines):
            if target_line in line:
                lines[i + 1] = replacement_line

        # Go back to the start of the file to overwrite it with modified lines
        file.seek(0)
        file.truncate()         # Clear the content of the file before writing
        file.writelines(lines)  # Write the modified lines

        
def modify_profile_dat(data_hydrus):
    # Modify the PROFILE.DAT file containing the initial conditions

    # Open the file for writing
    with open('PROFILE.DAT', 'w') as file:
        # Write the header
        file.write('Pcp_File_Version=4\n')
        file.write('3\n')
        file.write('1  0.000000e+000  1.000000e+000  1.000000e+000\n')
        file.write('2 -6.000000e+000  8.000000e+000  8.000000e+000\n')
        file.write('3 -1.000000e+002  8.000000e+001  1.000000e+000\n')
        file.write('81    0    0    1 x         h      Mat  Lay      Beta           Axz            Bxz            Dxz          Temp          Conc\n')      
        # Write the data matrix (data_hydrus['initial'])
        for i in range(len(data_hydrus['initial'])):
            row = data_hydrus['initial'][i]
            file.write(f"{int(row[0]):5} {row[1]:15.6e} {row[2]:15.6e} {int(row[3]):3} {int(row[4]):3} {row[5]:15.6e} {row[6]:15.6e} {row[7]:15.6e} {row[8]:15.6e}\n")
        
        # Write the three additional lines at the bottom
        file.write("\n")
        file.write("    1\n")
        file.write("   32\n")


def modify_atmosph_in(data_hydrus):
    # Modify the ATMOSPH.IN file containing atmospheric boundary condition

    # Open ATMOSPH.IN in read-write mode
    with open('ATMOSPH.IN', 'w') as file:
        # Write the header
        file.write('Pcp_File_Version=4\n')
        file.write('*** BLOCK I: ATMOSPHERIC INFORMATION  **********************************\n')
        file.write('MaxAL                    (MaxAL = number of atmospheric data-records)\n')
        file.write('6840\n')
        file.write('DailyVar  SinusVar  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy  lDummy\n')
        file.write('    f       f       f       f       f       t       f       f       f       f\n')
        file.write('hCritS                 (max. allowed pressure head at the soil surface)\n')
        file.write('    1\n')
        file.write('    tAtm        Prec       rSoil       rRoot      hCritA          rB          hB          ht    RootDepth\n')

        # Write the data matrix (data_hydrus['boundcon'])
        for i in range(len(data_hydrus['boundcon'])):
            row = data_hydrus['boundcon'][i]
            file.write(f"{int(row[0]):12} {row[1]:10.2e} {row[2]:17.4e} {int(row[3]):7} {int(row[4]):17} {int(row[5]):6} {row[6]:18.2e} {int(row[7]):2}\n")

        file.write('end*** END OF INPUT FILE ''ATMOSPH.IN'' **********************************\n')


def readobsnodeout():
    # Read the OBS_NODE.OUT file containing the time series of simulated soil water contents

    # Open OBS_NODE.OUT in read mode
    with open('OBS_NODE.OUT', 'r') as file:
        lines = file.readlines()

    # Skip lines the first 11 lines and do not read the last line as well
    lines_to_read = lines[11:-1]
    # Remove empty lines or lines with only whitespace
    lines_to_read = [line for line in lines_to_read if line.strip()]
    # Split each line into columns and ensure they're in the correct format
    clean_lines = []
    for line in lines_to_read:
        # Split the line by spaces and remove extra spaces
        clean_lines.append(' '.join(line.split()))

    # Now use np.loadtxt to read the cleaned lines
    data = np.loadtxt(clean_lines, delimiter = ' ')
    # Transpose the data to match the format in the MATLAB version
    data = data.T
    # Store simulated soil water contents in "sim_hydrus"
    sim_hydrus = {
        'hoy': data[0],   # 'hoy' corresponds to the first column (time)
        'water': data[2]  # 'water' corresponds to the third column (simulated water content)
    }

    return sim_hydrus


def Load_data():
    """
    Load observational data for Hydrus and define initial and boundary conditions.
    Returns a dictionary with the required data_hydrus.
    """
    # Load 'meas.mat'
    data = scipy.io.loadmat('meas.mat')
    meas = data['meas']
    data_hydrus = {}
    flat_waterind = meas['waterind'][0][0].flatten()
    flat_water = meas['water'][0][0]
    flat_hoy = meas['hoy'][0][0]
    data_hydrus['water'] = flat_water[flat_waterind==1]
    data_hydrus['hoy'] = flat_hoy[flat_waterind==1]
    data_hydrus['ObsNode'] = 1
    # Load 'ProfileDat.mat'
    profile_data = scipy.io.loadmat('ProfileDat.mat')
    data_hydrus['initial'] = profile_data['initial_conditions']  # Assuming it's an array or matrix
    # Load 'AtmosphIn.mat'
    atmosph_data = scipy.io.loadmat('AtmosphIn.mat')
    data_hydrus['boundcon'] = atmosph_data['boundcon']  # Assuming it's an array or matrix
    
    return data_hydrus