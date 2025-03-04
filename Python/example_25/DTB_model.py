import numpy as np

############################
### Case study 25
############################
def DTB_model(x):
    # ####################################################################### #
    # DTB_model Depth to bedrock landscape evoluton model of Gomez et al.     #
    # ####################################################################### #

    if not hasattr(DTB_model, "initialized"):                      # Store local variables in memory
        # Load input data
        data = np.loadtxt('Input_data.txt')
        # Python - reshape changed as MATLAB counts vertical, Python horizontal!
        DTB_model.Zx = np.reshape(data[1:, 0], (int(data[0, 1]), int(data[0, 0]))).T    # slope gradient x-direction
        # Python - reshape changed as MATLAB counts vertical, Python horizontal!
        DTB_model.Zy = np.reshape(data[1:, 1], (int(data[0, 1]), int(data[0, 0]))).T    # slope gradient y-direction
        # Python - reshape changed as MATLAB counts vertical, Python horizontal!
        Ld = np.reshape(data[1:, 2], (int(data[0, 1]), int(data[0, 0]))).T              # drainage distance
        DTB_model.Ld_n = (Ld - np.min(Ld)) / (np.max(Ld) - np.min(Ld))                  # normalized drainage distance
        DTB_model.Ld = Ld
        DTB_model.nabla_z = DTB_model.Zx ** 2 + DTB_model.Zy ** 2                       # gradient norm
        DTBdata = np.loadtxt('DTB_data.txt')
        id = DTBdata[:, 1].astype(int)                                                  # define training points
        # Turn these indices into valid Python rows and columns
        nr, nc = DTB_model.Ld.shape
        DTB_model.r = np.zeros_like(id)                                                 # Store row indices
        DTB_model.c = np.zeros_like(id)                                                 # Store column indices
        # Loop over the linear indices
        for i in range(len(id)):
            # Column index (0-based in Python)
            DTB_model.c[i] = np.ceil(id[i] / nr) - (1)                                  # Python's 0-based indexing
            # Row index (0-based in Python)
            DTB_model.r[i] = id[i] % nr - (1)                                           # Python 0-based indexing
        DTB_model.initialized = True                                                    # Flag to indicate that initialization is comple

    Lambda = np.exp(-x[1] * (1 - DTB_model.Ld_n) ** x[3])                               # bedrock valley shape term
    Psi = 1 - (DTB_model.nabla_z / x[2] ** 2)                                           # threshold angle of mass movement
    H = (Psi / Lambda) * np.sqrt(x[0] * DTB_model.Ld ** 2)                              # Modeled depth to bedrock
    H = H[DTB_model.r,DTB_model.c]                                                      # Isolate the measurement data points of id
    H = H.flatten()                                                                     # Flatten the matrix to a 1D array

    return H
