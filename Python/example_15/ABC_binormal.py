import numpy as np

############################
### Case study 15
############################
def ABC_binormal(x):
    # Persistent variables (stored globally for reuse)
    global n_pairs, sigma, id_1, id_2

    # Load the data and define local variables - only once
    if 'n_pairs' not in globals():
        n_pairs = 10                            # Number of x,y pairs
        sigma = 0.01**2 * np.eye(2)             # Define the sigma of the bivariate distribution
        # Reorganize the x-vector to get bivariate mu values
        id_1 = np.arange(0, n_pairs)            # Indices for x-coordinates
        id_2 = np.arange(n_pairs, 2 * n_pairs)  # Indices for y-coordinates

    # Check whether delta has been specified
    if len(x) == 2 * n_pairs:
        # Set delta to zero if it's not specified
        x = np.hstack([x, 0])

    # Define mu (mean of the bivariate normal distribution)
    mu = np.column_stack([x[id_1], x[id_2]])
    # Loop over each bivariate normal distribution
    X = []
    for i in range(n_pairs):
        # Draw 50 points from each bivariate normal distribution
        mvn_samples = np.random.multivariate_normal(mu[i, :2], sigma, 50)
        # Add a perturbation with std equal to delta [= 0]
        normal_samples = np.random.normal(0, x[-1], (50, 2))
        X.append(np.mean(mvn_samples + normal_samples, axis=0))

    # Return vector with model simulated values
    S = np.array(X).T.flatten()    
    # --> mu[0,0:2] ends up as S[0] and S[n_pairs]
    # --> mu[1,0:2] ends up as S[1] and S[n_pairs+1], etc.
    # --> this matches x coordinates
    # S = np.concatenate(X).flatten()

    return S
