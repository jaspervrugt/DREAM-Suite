import numpy as np

############################
### Case study 33
############################
def rectangle_lik(x):

    if x.ndim == 1:
        N = 1
        d = len(x)
        x = np.array([x])

    N, d = x.shape  # Get the number of proposals and the dimensionality of each proposal

    if d != 2:
        raise ValueError("Rectangle_lik: Candidate vector should be two-dimensional")

    lik = np.ones(N)  # Initialize likelihood as 1 for all proposals
    # Check if the candidates are within the domain [-0.5, 0.5] x [-3, 3]
    s = (x[:, 0] >= -0.5) & (x[:, 0] <= 0.5) & (x[:, 1] >= -3) & (x[:, 1] <= 3)
    # Set likelihood to 36 for candidates in the rectangle
    lik[s] = 36
    
    return lik
