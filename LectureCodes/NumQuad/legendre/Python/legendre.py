import numpy as np


def legendre(n, x):
    v = np.vstack([np.ones(x.shape), x])
    for j in range(1, n):
        v = np.vstack([v, ((2 * j + 1) / (j + 1.0)) * x * v[-1, :] -
                       j / (j + 1.0) * v[-2, :]])
    return v
