import numpy as np


def gaussquad(n):
    """ Computation of weights and nodes of n-point Gaussian quadrature rule
        on the interval [-1, 1]
    """
    if n == 1:
        return 0, 2
    i = np.arange(1, n, dtype='complex')
    b = i / np.sqrt(4 * i**2 - 1)
    J = np.diag(b, -1) + np.diag(b, 1)
    ew, ev = np.linalg.eig(J)
    x = np.real(ew)
    w = 2 * np.real(ev[0, :])**2
    return x, w
