import numpy as np


def quadformula(f, c, w):
    """ Generic numerical quadrature routine implementing quadfrom """
    n = len(c)
    assert len(w) == n, '#weight != #nodes'
    s = 0
    for j in range(n):
        s += w[j] * f(c[j])
    return s
