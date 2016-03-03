import numpy as np


def bisect(f, a, b, tol):
    """ Searching zero of f in [a, b] by bisection """
    if a > b:
        a, b = b, a
    fa = f(a)
    fb = f(b)
    if fa * fb > 0:
        raise ValueError('f(a) and f(b) have the same sign')
    v = -1 if fa > 0 else 1
    x = 0.5 * (a + b)
    while b - a > tol and a < x < b:
        if v * f(x) > 0:
            b = x
        else:
            a = x
        x = 0.5 * (a + b)
    return x


def example():
    def f(x): return np.sin(x) * np.cos(x)
    print('sin(x) * cos(x) has zero at {}'.format(bisect(f, 1, 2, 1e-6)))


if __name__ == '__main__':
    example()
