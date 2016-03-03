import numpy as np
from scipy import linalg


def dampnewton(F, DF, x, rtol, atol, lmin=1e-3):
    lu = linalg.lu_factor(DF(x))
    s = linalg.lu_solve(lu, F(x))
    sn = linalg.norm(s)

    xn = x - s
    l = 1.0

    f = F(xn)
    st = linalg.lu_solve(lu, f)
    stn = linalg.norm(st)

    while stn > rtol * linalg.norm(xn) and stn > atol:
        while linalg.norm(st) > (1.0 - l / 2.0) * sn:
            l /= 2.0
            if l < lmin:
                raise ValueError('no convergence: lambda -> 0')
            xn = x - l * s
            f = F(xn)
            st = linalg.lu_solve(lu, f)
        x = xn
        lu = linalg.lu_factor(DF(x))
        s = linalg.lu_solve(lu, f)
        sn = linalg.norm(s)
        l = min(2.0 * l, 1.0)
        xn = x - l * s
        f = F(xn)
        st = linalg.lu_solve(lu, f)
        stn = linalg.norm(st)
    return xn


def example():
    def F(x):
        return np.array([x[0]**2 - 2 * x[0] - x[1] + 1, x[0]**2 + x[1]**2 - 1])

    def DF(x):
        return np.array([[2 * x[0] - 2, -1], [2 * x[0], 2 * x[1]]])

    x = dampnewton(F, DF, np.array([2.0, 3.0]), 1e-6, 1e-8)
    print('||F(x)|| = {}'.format(linalg.norm(F(x))))


if __name__ == '__main__':
    example()
