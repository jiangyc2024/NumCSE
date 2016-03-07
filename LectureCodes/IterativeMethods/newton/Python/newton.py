import numpy as np
from scipy import linalg


def newton(F, DFinv, x0, rtol, atol):
    x = np.copy(x0)
    while True:
        s = DFinv(x, F(x))
        x -= s
        sn = linalg.norm(s)
        if sn <= rtol * linalg.norm(x) or sn <= atol:
            return x


def example():
    def F(x):
        return np.array([x[0]**2 - 2 * x[0] - x[1] + 1, x[0]**2 + x[1]**2 - 1])

    def DFinv(x, f):
        J = np.array([[2 * x[0] - 2, -1], [2 * x[0], 2 * x[1]]])
        return linalg.solve(J, f)

    x = newton(F, DFinv, np.array([2.0, 3.0]), 1e-6, 1e-8)
    print('||F(x)|| = {}'.format(linalg.norm(F(x))))


if __name__ == '__main__':
    example()
