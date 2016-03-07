import numpy as np


def secant(x0, x1, F, rtol, atol, maxit):
    fo = F(x0)
    for i in range(maxit):
        fn = F(x1)
        s = fn * (x1 - x0) / (fn - fo)
        x0 = x1
        x1 = x1 - s
        if abs(s) < max(atol, rtol * min(abs(x0), abs(x1))):
            return x1
        fo = fn


def example():
    def F(x):
        return np.sin(x) * np.cos(x)

    print('sin(x) * cos(x) has a zero at {}'.format(secant(1, 2, F, 1e-6, 1e-6,
                                                           100)))


if __name__ == '__main__':
    example()
