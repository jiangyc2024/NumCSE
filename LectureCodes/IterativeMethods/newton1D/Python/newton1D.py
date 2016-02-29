import numpy as np


def newton1D(F, DF, x0, rtol, atol):
    x = x0
    while True:
        s = F(x) / DF(x)
        x -= s
        if abs(s) <= rtol * abs(x) or abs(s) <= atol:
            return x


def example():
    def F(x):
        return np.sin(x) * np.cos(x)

    def DF(x):
        return np.cos(2 * x)

    print('sin(x) * cos(x) has a zero at {}'.format(newton1D(F, DF, 2, 1e-6,
                                                             1e-6)))


if __name__ == '__main__':
    example()
