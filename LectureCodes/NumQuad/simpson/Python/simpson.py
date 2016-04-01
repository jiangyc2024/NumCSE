import numpy as np


def simpson(f, a, b, n):
    """ Numerical quadrature using Simpson rule """
    x = np.linspace(a, b, 2 * n + 1)
    h = (b - a) / n
    fv = f(x)
    return np.sum(h * (fv[:-2:2] + 4 * fv[1::2] + fv[2::2])) / 6


def main():
    def f(x):
        return np.sqrt(1 - x**2)

    print('pi / 2 - Integral of sqrt(1 - x^2) in domain [-1, 1] = {}'.format(
          np.pi / 2 - simpson(f, -1.0, 1.0, 1000)))


if __name__ == '__main__':
    main()
