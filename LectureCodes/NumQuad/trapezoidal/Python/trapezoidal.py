import numpy as np


def trapezoidal(f, a, b, n):
    """ Numerical quadrature based on trapezoidal rule """
    x = np.linspace(a, b, n + 1)
    w = np.ones(n + 1)
    w[0] = w[-1] = 0.5
    h = (b - a) / n
    return h * np.dot(w, f(x))


def main():
    def f(x):
        return np.sqrt(1 - x**2)

    print('pi / 2 - Integral of sqrt(1 - x^2) in domain [-1, 1] = {}'.format(
          np.pi / 2 - trapezoidal(f, -1.0, 1.0, 1000)))


if __name__ == '__main__':
    main()
