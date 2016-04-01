import numpy as np


def intpolyval(t, y, x):
    n = len(t)  # number of interpolation nodes = degree of polynomial - 1
    N = len(x)  # number of evaluation points
    l = np.empty(n)
    for k in range(n):
        l[k] = 1 / np.prod(t[k] - t[np.arange(n) != k])
    p = np.empty(N)
    for i in range(N):
        # Compute quotient of weighted sums
        z = x[i] - t
        j = np.flatnonzero(z == 0)
        if j.size != 0:
            p[i] = y[j]
        else:
            mu = l / z
            p[i] = np.sum(mu * y) / np.sum(mu)
    return p


def main():
    from matplotlib import pyplot as plt

    t = np.linspace(-5, 5, 10)
    y = 1 / (1 + t**2)
    x = np.linspace(-5, 5, 100)

    p = intpolyval(t, y, x)

    plt.plot(t, y, 'o', label='Data')
    plt.plot(x, p, '-', label='Interpolant')
    plt.plot(x, 1 / (1 + x**2), '-', label='Function')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main()
