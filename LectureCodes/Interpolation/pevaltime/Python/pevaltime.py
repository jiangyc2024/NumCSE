import numpy as np
from matplotlib import pyplot as plt
import timeit
import warnings


def anipoleval(t, y, x):
    # scalar-only version
    # y = np.copy(y)
    # for i in range(len(y)):
    #     for k in range(i - 1, -1, -1):
    #         y[k] = y[k + 1] + (y[k + 1] - y[k]) * (x - t[i]) / (t[i] - t[k])
    # return y[0]

    # generic version
    y = np.repeat(y[:, np.newaxis], np.asanyarray(x).size, axis=1)
    for i in range(len(y)):
        for k in range(i - 1, -1, -1):
            y[k, :] = y[k + 1, :] + (y[k + 1, :] - y[k, :]) * ((x - t[i]) /
                                                               (t[i] - t[k]))
    return y[0, :]


def lagrangepoly(x, index, nodes):
    l = 1
    r = np.arange(len(nodes))
    for j in r[r != index]:
        l *= (x - nodes[j]) / (nodes[index] - nodes[j])
    return l


def intpolyval_lag(t, y, x):
    p = np.zeros_like(x)
    for k in range(len(t)):
        p += y[k] * lagrangepoly(x, k, t)
    return p


def intpolyval(t, y, x):
    l = np.empty_like(t)
    p = np.empty_like(x)
    r = np.arange(len(t))
    for k in r:
        l[k] = 1.0 / np.prod(t[k] - t[r != k])
    for i in range(len(x)):
        z = x[i] - t
        j = z == 0
        if np.any(j):
            p[i] = y[j]
        else:
            mu = l / z
            p[i] = np.sum(mu * y) / np.sum(mu)
    return p


def ipoleval(t, y, x):
    p = np.polyfit(t, y, len(y) - 1)
    return np.polyval(p, x)


def main(vector_version=False):
    def run(p, *args):
        return min(timeit.repeat(lambda: p(*args), repeat=10, number=1))
    f = np.sqrt
    res = []
    for n in np.mgrid[3:100]:
        print('Degree = {}'.format(n))
        t = np.mgrid[1:n + 1:1.0]
        y = f(t)
        x = np.random.uniform(high=n, size=n if vector_version else 1)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', np.RankWarning)
            t1 = run(anipoleval, t, y, x)
            t2 = run(ipoleval, t, y, x)
            t3 = run(intpolyval, t, y, x)
            t4 = run(intpolyval_lag, t, y, x)
        res.append((n, t1, t2, t3, t4))
    ns, t1s, t2s, t3s, t4s = np.transpose(res)

    plt.figure()
    plt.semilogy(ns, t1s, '-', label='Aitken-Neville scheme')
    plt.semilogy(ns, t2s, '-', label='Numpy polyfit')
    plt.semilogy(ns, t3s, '-', label='Barycentric formula')
    plt.semilogy(ns, t4s, '-', label='Lagrange polynomials')
    plt.legend()
    plt.xlabel('Polynomial degree')
    plt.ylabel('Computational time')
    plt.title('Polyeval timings (' + ('vector' if vector_version else
                                      'scalar') + 'version)')
    plt.show()


def example_plot():
    t = np.linspace(0, 2 * np.pi, 10)
    y = np.sin(t)
    x = np.linspace(0.1, 1.9 * np.pi, 101)

    plt.figure()
    plt.plot(t, y, 'o')
    plt.plot(x, anipoleval(t, y, x), '-')
    plt.plot(x, ipoleval(t, y, x), '-')
    plt.plot(x, intpolyval(t, y, x), '-')
    plt.plot(x, intpolyval_lag(t, y, x), '-')
    plt.show()


if __name__ == '__main__':
    main(vector_version=False)
    main(vector_version=True)
