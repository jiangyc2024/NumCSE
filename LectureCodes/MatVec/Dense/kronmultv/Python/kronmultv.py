import numpy as np


def kron_a(A, B, x):
    """ Computes kron(A, B) * x (matrix version). """
    M = np.kron(A, B)
    y = np.dot(M, x)
    return y


def kron_b(A, B, x):
    """ Computes kron(A, B) * x (smart version). """
    m, n = A.shape
    l, k = B.shape
    assert x.size == n * k, 'size mismatch'

    # kron gives a matrix with n x n blocks, block i, j is A[i, j] * B
    # => y = M * x can be done block-wise so that we resue B * x[...]

    # init
    y = np.zeros(m * l)
    # loop first over columns and then (!) over rows
    for j in range(n):
        # reuse B * x[...] part (constant in given column) => O(n^2)
        z = np.dot(B, x[j * k:(j + 1) * k])
        # add to result vector (need to go through full vector) => O(n^2)
        for i in range(m):
            y[i * l:(i + 1) * l] += A[i, j] * z
    return y


def kron_c(A, B, x):
    """ Computes kron(A, B) * x (smart version with reshapes). """
    m, n = A.shape
    l, k = B.shape
    assert x.size == n * k, 'size mismatch'

    # init
    yy = np.zeros((m, l))
    xx = np.reshape(x, (n, k))

    # precompute the multiplication of B with the parts of vector x
    Z = np.dot(B, xx.T)
    for j in range(n):
        yy += np.outer(A[:, j], Z[:, j])
    return np.ravel(yy)


def kron_d(A, B, x):
    """ Computes kron(A, B) * x (smart version without loop). """
    m, n = A.shape
    l, k = B.shape
    assert x.size == n * k, 'size mismatch'

    xx = np.reshape(x, (n, k))
    Z = np.dot(xx, B.T)
    yy = np.dot(A, Z)
    return np.ravel(yy)


def main():
    """ Main function, performs a scaling analysis and plots the results. """
    import timeit
    from matplotlib import pyplot as plt

    nruns = 3
    res = []
    nmax_a = 2**7
    for n in 2**np.mgrid[2:10]:
        print('n = {}:'.format(n))
        A = np.random.normal(size=(n, n))
        B = np.random.normal(size=(n, n))
        x = np.random.normal(size=n**2)

        # compute and measure
        if n < nmax_a:
            ya = kron_a(A, B, x)
            ta = min(timeit.repeat(lambda: kron_a(A, B, x),
                                   repeat=nruns, number=1))
        else:
            ta = None
        yb = kron_b(A, B, x)
        tb = min(timeit.repeat(lambda: kron_b(A, B, x),
                               repeat=nruns, number=1))
        yc = kron_c(A, B, x)
        tc = min(timeit.repeat(lambda: kron_c(A, B, x),
                               repeat=nruns, number=1))
        yd = kron_d(A, B, x)
        td = min(timeit.repeat(lambda: kron_d(A, B, x),
                               repeat=nruns, number=1))
        res.append((n, ta, tb, tc, td))

        # print errors
        print(' Errors:')
        if n < nmax_a:
            print('  a vs. b: {}'.format(np.linalg.norm(ya - yb)))
        print('  b vs. c: {}'.format(np.linalg.norm(yb - yc)))
        print('  b vs. d: {}'.format(np.linalg.norm(yb - yd)))

        # print timings
        print(' Timings [s]:')
        if n < nmax_a:
            print('  a: {}'.format(ta))
        print('  b: {}'.format(tb))
        print('  c: {}'.format(tc))
        print('  d: {}'.format(td))

    # plot results
    ns, tas, tbs, tcs, tds = np.transpose(res)
    c2 = np.sum(tbs) / np.sum(ns**2)
    c3 = np.sum(tbs) / np.sum(ns**3)
    mask = ns < nmax_a
    c4 = np.sum(tas[mask]) / np.sum(ns[mask]**4)
    plt.figure()
    plt.loglog(ns[mask], tas[mask], 'o-', label='slow evaluation')
    plt.loglog(ns, tbs, 'o-', label='efficient evaluation')
    plt.loglog(ns, tcs, 'o-', label='efficient ev. + reshape')
    plt.loglog(ns, tds, 'o-', label='efficient ev. + no loop')
    plt.loglog(ns, c2 * ns**2, 'k--', label=r'$\mathcal{O}(n^2)$')
    plt.loglog(ns, c3 * ns**3, 'k-.', label=r'$\mathcal{O}(n^3)$')
    plt.loglog(ns, c4 * ns**4, 'k:', label=r'$\mathcal{O}(n^4)$')
    plt.legend(loc='upper left')
    plt.show()


if __name__ == '__main__':
    main()
