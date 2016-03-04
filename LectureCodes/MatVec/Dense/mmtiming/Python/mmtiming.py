# script for timing different implementations of matrix multiplications
import numpy as np
from matplotlib import pyplot as plt
import timeit


def mm_loop_based(A, B, C):
    m, n = A.shape
    _, p = B.shape
    for i in range(m):
        for j in range(p):
            for k in range(n):
                C[i, j] += A[i, k] * B[k, j]
    return C


def mm_blas1(A, B, C):
    m, n = A.shape
    _, p = B.shape
    for i in range(m):
        for j in range(p):
            C[i, j] = np.dot(A[i, :], B[:, j])
    return C


def mm_blas2(A, B, C):
    m, n = A.shape
    _, p = B.shape
    for i in range(m):
        C[i, :] = np.dot(A[i, :], B)
    return C


def mm_blas3(A, B, C):
    C = np.dot(A, B)
    return C


def main():
    nruns = 3
    res = []
    for n in 2**np.mgrid[2:11]:
        print('matrix size n = {}'.format(n))
        A = np.random.uniform(size=(n, n))
        B = np.random.uniform(size=(n, n))
        C = np.random.uniform(size=(n, n))

        tloop = min(timeit.repeat(lambda: mm_loop_based(A, B, C),
                                  repeat=nruns, number=1))
        tblas1 = min(timeit.repeat(lambda: mm_blas1(A, B, C),
                                   repeat=nruns, number=1))
        tblas2 = min(timeit.repeat(lambda: mm_blas2(A, B, C),
                                   repeat=nruns, number=1))
        tblas3 = min(timeit.repeat(lambda: mm_blas3(A, B, C),
                                   repeat=nruns, number=1))
        res.append((n, tloop, tblas1, tblas2, tblas3))

    ns, tloops, tblas1s, tblas2s, tblas3s = np.transpose(res)
    plt.figure()
    plt.loglog(ns, tloops, 'o-', label='loop implementation')
    plt.loglog(ns, tblas1s, 'o-', label='dot product implementation')
    plt.loglog(ns, tblas2s, 'o-', label='matrix-vector implementation')
    plt.loglog(ns, tblas3s, 'o-', label='BLAS gemm (np.dot)')
    plt.legend(loc='upper left')
    plt.show()


if __name__ == '__main__':
    main()
