import numpy as np
import timeit
from matplotlib import pyplot as plt


def col_wise(A):
    for j in range(A.shape[1] - 1):
        A[:, j + 1] -= A[:, j]


def row_wise(A):
    for i in range(A.shape[0] - 1):
        A[i + 1, :] -= A[i, :]


# Timing for row/column-wise operations on matrix, we conduct several runs in
# order to reduce risk of skewed measurements due to OS activity during run.

nruns = 3
res = []
for n in 2**np.mgrid[4:14]:
    A = np.random.normal(size=(n, n))

    t1 = min(timeit.repeat(lambda: col_wise(A), repeat=nruns, number=1))
    t2 = min(timeit.repeat(lambda: row_wise(A), repeat=nruns, number=1))

    res.append((n, t1, t2))

# plot runtime versus matrix sizes
ns, t1s, t2s = np.transpose(res)

plt.figure()
plt.plot(ns, t1s, 'o-', label='A[:, j + 1] -= A[:, j]')
plt.plot(ns, t2s, 'o-', label='A[i + 1, :] -= A[i, :]')
plt.xlabel(r'$n$')
plt.ylabel(r'runtime [s]')
plt.legend(loc='upper left')

plt.figure()
plt.loglog(ns, t1s, 'o-', label='A[:, j + 1] -= A[:, j]')
plt.loglog(ns, t2s, 'o-', label='A[i + 1, :] -= A[i, :]')
plt.xlabel(r'$n$')
plt.ylabel(r'runtime [s]')
plt.legend(loc='upper left')

plt.show()
