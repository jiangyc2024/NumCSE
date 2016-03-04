import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt
import timeit

# test script: assessing the gain from using triangular solver instead of
# the default one

nruns = 3
res = []
for n in 2**np.mgrid[2:13]:
    # create test matrix
    A = np.triu(np.diag(np.mgrid[:n]) + np.ones((n, n)))
    # slight perturbation below the diagonal
    A += np.flipud(np.diag(np.finfo(float).eps * np.random.uniform(size=n)))
    b = np.random.uniform(size=n)

    ta = min(timeit.repeat(lambda: linalg.solve(A, b),
                           repeat=nruns, number=1))
    tb = min(timeit.repeat(lambda: linalg.solve_triangular(A, b),
                           repeat=nruns, number=1))
    res.append((n, ta, tb))

ns, tas, tbs = np.transpose(res)
plt.figure()
plt.loglog(ns, tas, 'o-', label='solve')
plt.loglog(ns, tbs, 'o-', label='solve_triangular')
plt.xlabel('matrix size n')
plt.ylabel('runtime for direct solve [s]')
plt.legend(loc='best')
plt.show()
