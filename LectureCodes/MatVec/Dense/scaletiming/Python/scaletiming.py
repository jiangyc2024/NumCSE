import numpy as np
from matplotlib import pyplot as plt
import timeit

# script for timing a smart and foolish way to carrry out
# multiplication with a scaling matrix

nruns = 3
res = []
for n in 2**np.mgrid[2:15]:
    d = np.random.uniform(size=n)
    x = np.random.uniform(size=n)

    tbad = min(timeit.repeat(lambda: np.dot(np.diag(d), x),
                             repeat=nruns, number=1))
    tgood = min(timeit.repeat(lambda: d * x, repeat=nruns,
                              number=1))

    res.append((n, tbad, tgood))

ns, tbads, tgoods = np.transpose(res)
plt.figure()
plt.loglog(ns, tbads, 'o-', label='using np.diag')
plt.loglog(ns, tgoods, 'o-', label='using *')
plt.legend(loc='best')
plt.title('Timing for different ways to do scaling')
plt.show()
