import numpy as np
from scipy import linalg


def F(x):
    return np.array([x[0]**2 - x[1]**4, x[0] - x[1]**3])


def DF(x):
    return np.array([[2 * x[0], -4 * x[1]**3], [1, - 3 * x[1]**2]])


x = np.array([0.7, 0.7])
x_ast = np.array([1, 1])
tol = 1e-10

errors = [linalg.norm(x - x_ast)]
while True:
    s = linalg.solve(DF(x), F(x))
    x -= s
    errors.append(linalg.norm(x - x_ast))
    if linalg.norm(s) <= tol * linalg.norm(x):
        break
ld = np.diff(np.log(errors))
rates = ld[1:] / ld[:-1]
print('Errors:\n{}'.format(np.array(errors)))
print('Rates:\n{}'.format(rates))
