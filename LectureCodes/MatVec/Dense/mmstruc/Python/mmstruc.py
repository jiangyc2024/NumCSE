import numpy as np
from matplotlib import pyplot as plt


n = 100
A = np.diag(np.mgrid[:n])
A[:, -1] = np.mgrid[:n]
A[-1, :] = np.mgrid[:n]
B = A[::-1, :]
C = np.dot(A, A)
D = np.dot(A, B)

plt.figure()
plt.spy(A)
plt.savefig('../PYTHON_PICTURES/Aspy.eps')
plt.figure()
plt.spy(B)
plt.savefig('../PYTHON_PICTURES/Bspy.eps')
plt.figure()
plt.spy(C)
plt.savefig('../PYTHON_PICTURES/Cspy.eps')
plt.figure()
plt.spy(D)
plt.savefig('../PYTHON_PICTURES/Dspy.eps')
plt.show()
