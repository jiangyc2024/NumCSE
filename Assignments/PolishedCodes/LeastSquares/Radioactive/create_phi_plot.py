import numpy as np
import matplotlib.pyplot as plt

lmda_1 = 1.
lmda_2 = 0.5
A_0 = 1.
B_0 = 0.

x = np.linspace(0, 5, 100)
phi_A = np.exp(-lmda_1 * x) * A_0
phi_B = np.exp(-lmda_2 * x) * B_0 + lmda_1 / (lmda_2 - lmda_1) * (np.exp(-lmda_1 * x) - np.exp(-lmda_2 * x)) * A_0

plt.figure()
plt.plot(x, phi_A, label='$\Phi_A$')
plt.plot(x, phi_B, label='$\Phi_B$')
plt.xlabel('time t')
plt.ylabel('amount of substances')
plt.title('Radioactive decay chain')
plt.legend()
plt.savefig('PhiAPhiB.eps')