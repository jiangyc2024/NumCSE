import numpy as np

y = []
x = 0.4
for _ in range(15):
    x += (np.cos(x) + 1) / np.sin(x)
    y.append(x)
err = np.array(y) - x
rate = err[1:] / err[:-1]
print(rate)
