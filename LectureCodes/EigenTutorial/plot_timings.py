# Data visualization script for ncse_eigenopstiming
import numpy as np
import matplotlib.pyplot as plt

print("Reading data from file eigenopstimings.csv");
data = np.genfromtxt("eigenoptimings.csv", delimiter = ',')
n = data[:,0];
print("Timings for n = ",n);

fig = plt.figure()
plt.loglog(n, data[:,2], 'b+-', linewidth=0.5,label="matrix-vector product")
plt.loglog(n, data[:,3], 'r*-', linewidth=0.5,label="matrix-matrix product")
plt.loglog(n,1.0E-3*n**2,'k-',linewidth=0.5,label="O(n^2)");
plt.loglog(n,1.0E-3*n**3,'k--',linewidth=0.5,label="O(n^3)");
plt.title("Measured runtimes, code eigenopstiming.cpp")
plt.xlabel('problem size parameter n')
plt.ylabel('runtime (microseconds)')
plt.legend()
plt.savefig("eigenopstiming.eps")

