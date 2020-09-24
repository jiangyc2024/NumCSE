import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

x = pd.read_fwf( "timings.txt", widths = [10, 10, 10, 10]).to_numpy( )
print( x )
plt.figure( )
plt.title( "initialization of a banded matrix\nperformance vs. matrix size" )
plt.loglog( x[:,0], x[:,1], "g-", label = "triplets" )
plt.loglog( x[:,0], x[:,2], "b-", label = "coeffRef (sufficient reservation)" )
plt.loglog( x[:,0], x[:,3], "r-", label = "coeffRef (insufficient reservation)" )
plt.loglog( x[:,0], 1e-2 * np.power(x[:,0], 2), "k--", label = "O(n²)" )
plt.ylabel( "time [µs]" )
plt.xlabel( "matrix size (n)" )
plt.legend( )
plt.savefig( "./sparsetiming.eps" )
