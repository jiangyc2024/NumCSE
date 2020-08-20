import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.io import ascii

x = pd.read_fwf( "timings.txt", widths = [5, 10, 10, 10, 10]).to_numpy( )
print( x )
plt.figure( )
plt.title( "inversion of diagonal matrices\nperformance vs. matrix size" )
plt.loglog( x[:,0], x[:,1], "r-", label = "explicit" )
plt.loglog( x[:,0], x[:,2], "g-", label = "auto" )
plt.loglog( x[:,0], 1e-3 * np.power(x[:,0], 3), "k--", label = "O(n³)" )
plt.ylabel( "time [µs]" )
plt.xlabel( "matrix size (n)" )
plt.legend( )
plt.savefig( "./invertPlot.eps" )

plt.figure( )
plt.title( "basic arithmetic on square matrices\nperformance vs. matrix size" )
plt.plot( x[:,0], x[:,3], "r-", label = "sequential" )
plt.plot( x[:,0], x[:,4], "g-", label = "vectorized" )
plt.ylabel( "time [µs]" )
plt.xlabel( "matrix size (n)" )
plt.legend( )
plt.savefig( "./arithPlot.eps" )