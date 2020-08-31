import pandas as pd
import numpy as np
import math
import matplotlib
matplotlib.use('ps')
import matplotlib.pyplot as plt
from matplotlib import rc

x = pd.read_fwf( "timings.txt", widths = [10, 10, 10, 10, 10]).to_numpy( )
print( x )
plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble" : r'\usepackage{color}',
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

plt.figure( )
plt.title( "access to sparse matrices in " + r"\textcolor{blue}{CCS}"+ " and " + r"\textcolor{red}{CRS}" + "\nperformance vs. matrix size" )
n = np.power( 2, x[:,0])
plt.loglog( n, x[:,1], "b--", label = r"\textcolor{blue}{CCS}"+ " row access" )
plt.loglog( n, x[:,2], "bd-", label = r"\textcolor{blue}{CCS}"+ " column access" )
plt.loglog( n, x[:,3], "r--", label = r"\textcolor{red}{CRS}" + " row access" )
plt.loglog( n, x[:,4], "rd-", label = r"\textcolor{red}{CRS}" + " column access" )
plt.loglog( n, 1e-1 * np.power(n, 1), "k--", label = "$O(n)$" )
plt.ylabel( "time [µs]" )
plt.xlabel( "matrix size (n)" )
plt.legend( )
plt.savefig( "./formataccess.eps" )
