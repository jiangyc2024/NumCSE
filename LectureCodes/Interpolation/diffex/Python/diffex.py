# Pyhton script for generating visuals for the section on "Extrapolation to zero" of the
# Course "Numerical methods for CSE"
# Author: Ralf Hiptmair, SAM, D-MATH
# Date: October 2020

import sys
import math
import numpy as np
import numpy.linalg as la 
import matplotlib.pyplot as plt

def diffexvisual():
    """
    Pyhton function creating som visuals for section of extrapolation to zero 
    """
    # Data points sampled from arctan 
    h = np.array([1.0,1.0/2,1.0/3,1.0/4])
    y = np.arctan(2*h)
    # Compute and evaluate polynomial
    x = np.linspace(-0.25,1.0,200)
    # Create figure 
    fig, ax = plt.subplots()
    plt.plot(h,y,'r*',label="$(h_i,\psi(h_i))$")
    plt.plot(np.array([-0.25,1.0]),np.array([0.0,0.0]),'k-')
    plt.plot(np.array([0.0,0.0]),np.array([-0.75,1.0]),'k-')
    npts = 2
    plt.plot(x,np.polyval(np.polyfit(h[0:npts],y[0:npts],npts-1),x),'r-',linewidth=0.5,label="degree = %d" % (npts-1))
    npts = 3
    plt.plot(x,np.polyval(np.polyfit(h[0:npts],y[0:npts],npts-1),x),'m-',linewidth=0.5,label="degree = %d" % (npts-1))
    npts = 4
    plt.plot(x,np.polyval(np.polyfit(h[0:npts],y[0:npts],npts-1),x),'b-',linewidth=0.5,label="degree = %d" % (npts-1))
    plt.xlabel("t")
    plt.ylabel("p")
    plt.legend(loc="best")
    plt.savefig("extpatan.eps", bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)

if __name__ == "__main__":
    plt.rcParams.update({"text.usetex": True})
    print("Visualization of extrapolation method\n") 
    diffexvisual();
