import sys
import math
import numpy as np
import numpy.linalg as la 
import matplotlib.pyplot as plt
import scipy.interpolate as intp
from scipy.spatial import ConvexHull
import bezier 

def plotheart(filename ="heartcurve.eps"):
    t = np.linspace(0,2*math.pi,200)
    x = np.array([16*math.sin(x)**3 for x in t])
    y = np.array([13*math.cos(x)-5*math.cos(2*x) - 2*math.cos(3*x) - math.cos(4*x) for x in t])
    fig, ax = plt.subplots()
    plt.plot(x,y,'r-',label="curve")
    plt.text(0.0,0.0,"I love numerics",color='r',horizontalalignment='center',fontsize=14)
    plt.xlabel("$x_1$")
    plt.ylabel("$x_2$")
    plt.legend(loc='best')
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)

def plotparabola(filename = "parabolacurve.eps"):
    t = np.linspace(-1.0,1.0,200)
    fig, ax = plt.subplots()
    plt.plot(t,t**2,'r-',label="parabola")
    plt.xlabel("$x_1$")
    plt.ylabel("$x_2$")
    plt.legend(loc='best')
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)
    
def knotsfrompoints(points):
    """
    Computes knots from control points using arclength of control polygon 
    For n points also n knots are returned (not normalized to [0,1]!
    """
    if points.shape[0] != 2:
        raise RuntimeError("Points must have two coordinates")
    n = points.shape[1]; # number of points
    if n < 2:
        raise RuntimeError("At least two points required!")
    # Obtain interpolation nodes from distances of points 
    return np.cumsum(np.hstack((np.array([0.0]),la.norm(np.diff(points),axis=0))))

# Polynomial curve
def polycurve(points, N = 100):
    """
    Interpolating polynomial curve through given points in the plane.
    The 'points' argument must pass a 2xn matrix whose columns contain 
    the point coordinates.
    The 'N' argument tells the number of evaluation points 
    The function returns points on the curve correspondig to N evenly spaced 
    points in the parameter interval 
    """
    # Obtain interpolation nodes from distances of points 
    knots = knotsfrompoints(points)
    # Abscissas for evaluation
    t = np.linspace(knots[0],knots[-1],N)
    # Barycentric polynomial interpolation
    return intp.BarycentricInterpolator(knots,points,axis=1)(t)

# Curve generated by cubic spline interpolation
def cubicsplinecurve(points, N = 100):
    """
    Use library routine to compute cubic spline interpolant running through control points
    """
    # Obtain interpolation nodes from distances of points 
    knots = knotsfrompoints(points)
    # Abscissas for evaluation
    t = np.linspace(knots[0],knots[-1],N)
    return intp.CubicSpline(knots,points,axis=1,bc_type='natural')(t)

# Curve generated by PCHIP interpolation
def pchipcurve(points, N = 100):
    """
    Use library routine to compute monotonic Hermite interpolant running through control points
    """
    # Obtain interpolation nodes from distances of points 
    knots = knotsfrompoints(points)
    # Abscissas for evaluation
    t = np.linspace(knots[0],knots[-1],N)
    return intp.PchipInterpolator(knots,points,axis=1)(t)

# Computation of Bernstein polynomials (on reference interval)
def BernsteinPolynomials(n,t):
    """Computes a matrix whose i-th row contains the point values of the Bernstein 
    polynomials of degree n with index i=0,...,N evaluated at the values passed in t.
    """
    if len(t.shape) > 1:
        raise RuntimeError("t must be a 1D array!")
    N = len(t) # Number of evaluation points
    if n == 0:
        return np.ones((1,N))
    else:
        t = t.flatten()
        X = np.vstack((np.ones((1,N)),np.zeros((n,N))))
        # Evaluation of recursion formula
        for d in range(1,n+1):
            for i in range(d,0,-1):
                X[i,:] = t*X[i-1,:] + (1-t)*X[i,:]
            X[0,:] = (1-t)*X[0,:]
        return X


# Computations of B-splines
def BSplines(tau,d,t):
    """
    Computes the values of all B-splines of a given degree d and for the knot set tau 
    for the evaluation points passed in t. Returns the values in an (n+1)xN-matrix, 
    where N is the number of elements in t
    """
    if len(t.shape) > 1:
        raise RuntimeError("t must be a 1D array!")
    N = len(t) # Number of evaluation points
    if len(tau.shape) > 1:
        raise RuntimeError("tau must be a 1D array!")
    m = len(tau)-1 # knots numbered 0, ... , m 
    if m < d:
        raise RuntimeError("degree %d >= number of nodes %d" % (d,m+1))
    if (t[0] < tau[0]) or (t[-1] > tau[-1]):
        raise RuntimeError("Evaluation points out of knot range!")
    # First initialize the characteristic functions
    # print(("tau = ",tau))
    # print(("t = ",t))
    X = np.zeros((m,N))
    l = 0;
    k = 0;
    while (k<N) and (l <= m):
        if t[k] < tau[l]:
            X[l-1,k] = 1.0
            k += 1
        else:
            l += 1
    # Recursion
    for k in range(1,d+1):
        for i in range(0,m-k):
            if tau[i+k] > tau[i]:
                X[i,:] = (t-tau[i])/(tau[i+k]-tau[i])*X[i,:]
            if tau[i+k+1] > tau[i+1]:
                X[i,:] += (tau[i+k+1]-t)/(tau[i+k+1]-tau[i+1])*X[i+1,:]
    return X[0:m-d,:]
  
            
# Inefficient implementation of evaluation of Bezier curve
def Beziercurve(points, N = 100):
    """Computes points on a Bezier curve given control points 
    and for equidistant evaluation points
    """
    # print(("points = ",points))
    if len(points.shape) != 2:
        raise RuntimeError("Points must be a 2xn array of floats")
    if points.shape[0] != 2:
        raise RuntimeError("Points must have two coordinates")
    n = points.shape[1]; # number of points
    # Evaluations points in [0,1]
    t = np.linspace(0,1,N)
    # Compute all Bernstein polynomials of degree n-1
    X = BernsteinPolynomials(n-1,t)
    # print(("X.shape = ",X.shape))
    return points.dot(X)

# Inefficient implementation of B-spline curve 
def BSplinecurve(points, d, N = 100):
    """
    Given a 2xn matrix of point coordinates (in its columns) and the degree d, 
    compute knots from point distances and then the corresponding B-spline curve
    with those control points. 
    """
    if points.shape[0] != 2:
        raise RuntimeError("Points must have two coordinates")
    m = points.shape[1]; # number of points
    if m < 2:
        raise RuntimeError("At least two points required!")
    # For m control points we need n = m - d + 1 knots
    # Choose them equispaced in [0,1]
    n = m-d+1
    knots = np.linspace(0.0,1.0,n)
    # Pad by d more knots at the ends of the interval
    tau = np.hstack((np.ones(d)*knots[0],knots,np.ones(d)*knots[-1]));
    # Abscissas for evaluation
    epsval = 100*np.finfo(float).eps
    t = np.linspace(knots[0]+epsval,knots[-1]-epsval,N)
    # Matrix X of B-spline values, size (n+2d-1-d = m)xN
    X = BSplines(tau,d,t)
    print(("Last column of X",X[:,-1]))
    return points.dot(X)

def plotControlPoints(points, filename = "controlpoints.eps"):
    """
    Simple plotting of control points 
    """
    if points.shape[0] != 2:
        raise RuntimeError("Points must have two coordinates")
    m = points.shape[1]; # number of points
    if m < 2:
        raise RuntimeError("At least two points required!")
    fig, ax = plt.subplots()
    ax.set_ylim([-0.5,6.5])
    ax.set_xlim([-4,16])
    ax.set_aspect('equal','box')
    plt.plot(points[0,:],points[1,:],'ro',label="control points")
    plt.plot(points[0,:],points[1,:],'k-',linewidth=0.5,label="control polygon")
    for i in range(0,m):
        plt.text(points[0,i],points[1,i],'%d' % i,color = 'b')
    plt.grid()
    plt.xlabel("$x_1$")
    plt.ylabel("$x_2$")
    plt.legend(loc='best')
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)
    

def polycurve_plot(points,N=100):
    """Plotting of polynomial interpolating curve"""
    y = polycurve(points,N)
    fig, ax = plt.subplots()
    ax.set_aspect('equal','box')
    plt.plot(points[0,:],points[1,:],'r-+',linewidth=0.5,label="polygon")
    plt.plot(y[0,:],y[1,:],'b-',linewidth=1,label="polynomial curve")
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.legend(bbox_to_anchor=(0.1,1.1))
    plt.show(block = False)
    
def interpolators_plot(points,filename = "interpolatingcurves.eps", N=100):
    """Plotting of various interpolating curves"""
    y_poly = polycurve(points,N)
    y_spline = cubicsplinecurve(points,N)
    y_pchip = pchipcurve(points,N)
    fig, ax = plt.subplots()
    ax.set_ylim([-0.5,6.5])
    ax.set_xlim([-4,16])
    ax.set_aspect('equal','box')
    plt.plot(points[0,:],points[1,:],'r-+',linewidth=0.5,label="control polygon")
    plt.plot(y_poly[0,:],y_poly[1,:],'b--',linewidth=1,label="polynomial interpolant")
    plt.plot(y_spline[0,:],y_spline[1,:],'m--',linewidth=1,label="$C^2$ cubic spline interpolant")
    plt.plot(y_pchip[0,:],y_pchip[1,:],'g--',linewidth=1,label="$C^1$ Hermite interpolant")
    plt.xlabel("$x_1$")
    plt.ylabel("$x_2$")
    plt.legend(loc='best')
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)

def plotBernstein(n, filename = "bernstein.eps", N = 100):
    t = np.linspace(0,1,N)
    X = BernsteinPolynomials(n,t)
    fig, ax = plt.subplots()
    for i in range(0,n+1):
        plt.plot(t,X[i,:],label="${B}^{%d}_{%d}$" % (n,i))
    plt.xlabel('$t$')
    plt.ylabel('${B}_{i}^{d}(t)$')
    plt.title("Bernstein polynomials, degree = %d" % n)
    plt.legend(loc='best')
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)

def plotBSplines(tau,d,filename="bspline.eps",N=100):
    """
    Plotting non-zero B-splines for a knot set tau, which is assumed to be sorted
    """
    epsval = 100*np.finfo(float).eps
    t = np.linspace(tau[0]+epsval,tau[-1]-epsval,N)
    m = len(tau)-1
    X = BSplines(tau,d,t)
    print(("Column sum = ",X.sum(axis=0)))
    fig, ax = plt.subplots()
    for i in range(0,m-d):
        plt.plot(t,X[i,:],label="$N^{%d}_{%d}$" % (d,i))
    plt.plot(tau,np.zeros(len(tau)),'r+',label="knots")
    plt.xlabel('$t$')
    plt.ylabel('$N_{i}^{d}(t)$')
    plt.title("%d B-splines, degree = %d, %d knots" % (m-d, d, len(tau)) )
    plt.legend(loc='best')
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)  # ("d=%d%s" % (d,filename, bbox_inches='tight',pad_inches=0.0))
    plt.show(block=False)
    
def plotBezier(points, filename = "beziercurve.eps", N = 100):
    """
    Plotting a simple Bezier curve
    """
    y = Beziercurve(points,N)
    fig, ax = plt.subplots()
    ax.set_ylim([-0.5,6.5])
    ax.set_xlim([-3.5,16])
    ax.set_aspect('equal','box')
    plt.plot(points[0,:],points[1,:],'r-+',linewidth=0.5,label="polygon")
    plt.plot(y[0,:],y[1,:],'b-',linewidth=1,label="Bezier curve")
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.title("Bezier curve, %d control points" % points.shape[1])
    plt.legend(loc='best')
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)    

def plotSplineCurve(points, d, filename = "splinecurve.eps", N = 100):
    """
    Plotting a spline curve determined by control points 
    """
    y = BSplinecurve(points,d,N)
    fig, ax = plt.subplots()
    ax.set_aspect('equal','box')
    ax.set_ylim([-0.5,6.5])
    ax.set_xlim([-3.5,16])
    plt.plot(points[0,:],points[1,:],'r-+',linewidth=0.5,label="control polygon")
    plt.plot(y[0,:],y[1,:],'b-',linewidth=1,label="spline curve")
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.title("Spline curve, degree %d, %d control points" % (d, points.shape[1]))
    plt.legend(loc='best')
    plt.savefig(filename,bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)    
    
def plotBezierLib(points, filename = "beziercurvelib.eps", N = 100):
    """
    Plotting a Bezier curve using the library function 
    """
    if len(points.shape) != 2:
        raise RuntimeError("Points must be a 2xn array of floats")
    if points.shape[0] != 2:
        raise RuntimeError("Points must have two coordinates")
    n = points.shape[1]; # number of points
    c = bezier.Curve(points,degree=n-1)
    y = c.evaluate_multi(np.linspace(0,1,N))
    fig, ax = plt.subplots()
    ax.set_aspect('equal','box')
    plt.plot(points[0,:],points[1,:],'r-+',linewidth=0.5,label="polygon")
    plt.plot(y[0,:],y[1,:],'b-',linewidth=1,label="Bezier curve")
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.title("Bezier curve, %d control point, library function" % n)
    plt.legend(bbox_to_anchor=(0.4,-0.2))
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)    

def Bezierconstruct(filename = "bezierconstruct.eps", N = 100):
    """
    Demonstration of convex combinations
    """
    nodes = np.array([[0,1],[0.5,0],[1,1],[2,1]]).T
    fig, ax = plt.subplots()
    plt.plot(nodes[0,:],nodes[1,:],'b-+',linewidth=0.5,label="polygon, $d=1$")
    y2_left = Beziercurve(nodes[:,0:3],N)
    y2_right = Beziercurve(nodes[:,1:4],N)
    plt.plot(y2_left[0,:],y2_left[1,:],'m-',linewidth=0.5,label="left part $\\mathbf{b}_0^2$, $d=2$")
    plt.plot(y2_right[0,:],y2_right[1,:],'c-',linewidth=0.5,label="right part $\\mathbf{b}_1^2$, $d=2$")
    y3 = Beziercurve(nodes,N)
    plt.plot(y3[0,:],y3[1,:],'r-',linewidth=0.5,label="Bezier curve $\\mathbf{b}^3$, $d=3$")
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.legend(loc='best')
    plt.title("Construction of cubic Bezier curve, %d control points" % 4)
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)

def plotConvexHull(filename = "convexhull.eps"):
    points = np.random.rand(30, 2)   # 30 random points in 2-D
    hull = ConvexHull(points)
    fig, ax = plt.subplots()
    plt.plot(points[:,0], points[:,1], 'bo')
    extvt = np.hstack((hull.vertices,[hull.vertices[0]])) 
    plt.fill(points[extvt,0], points[extvt,1], 'xkcd:light yellow')
    plt.plot(points[extvt,0], points[extvt,1], 'ro-', lw=2)
    plt.savefig(filename, bbox_inches='tight',pad_inches=0.0)
    plt.show(block = False)
    
if __name__ == "__main__":
    plt.rcParams.update({
        "text.usetex": True})
    # The "car shape"
    nodes = np.array([[-1,0],[0,2],[4,2],[5.5,4],[9.5,4],[11,2],[15,2],[15.5,3],[15.5,0]]).T
    print("***** Validation *****")
    n = nodes.shape[1]; # number of points
    c = bezier.Curve(nodes, degree=n-1)
    y = c.evaluate_multi(np.linspace(0.0,1.0,10))
    print(("y = ",y))
    print("***** Plotting nodes *****")
    plotControlPoints(nodes)
    print("**** Plotting convex hull ****")
    plotConvexHull()
    print("****** Plotting heart curve ******")
    plotheart()
    print("****** Plotting parabola ********")
    plotparabola()
    print("****** Plotting interpolating curves *****")
    interpolators_plot(nodes)
    print("Visualizing construction of Bezier curves")
    Bezierconstruct()
    print("***** Bezier curve for car shape *****")
    plotBezier(nodes)
    print("****** Plotting Bernstein polynomials *****")
    plotBernstein(8)
    print("*** Plotting B-splines ****")
    tau = np.hstack((np.array([0.0, 0.0]),np.linspace(0.0,1.0,8),np.array([1.0, 1.0])))
    plotBSplines(tau,2, filename = "bspline2.eps")
    tau = np.hstack((np.array([0.0, 0.0, 0.0]),np.linspace(0.0,1.0,7),np.array([1.0, 1.0, 1.0])))
    plotBSplines(tau,3, filename = "bspline3.eps")
    print("**** Plotting car spline curve ****")
    plotSplineCurve(nodes,3,filename="cubicsplinecar.eps")
    plotSplineCurve(nodes,2,filename="quadraticsplinecar.eps")
    print("""
    The following files should have been created:
     bernstein.eps          
     bezierconstruct.eps    
     beziercurve.eps        
     bspline2.eps           
     bspline3.eps           
     cubicsplinecar.eps     
     heartcurve.eps         
     interpolatingcurves.eps
     parabolacurve.eps      
     quadraticsplinecar.eps 
    """)
    
