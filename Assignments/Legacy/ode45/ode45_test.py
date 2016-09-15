#! /usr/bin/env python2

from numpy import *
#from matplotlib.pyplot import *
from ode45 import ode45

def test1():
    # Parameter und Anfangsbedingungen
    T = 10.
    a1 = 3
    a2 = 2
    b1 = 0.1
    b2 = 0.1
    f = lambda t, y: array([(a1 - b1*y[1])*y[0], (b2*y[0] - a2)*y[1]])
    y0 = array([100., 5.])
    tspan = [0., T]

    # RK45
    t45, y45 = ode45(f, tspan, y0) 

    #print t45
    #print y45
    print len(t45)

def test2():
    # Parameter und Anfangsbedingungen
    T = 10000.
    f = lambda t, y: [1 / y]
    y0 = array([0.2])
    tspan = [0., T]

    # RK45
    t45, y45 = ode45(f, tspan, y0) 

    #print t45
    #print y45
    print len(t45)

if __name__ == '__main__':
    #test1()
    test2()
