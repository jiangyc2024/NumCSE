# python script for testing the computation of the zeros of a parabola
import numpy as np


def zerosquadpol(alpha, beta):
    """ Computes the zeros of a quadratic polynomial
    p(x) = x^2 + alpha * x + beta by means of the familiar discriminant formula
    x{1,2} = 0.5 * (-alpha +- sqrt(alpha^2 - 4 * beta)). However this
    implementation is vulnerable to round-off! """
    D = alpha**2 - 4 * beta
    if D < 0:
        return None
    wD = np.sqrt(D)
    return [0.5 * (-alpha - wD), 0.5 * (-alpha + wD)]


def zerosquadpolstab(alpha, beta):
    """ Computes the zeros of a quadratic polynomial
    p(x) = x^2 + alpha * x + beta by means of the familiar discriminant formula
    x{1,2} = 0.5 * (-alpha +- sqrt(alpha^2 - 4 * beta)). This is a stable
    implementation based on Vieta's theorem. """
    D = alpha**2 - 4 * beta
    if D < 0:
        return None
    wD = np.sqrt(D)
    # use discriminant formula only for zero far away from 0 in order to avoid
    # cancellation. For the other zeros use Vieta's formula.
    if alpha >= 0:
        t = 0.5 * (-alpha - wD)
        return [t, beta / t]
    else:
        t = 0.5 * (-alpha + wD)
        return [beta / t, beta]


def main():
    from matplotlib import pyplot as plt

    res = []  # array for storing the results
    gammavec = np.mgrid[2:1000:10]  # test cases

    for gamma in gammavec:
        # polynomial p(x) = (x - gamma) * (x - 1 / gamma)
        # compute coefficients
        alpha = -(gamma + 1.0 / gamma)
        beta = 1.0
        z1 = zerosquadpol(alpha, beta)  # unstable ways to compute zeros
        z2 = zerosquadpolstab(alpha, beta)  # stable implementation
        # compute relative errors of the left zero
        ztrue = 1.0 / gamma
        z2true = gamma
        res.append((gamma, np.abs((z1[0] - ztrue) / ztrue),
                    np.abs((z2[0] - ztrue) / ztrue),
                    np.abs((z1[1] - z2true) / z2true)))

    # graphical output of relative error of toors computed by unstable
    # implementation
    gammas, z1err_smallroot, z2err_smallroot, z1err_largeroot = zip(*res)
    plt.figure()
    plt.plot(gammas, z1err_smallroot, 'o-', label='small root')
    plt.plot(gammas, z1err_largeroot, 'o-', label='large root')
    plt.xlabel(r'$\gamma$')
    plt.ylabel(r'relative errors in $\xi_1$, $\xi_2$')
    plt.title('Roots of parabola computed in an unstable manner')
    plt.legend(loc='best')
    # graphical output of relative errors (comparison)
    plt.figure()
    plt.plot(gammas, z1err_smallroot, 'o-', label='unstable')
    plt.plot(gammas, z2err_smallroot, 'o-', label='stable')
    plt.xlabel(r'$\gamma$')
    plt.ylabel(r'relative error in $\xi_1$')
    plt.title('Roundoff in the computation of zeros of a parabola')
    plt.legend(loc='best')
    plt.show()


if __name__ == '__main__':
    main()
