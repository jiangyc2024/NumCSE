import numpy as np


def adaptquad(f, m, rtol, atol):
    """ Adaptive numerical quadrature """
    h = np.diff(m)  # Distances of quadrature nodes
    mp = 0.5 * (m[:-1] + m[1:])  # Positions of midpoints
    fx = f(m)
    fm = f(mp)
    trp_loc = h * (fx[:-1] + 2 * fm + fx[1:]) / 4  # trapezoidal rule
    simp_loc = h * (fx[:-1] + 4 * fm + fx[1:]) / 6  # Simpson rule
    est_loc = np.abs(simp_loc - trp_loc)  # Local error estimation
    err_tot = np.sum(est_loc)  # Estimate for quadrature error
    s = np.sum(simp_loc)  # Simpson approximation of integral value
    # Termination check
    if err_tot > rtol * np.abs(s) and err_tot > atol:
        refcells = np.flatnonzero(est_loc > 0.9 * np.sum(est_loc) / len(h))
        return adaptquad(f, np.sort(np.hstack([m, mp[refcells]])), rtol, atol)
    return s


def main():
    def f(x):
        return np.exp(-x**2)

    m = np.array([-100, 0.1, 0.5, 100])
    print('Sqrt(pi) - Integral of exp(-x^2) in domain [-100, 100] = {}'.format(
          np.sqrt(np.pi) - adaptquad(f, m, 1e-10, 1e-12)))


if __name__ == '__main__':
    main()
