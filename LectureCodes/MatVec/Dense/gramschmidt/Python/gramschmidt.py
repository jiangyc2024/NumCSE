import numpy as np
from scipy.linalg import hilbert


def gramschmidt(A):
    """ Gram-Schmidt orthogonalization of column-vectors. Matrix A passes
    vectors in its columns, orthonormal system is returned in columns of
    matrix Q. """
    _, k = A.shape

    # first basis vector
    Q = A[:, [0]] / np.linalg.norm(A[:, 0])
    for j in range(1, k):
        # orthogonal projection, loop-free implementation
        q = A[:, j] - np.dot(Q, np.dot(Q.T, A[:, j]))

        # check premature termination
        nq = np.linalg.norm(q)
        if nq < 1e-9 * np.linalg.norm(A[:, j]):
            break
        # add new basis vector as another column of Q
        Q = np.column_stack([Q, q / nq])
    return Q


def main():
    """ Main function, demonstrates roundoff on the result of the Gram-Schmidt
    orthogonalization. """
    # set print options to use lower precision
    printopt = np.get_printoptions()
    np.set_printoptions(formatter={'float': '{:8.2g}'.format}, linewidth=200)

    # create special matrix, the so-called Hilbert-matrix Aij = 1 / (i + j + 1)
    A = hilbert(10)
    Q = gramschmidt(A)

    # matrix according to theory should be unit matrix:
    I = np.dot(Q.T, Q)
    print('I = \n{}'.format(I))

    # numpy's internal orthogonaliztation by QR-decomposition
    Q1, R1 = np.linalg.qr(A)
    D = A - np.dot(Q1, R1)
    print('D = \n{}'.format(D))
    I1 = np.dot(Q1.T, Q1)
    print('I1 = \n{}'.format(I1))

    np.set_printoptions(**printopt)


if __name__ == '__main__':
    main()
