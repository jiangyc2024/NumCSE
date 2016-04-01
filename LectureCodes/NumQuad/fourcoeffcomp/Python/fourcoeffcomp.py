import numpy as np


def fourcoeffcomp(c, m, ovsmpl=2):
    """ Compute the Fourier coefficients of the function c """
    ovsmpl = np.ceil(ovsmpl)
    # Number of quadrature points
    n = (2 * m + 1) * ovsmpl
    # Inverse DFT
    y = np.fft.ifft(c(np.linspace(0, 1, n, endpoint=False)))
    # Undo oversampling and wrapping of Fourier coefficient array
    return np.hstack([y[n - m:n], y[:m + 1]])


def main():
    y = fourcoeffcomp(lambda x: np.exp(1j * x * np.pi), 4)
    print(y)


if __name__ == '__main__':
    main()
