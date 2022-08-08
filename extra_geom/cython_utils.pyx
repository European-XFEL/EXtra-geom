cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def condat_shift(double[:] signal, double tau):
    """Shift a 1D array of pixels by a fractional amount. Used to transform DSSC data.

    This implements a Type II filter for N = 2, as described in
    Condat et al. 2008, *Reversible, Fast, and High-Quality Grid
    Conversions*, IEEE Transactions on Image Processing.
    There is a simple implementation for N = 2 (p 686 in the paper);
    generalising to other values of N would be more complex.

    The data in signal is shifted in place.
    """
    cdef:
        double a_numerator, a_plus_tau, a_minus_tau
        int j, n

    if tau < -0.5 or tau > 0.5:
        raise ValueError(f"Condat filter delay is out of range [-0.5, 0.5]: {tau}")

    a_numerator = -4 + tau**2 + (12 - 3 * tau**2) ** 0.5
    a_plus_tau = a_numerator / (3 * tau + 2 + tau**2)
    a_minus_tau = a_numerator / (-3 * tau + 2 + tau**2)

    n = len(signal)
    # TODO - I think the limits of the iteration here are slightly wrong -
    # I would put range(1, n) and range(n, 0, -1). But for now, I'm aiming to
    # preserve exactly the behaviour of the original code. -TK
    for j in range(1, n - 3):
        signal[j] += a_plus_tau * (signal[j - 1] - signal[j + 1])
    for j in range(n - 3, 1, -1):
        signal[j] += a_minus_tau * (signal[j + 1] - signal[j - 1])
