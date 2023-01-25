# -*- coding: utf-8 -*-
"""

Function that returns the kernels for convolution used by pascal bouda
comming from http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/

"""

from __future__ import print_function, division, absolute_import, unicode_literals

import numpy as np
import scipy.special


def _sym_coeff(coeff, sgn=1):
    return np.hstack((coeff[1:][::-1], sgn * coeff))


def d0coeff(m):
    """
    Returns the kernel used for smoothing.

    :param m: half_size
    :type m: integer

    :return: the kernel
    :rtype: numpy array
    """
    if m == 0:
        m = 1

    k = np.arange(m + 1).astype(np.float)

    s = (3 * m - 1 - 2 * (k ** 2)) / (2 * m - 1)
    s *= scipy.special.comb(2 * m, m + k) / (2 ** (2 * m))

    return _sym_coeff(s)


def d1coeff(m):
    """
    Returns the kernel used for first order differentiation.

    :param m: half_size
    :type m: integer

    :return: the kernel
    :rtype: numpy array
    """

    if m == 0:
        m = 1

    k = np.arange(1, m + 1).astype(np.float)

    s = np.zeros(k.size + 1)

    s[1:] = 1.0 / 2 ** (2 * m - 1)
    s[1:] *= scipy.special.comb(2 * m - 2, m - k) - scipy.special.comb(
        2 * m - 2, m - k - 2
    )

    return _sym_coeff(s, sgn=-1)


def d2coeff(m):
    """
    Returns the kernel used for second order differentiation.

    :param m: half_size
    :type m: integer

    :return: the kernel
    :rtype: numpy array
    """

    if m == 0:
        m = 1

    n = 2 * m + 1
    s = np.zeros(m + 2)

    s[m] = 1
    for k in np.arange(1, m + 1):
        s[m - k] = (
            (2 * n - 10) * s[m - k + 1] - (n + 2 * (m - k) + 3) * s[m - k + 2]
        ) / np.float(n - 2 * (m - k) - 1)

    s = s[: m + 1] / 2 ** (n - 3)

    return _sym_coeff(s)
