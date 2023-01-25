# -*- coding: utf-8 -*-
"""
Here are some functions to perform smoothing and differentiation mainly using convolution.
"""
# ==============================================================================
# standard libraries
# ==============================================================================
import numpy as np
import scipy.signal as sg
import scipy.special as sp

# from sympy.utilities.lambdify import lambdify
from scipy.integrate import simps, cumtrapz

# ==============================================================================
# libraries in candy
# ==============================================================================
from pycurves.courbe import *
from pycurves.courbe_tools import *
from pycurves.print_courbe import *
import pycurves.robust_differentiators as pcrbd


try:
    from RBF_tools.window_size import *
    from RBF_tools.RBF_interpolate_corrected import (
        RBFInterpolant,
    )  # this one replace the equivalent in the RBF package
    import rbf.basis

except ImportError:  # no rbf basis.installed
    pass


# ==============================================================================
# support function
# ==============================================================================


def convolve_curve(C_in, kernel, extend_method="point_mirror"):
    """
    Compute the convolution of a curve with a kernel. In practice it calls pycurves.extend_curve before doing it to preserve the boundaries of the curve.

    :param C_in: input curve.
    :type C_in: pycurves.Courbe

    :param kernel: kernel used for the convolution
    :type kelnel: numpy array

    :return: convolution of the curve and kernel
    :rtype: pycurves.Courbe

    """
    C_out = copy.deepcopy(C_in)
    C_extend = extend_curve(C_in, kernel.size, direction="both", method=extend_method)

    C_extend.y = sg.convolve(C_extend.y, kernel, mode="same")

    C_out = C_extend.project(C_in)
    return C_out


def _Holoborodko_kernek(N):
    """
    return the kernel defined by Pavel Holoborodko on his website :
    http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
    for smoothing
    """
    if N % 2 == 0:
        print("N must be odd")
    else:
        kernel = []
        m = (N - 1) / 2
        for k in np.linspace(-m, m, N):
            kernel.append(
                (3 * m - 1 - 2 * k ** 2)
                / (2 * m - 1)
                * 1.0
                / (2 ** (2 * m))
                * sp.comb(2 * m, m + k)
            )

        kernel = np.asarray(kernel)
    return kernel


def _assert_kernel_size(N):
    """
    Check if the kernel size is an odd integer otherwise if it can it transform it
    """
    if isinstance(N, float):
        N = int(N)

    if isinstance(N, int):
        if N % 2 == 0:
            N += 1
    else:
        print("N type is not accepted")
    return N


def smoothing_holoborodko(C_in, N, *args, **kwargs):
    """
    Smooth a curve using the method described here :
    http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/

    :param C_in: input curve.
    :type C_in: pycurves.Courbe

    :param N: size of the kernel used for the convolution. Must be odd
    :type N: integer

    :return: Smoothed curve
    :rtype: pycurves.Courbe
    """
    kernel = _Holoborodko_kernek(N)

    return convolve_curve(C_in, kernel, *args, **kwargs)


def mobile_mean(C_in, N, Xtimes=1, *args, **kwargs):
    """
    Smooth a curve using a mobile mean method

    :param C_in: input curve.
    :type C_in: pycurves.Courbe

    :param N: size of the kernel used for the convolution.
    :type N: integer

    :param Xtimes: Number of time to apply the smoothing. default = 1
    :type Xtimes: integer

    :return: Smoothed curve
    :rtype: pycurves.Courbe
    """
    kernel = np.zeros(N)
    kernel += 1.0 / N
    C_out = copy.deepcopy(C_in)
    for i in range(Xtimes):
        C_out = convolve_curve(C_out, kernel, *args, **kwargs)
    return C_out


def gaussian_smoothing(C_in, sigma, nb_sigma=3, *args, **kwargs):
    """
    Smooth a curve using a gaussian kernel

    :param C_in: input curve.
    :type C_in: pycurves.Courbe

    :param sigma: number of point of the sigma
    :type sigma: integer

    :param sigma: number of sigma used to build the kernel, default = 3
    :type sigma: integer


    :return: Smoothed curve
    :rtype: pycurves.Courbe
    """

    largeur = int(nb_sigma * sigma)
    C = create_courbe(
        -largeur,
        1,
        2 * largeur + 1,
        lambda x: 1 / (sigma * 2 ** 0.5 * np.pi) * np.exp(-(x ** 2) / (2 * sigma ** 2)),
    )
    kernel = C.y / np.sum(C.y)
    kernel = kernel.astype("float32")

    C_out = copy.deepcopy(C_in)
    C_out = convolve_curve(C_out, kernel, *args, **kwargs)
    return C_out


def differentiate_holoborodko(c_in, n):
    """
    First order differentiation of a Curve object using holoborodko kernel
    only work if the curve has a constant dx

    :param c_in: Courbe object to differentiate
    :type c_in: pycurves.Courbe

    :param n: number of points in the kernel used in the convolution
    :type n: integer

    :return: First derivative of c_in
    :rtype: pycurves.Courbe

    """
    kernel = pcrbd.d1coeff(n)
    return convolve_curve(c_in, kernel) / c_in.dx


def smoothing_by_RBF(C_in, sigma=0.1, window_size=1, order=1):
    """
    function to automatiquelly smooth a curve by use of sparse RBF

    :param c_in: Courbe object to differentiate
    :type c_in: pycurves.Courbe

    :param sigma: control the amount of smoothing
    :type sigma: float

    :param window_size: number of points on one size of the window (I don't remeber)
    :type window_size: integer

    :param order: order of the added polynomial terms
    :type order: integer

    :return: Smoothed curve
    :rtype: pycurves.Courbe
    """
    x_obs = np.zeros((C_in.x.size, 1))
    x_obs[:, 0] = C_in.x
    u_obs = C_in.y

    # adaptation of the number of points to the window size
    window_size = int(window_size * C_in.dx)

    basis = rbf.basis.SparseRBF(
        sympy.exp(-(rbf.basis._EPS * rbf.basis._R) ** 2), supp=window_size
    )

    I = RBFInterpolant(
        x_obs,
        u_obs,
        sigma=sigma,
        eps=eps_value(basis, ws=window_size),
        basis=basis,
        order=order,
    )

    C_out = Courbe()
    C_out.x = copy.deepcopy(C_in.x)
    C_out.y = I(x_obs)

    return C_out


def integrate(C, *args, **kwargs):
    """
    Compute the integral of a curve
    uses the scipy.integrate.simps function

    :param C: Courbe object to integrate
    :type C: pycurves.Courbe

    :return: Integral of the curve
    :rtype: pycurves.Courbe
    """
    I = copy.deepcopy(C)
    for i in range(C.x.size):
        I.y[i] = simps(C.y[0 : i + 1], C.x[0 : i + 1], *args, **kwargs)

    return I


def integrate_2(C, *args, **kwargs):
    """
    Compute the integral of a curve
    uses the scipy.integrate.simps function

    :param C: Courbe object to integrate
    :type C: pycurves.Courbe

    :return: Integral of the curve
    :rtype: pycurves.Courbe
    """
    I = copy.deepcopy(C)
    I.y[0] = 0
    I.y[1::] = cumtrapz(C.y, C.x, *args, **kwargs)

    return I


def grad(c_in, *args, **kwargs):
    """
    First order differentiation of a Curve object using the numpy gradient function
    only work if the curve has a constant dx

    :param c_in: Courbe object to differentiate
    :type c_in: pycurves.Courbe

    :return: First derivative of c_in
    :rtype: pycurves.Courbe

    """
    G = copy.deepcopy(c_in)
    G.y = np.gradient(c_in.y) / np.gradient(c_in.x)

    return G


def savgol_differentiation(C, window_length, polyorder, deriv=1, *args, **kwargs):
    """
    Differentiation of a Curve object using the Savitzky-Golay filter

    it assumes a constant step in the data

    :param c_in: Courbe object to differentiate
    :type c_in: pycurves.Courbe

    :param window_length: The length of the filter window (i.e. the number of coefficients).
                            window_length must be a positive odd integer.
    :type window_length: integer

    :param poly_order: The order of the polynomial used to fit the samples.
                    polyorder must be less than window_length.
    :type poly_order: integer

    :param deriv: The order of the derivative to compute. Default is 1.
    :type deriv: integer

    :return: Derivative of input curve
    :rtype: pycurves.Courbe

    .. Note::   all other scipy.signal.savgol_filter parameters can be passed

    """

    C_out = copy.deepcopy(C)
    C_out.y = scipy.signal.savgol_filter(
        C.y,
        window_length,
        polyorder,
        deriv=deriv,
        delta=C.x[1] - C.x[0],
        *args,
        **kwargs
    )

    return C_out
