# -*- coding: utf-8 -*-
"""
Some function to smooth and differentiate curves using splines.
"""
from scipy.interpolate import LSQUnivariateSpline

from scipy.interpolate import UnivariateSpline

import copy
import numpy as np


def spline_fitting(curve, N, *args, **kwargs):
    """
    Fit a spline on a Courbe object.
    It is a wrapper around the  scipy.interpolate.LSQUnivariateSpline

    :param curve: Courbe object to integrate
    :type curve: pycurves.Courbe

    :param N: step for the number of knot to use. For exemple if N=100, one knot is used every 100 points.
    The larger the value the larger the smoothing.
    :type N: pycurves.Courbe

    :return: Smoothed curve
    :rtype: pycurves.Courbe

    """
    new_curve = copy.deepcopy(curve)

    x = np.arange(curve.size)
    knots = x[:: 2 * N + 1][1:]
    clspline = LSQUnivariateSpline(x, curve.y, knots)
    new_curve.y = clspline(x)

    return new_curve


def spline_derivative(curve, SF=None, order=1, *args, **kwargs):
    """
    Fit a spline on a Courbe object and return the derivative of the curve.
    It is a wrapper around the  scipy.interpolate.UnivariateSpline

    :param curve: Courbe object to integrate
    :type curve: pycurves.Courbe

    :param SF: Smoothing factor. The larger the value the larger the smoothing.
    :type SF: pycurves.Courbe

    :return: Smoothed curve
    :rtype: pycurves.Courbe

    """
    new_curve = copy.deepcopy(curve)
    w = np.isnan(new_curve.y)
    new_curve.y[w] = 0.0
    spl = UnivariateSpline(new_curve.x, new_curve.y, w=~w, *args, **kwargs)

    if SF != None:
        spl.set_smoothing_factor(SF)
    new_spl = spl.derivative(n=order)
    new_curve.y = new_spl(new_curve.x)

    return new_curve


if __name__ == "__main__":
    from pycurves import *

    plt.close("all")
    c = create_courbe(0, 2, 200, lambda x: 2 * x - 0.002 * x ** 2)
    c.plot()
    noise = np.random.normal(0, 20, 200)
    c2 = copy.deepcopy(c)
    c2.y += noise
    c2.label = "C + noise"
    plotXcourbes([c, c2])

    c3 = spline_fitting(c2, 1000000)
    c3.label = "spline_interp"
    plotXcourbes([c, c2, c3])

    c_deriv = create_courbe(0, 2, 200, lambda x: 2 - 0.004 * x)
    c_deriv.label = "c_deriv"

    c3_deriv = spline_derivative(c2, 1000000)
    c3_deriv.label = "spline_deriv"
    plotXcourbes([c_deriv, c3_deriv])
#    c=Courbe()
#    c.x = np.linspace(-3, 3, 50)
#    c.y = np.exp(-c.x**2) + 0.1 * np.random.randn(50)
#    c.plot()
#
#    c2 = spline_fitting(c,0.1)
#    plotXcourbes([c,c2])
