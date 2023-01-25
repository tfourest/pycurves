# -*- coding: utf-8 -*-
"""
Define usefull function for the analysis of pycurves.courbe.Courbe instances.
"""

import pycurves as pyc

from scipy.optimize import minimize, rosen, rosen_der
import numpy as np
import pycurves.math_courbe as mc

import scipy.signal as sg

from scipy import stats
import copy


def offset_x(courbe, offset, *agrs, **kwargs):
    """
    return an offseted Courbe()
    """
    c_modif = copy.deepcopy(courbe)
    c_modif.offset_x(offset, *agrs, **kwargs)
    return c_modif


def create_courbe(x0, dx, nbstep, fonction):
    """
    Function to generate a pycurves.courbe.Courbe with a regular x step using a lambda function.

    :param x0:    x of the first point.
    :type x0: float

    :param dx: x step.
    :type dx: float

    :param nbstep: number of steps.
    :type nbstep: int

    :param fonction: function y=f(x), typically a lambda function.
    :type fonction: funtion

    :return: The curve
    :rtype: pycurves.Courbe
    """
    Ctemp = pyc.Courbe()
    Ctemp.x = np.linspace(x0, x0 + dx * (nbstep - 1), nbstep)
    Ctemp.y = copy.deepcopy(fonction(Ctemp.x))

    return Ctemp


def synchronize(C1, C2, over_sampling=1, delta_max=None):
    """
    Function to compute a likely delay between two curves using the correlation of the curves.

    :param C1: First curve
    :type C1: pycurves.Courbe

    :param C2: Second curve
    :type C2: pycurves.Courbe

    :param over_sampling: factor of the resampling of the curves to increase precision in the computation of the delay. Default is 1, which means none.
    :type over_sampling: int

    :param delta_max: Delay max, Defautl is None.
    :type delta_max: funtion

    :return: The delay of C2 versus C1
    :rtype: float

    """

    new_x = np.linspace(C1.x[0], C1.x[-1], C1.x.size * over_sampling)
    A = pyc.Courbe()
    A.x = new_x
    C1_ = C1.project(A)
    C2_ = C2.project(A)

    conv = pyc.Courbe()
    conv.y = sg.correlate(C1_.y, C2_.y, mode="full")

    dx = new_x[1] - new_x[0]

    x0 = new_x[0] - (new_x.size - 1) * dx
    conv.x = np.linspace(x0, x0 + (conv.y.size - 1) * dx, conv.y.size)
    if delta_max != None:
        conv.y[np.where(np.abs(conv.x) > delta_max)] = 0.0

    indice_max = np.where(conv.y == np.nanmax(conv.y))
    delay = conv.x[indice_max[0][0]]

    return -delay


def has_same_x(l_C):
    """
    Tests if all curves in l_C have the same .x.

    :param l_C: list of pycurves.Courbe
    :type l_C: list

    :return: True if all curves have the same .x, False otherwise
    :rtype: boolean
    """
    rt = True
    for courbe in l_C:
        rt = rt and (l_C[0].x == courbe.x).all
    return rt


def extend_curve(C, nb_points, direction="both", method="flat"):
    """
    this function expand a curves in one or two directions using the extremal values of the functions

    :param C: the curve to extend
    :type C: pycurves.Courbe

    :param nb_point:  the number of points to add on each side
    :type nb_point: int

    :param direction: direction of the extend. Can be 'right', 'left' or 'both' the default is 'both'
    :type direction: string

    :param method:     method of extension = 'flat', 'gradient',  'line_mirror', 'point_mirror'
    :type method: string

    :return: An extended version of the curve
    :rtype: pycurves.Courbe
    """

    C_extended = copy.deepcopy(C)
    step = C.x[1] - C.x[0]
    if direction == "left" or direction == "both":
        #        array = np.arange(start= C.x[0] -nb_points*step  ,stop = C.x[0],step = step)
        array = np.linspace(
            start=C.x[0] - nb_points * step, stop=C.x[0] - 1 * step, num=nb_points
        )
        C_extended.x = np.append(array, C.x)
        if method == "flat":
            C_extended.y = np.append(array * 0 + C.y[0], C.y)
        elif method == "gradient":
            step_y = C.y[1] - C.y[0]
            array = np.arange(
                start=C.y[0] - nb_points * step_y, stop=C.y[0], step=step_y
            )
            C_extended.y = np.append(array, C.y)
        elif method == "line_mirror":
            array = C.y[1 : nb_points + 1]
            C_extended.y = np.append(array[::-1], C.y)
        elif method == "point_mirror":
            array = 2 * C.y[0] - C.y[1 : nb_points + 1]
            C_extended.y = np.append(array[::-1], C.y)

    if direction == "right" or direction == "both":
        array = np.linspace(
            start=C.x[-1] + step, stop=C.x[-1] + (nb_points + 1) * step, num=nb_points
        )
        C_extended.x = np.append(C_extended.x, array)
        if method == "flat":
            C_extended.y = np.append(C_extended.y, array * 0 + C.y[-1])
        elif method == "gradient":
            step_y = C.y[-1] - C.y[-2]
            array = np.arange(
                C.y[-1] + step_y, stop=C.y[-1] + (nb_points + 1) * step_y, step=step_y
            )
            C_extended.y = np.append(C_extended.y, array)
        elif method == "line_mirror":
            array = C.y[-nb_points - 1 : -1]
            C_extended.y = np.append(C_extended.y, array[::-1])
        elif method == "point_mirror":
            array = 2 * C.y[-1] - C.y[-nb_points - 1 : -1]

            C_extended.y = np.append(C_extended.y, array[::-1])

    return C_extended


def compute_non_linearity_threshold(
    curve, diff_threshold=5, pc_ignored_pts=5, return_total_error=False
):
    """
    Method to compute the non-linearity threshold of a curves, i.e. the limit from which it is non-linear
    It is based on :
    M. Castres, J. Berthe, E. Deletombe, M. Brieu
    Experimental evaluation of the elastic limit of carbon-fibre reinforced epoxy
    composites under a large range of strain rate and temperature conditions
    Strain. 2017; 53:e12248.

    :param curve:                 Curve on which the non-linearity threshold is search
    :type curve:     Courbe

    :param diff_threshold:  Difference in percentage for which the non-linearity threshold will be considered (default = 5)
    :type diff_threshold:
           float
    :param pc_ignored_pts:  percentage of ingnored points at the begining (default = 5)
    :type pc_ignored_pts: float

    :param return_total_error:     True if a curve of the error for each points should be return
    :type return_total_error:           boolean

    :return :   x and y values of the threshold, Curve of the error computed for all points
    :rtype:     array [x,y], pycurves.Courbe

    """

    nb_pts_ini = int(curve.x.size * (pc_ignored_pts / 100))

    total_error = pyc.Courbe()
    total_error.x = copy.deepcopy(curve.x)
    total_error.y = np.zeros_like(curve.x)

    i = 0
    for nb_pts in range(2, curve.x.size - 1):
        new_curve = copy.deepcopy(curve)
        new_curve.x = new_curve.x[0:nb_pts]
        new_curve.y = new_curve.y[0:nb_pts]

        slope, intercept, r_value, p_value, std_err = stats.linregress(
            new_curve.x, new_curve.y
        )
        f = lambda x: slope * x + intercept
        numerator = np.abs(curve.y[nb_pts + 1] - f(curve.x[nb_pts + 1]))
        denominator = curve.y[nb_pts + 1]
        total_error.y[nb_pts + 1] = np.sum(
            np.divide(
                numerator,
                denominator,
                where=denominator != 0,
                out=np.zeros_like(numerator),
            )
        )

    indices = np.where(total_error.y[nb_pts_ini::] > diff_threshold / 100)
    if indices != np.array([]):
        valeurs_x_y = np.array(
            [curve.x[nb_pts_ini + indices[0][0]], curve.y[nb_pts_ini + indices[0][0]]]
        )
    else:
        valeurs_x_y = np.array([])

    if return_total_error == False:
        return valeurs_x_y
    else:
        return valeurs_x_y, total_error
