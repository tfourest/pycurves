# -*- coding: utf-8 -*-
"""
A function to easily compute the time of images.
"""

import pycurves as pyc
import numpy as np
import scipy.misc
from pycurves.os_tools import *
import copy


def nearest_interp(xi, x, y):
    idx = np.abs(x - xi[:, None])
    return y[idx.argmin(axis=1)]


def fast_nearest_interp(xi, x, y):
    """Assumes that x is monotonically increasing!!."""
    # Shift x points to centers
    spacing = np.diff(x) / 2
    x = x + np.hstack([spacing, spacing[-1]])
    # Append the last point in y twice for ease of use
    y = np.hstack([y, y[-1]])
    return y[np.searchsorted(x, xi)]


def temps_Camera(inp, name=None):
    """
    Compute a numpy array which contain the time for each image. Here we assume that the frequency is constant.

    the input is a dictionnary containing :
        -t_trig = the time of the trigger
        -freq = the frequency of the camera
        -im_ini = the number of the first image used
        -im_fin ou nb_images  the number of the final image used or the number of images used
        -skip = the skip factor used to resample the images
        -start = 0 ou 1  the number of the image taken at t = t_trig (0 for fastcam, 1 for CEDIP)

    name is a dictionnary for the correspondance if other names have been used (better not to use it)

    return : numpy array with the time corresponding to each image
    """

    inputs = copy.deepcopy(inp)
    if name != None:
        alias(inputs, name)

    #    if inputs.has_key('nb_images') :
    #        print('nb_images : ',inputs['nb_images'])

    if not "skip" in inputs:
        inputs["skip"] = 1
    if "im_fin" in inputs and not "nb_images" in inputs:
        inputs["nb_images"] = inputs["im_fin"] - inputs["im_ini"] + 1
    #    print('nb_images : ',inputs['nb_images'])

    t_ini = inputs["t_trig"] + (inputs["im_ini"] - inputs["start"]) / inputs["freq"]
    t_fin = (
        inputs["t_trig"]
        + (inputs["im_ini"] - inputs["start"]) / inputs["freq"]
        + (inputs["nb_images"] - 1) / (inputs["freq"] / inputs["skip"])
    )
    nb_pts = inputs["nb_images"]

    temps_camera = np.linspace(t_ini, t_fin, nb_pts)

    return temps_camera
