# -*- coding: utf-8 -*-
"""
This module provide wrappers of scipy.fftpack fft and ifft functions for direct use on pycurves.courbe instances.
"""


import pycurves as pyc

import math

from numpy import zeros, argmax
from numpy import array
import numpy as np

import scipy.fftpack


########################################################################
def fft(courbe, complexe=False, remove_mean=False, normalized=False, length_out="half"):
    """
    return the fast fourrier transfert funtion of a pycurves.Courbe
    
    :param courbe: input curve of the FFT
    :type courbe: pycurves.Courbe

    :param complexe: indicate if the curve.y are complex numbers. Default is False.
    :type complexe: boolean

    :param remove_mean: indicate if the curve mean value should be removed. Default is False.
    :type remove_mean: boolean

    :param normalized: indicate if output must be normalised. Default is False.
    :type normalized: boolean    

    :param length_out: indicate if the 'half or 'full' output should be returned. Default is 'half''
    :type length_out: string
    
    :return: magnitude, phase of the input curve
    :rtype: pycurves.Courbe 
        
    """

    n = courbe.y.size

    timestep = courbe.x[1] - courbe.x[0]
    xf = scipy.fftpack.fftfreq(n, d=timestep)
    yf = scipy.fftpack.fft(courbe.y)

    c_magnitude_fft = pyc.Courbe()
    c_phase_fft = pyc.Courbe()

    if complexe == True:
        c_fft = pyc.Courbe()
        c_fft.x = xf
        c_fft.y = yf
        return c_fft
    else:

        if length_out == "half":
            c_magnitude_fft.x = xf[: n // 2]
            c_magnitude_fft.y = np.absolute(yf)[: n // 2]

            c_phase_fft.x = xf[: n // 2]
            c_phase_fft.y = np.unwrap(np.angle(yf)[: n // 2])

        elif length_out == "full":
            c_magnitude_fft.x = xf
            c_magnitude_fft.y = np.absolute(yf)

            c_phase_fft.x = xf
            c_phase_fft.y = np.unwrap(np.angle(yf))

        else:
            print("no option " + length_out + " for length_out")

        return c_magnitude_fft, c_phase_fft


########################################################################
def ifft(c_mag, c_phase, complexe=False, length_in="half"):
    """
    return the inverse fast fourrier transfert funtion of a pycurves.Courbe() instance


    :param c_mag: magnitude of the signal on which to perform ifft
    :type c_mag: pycurves.Courbe

    :param c_phase: phase of the signal on which to perform ifft
    :type c_phase: pycurves.Courbe

    :param length_in: indicate if the 'half or 'full' input is given to the function. Default is 'half''
    :type length_in: string
    
    :return: inverse FFT of the signal
    :rtype: pycurves.Courbe 
    
    .. Warning:: When performing the ifft or a fft result, if length_in='half' is used and the initial signal has an odd lenght a small error is done.

    """
    if length_in == "half":
        # ajout de la seconde partie de la fft que est redondante mais permet de conserver
        # le même nombre de points que le signal original une fois ifft fait
        #
        mag_complet = copy.deepcopy(c_mag.y)
        phase_complet = copy.deepcopy(c_phase.y)
        complex_half = mag_complet * np.exp(phase_complet * 1j)

        complex_fft_terms = np.concatenate(
            (complex_half, array([complex_half[-1]]), np.conj(complex_half[::-1]))
        )
        complex_fft_terms = np.delete(complex_fft_terms, -1)

        n = complex_fft_terms.size

    elif length_in == "full":
        complex_fft_terms = c_mag.y * np.exp(c_phase.y * 1j)
        n = c_mag.y.size

    else:
        print("no option " + length_in + " for length_in")

    c_ifft = pyc.Courbe()
    c_ifft.y = scipy.fftpack.ifft(complex_fft_terms)

    f_max = np.max(c_mag.x)
    sampling_freq = 2 * f_max
    if n % 2 == 0:
        sampling_freq = f_max * 2 * n / (n - 2)
    else:
        sampling_freq = f_max * 2 * n / (n - 1)

    temporaire = create_courbe(0, 1.0 / sampling_freq, n, lambda x: x)
    c_ifft.x = temporaire.x

    return c_ifft


########################################################################
def correlate(c1, c2, mode="full"):
    """
    Calculate the cross-correlation of two picurves.Courbe instance.
    the chosen mode is same for scipy.correlate is 'same' it impose to have the same x for the two curves
    return a correlation Courbe with c_out.x = c1.x
    
    mode = full or same
    """
    c_out = pyc.Courbe()
    #    if normalize:
    a = (c1.y - np.mean(c1.y)) / (np.std(c1.y) * len(c1.y))
    v = (c2.y - np.mean(c2.y)) / np.std(c2.y)
    c_out.y = scipy.signal.correlate(a, v, mode=mode)
    #    x1 = x1/x1.std()

    # x2 = x2/x2.std() and then as you did it
    #    print(c1.x)
    dt = c1.x[1] - c1.x[0]
    #    print('dt = ',dt)

    if mode == "full":
        X = np.linspace(-c_out.y.size / 2 * dt, c_out.y.size / 2 * dt, c_out.y.size)
        c_out.x = X

    elif mode == "same":
        c_out.x = c1.x

    else:
        print("no other mode available than full or same ")

    return c_out


########################################################################
def autocorrelate(c1):
    """
    calculate the autocorrelation of a Courbe() object
    In practice it calls correlate(c1,c1)
    """
    return correlate(c1, c1)


########################################################################
def covariance(c1, c2):
    """
    Return the covariance of two Courbe(),
    the relation cross_correlation(c1,c2) = 1/(mean(c1)*mean(c2)) * covariance(c1,c2) is used
    """
    c_out = correlate(c1, c2) * np.mean(c1.y) * np.mean(c2.y)
    return c_out


if __name__ == "__main__":

    from useful_tools import *

    N = 100.0
    E01 = pyc.Courbe()
    E01.x = np.linspace(0.0, 1, 100)
    #    E01.y = np.sin(50.0 * 2.0*np.pi*E01.x)
    E01.y = np.concatenate((np.zeros(40), np.ones(20), np.zeros(40)))
    E01.label = "E01"

    E02 = pyc.Courbe()
    E02.x = np.linspace(0.0, 1, 100)
    #    E01.y = np.sin(50.0 * 2.0*np.pi*E01.x)
    E02.y = np.concatenate((np.zeros(40), np.linspace(0, 1, 20)[::-1], np.zeros(40)))

    E02.label = "E02"

    #    plotXcourbes([E01,E02 ])

    CC = correlate(E01, E02)
    CC.label = "cross_correlate_E01*E02"

    AC = autocorrelate(E01)
    AC.label = "auto_correlate_E01"

    CCovar = covariance(E01, E02)
    CCovar.label = "covariance_E01_E02"
    plotXcourbes([E01, E02])

    plotXcourbes([E01, E02, CC, AC, CCovar])

#
#
