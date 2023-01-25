# -*- coding: utf-8 -*-
"""
Some usefull functions when using frequency filters or looking at the freqeuncy signal of curves.
"""

import pycurves as pyc
import numpy as np
import copy
import scipy.stats as stats

def pow2up(x):
    """
    Return the power of two superior to x
    """
    return 1 << (x - 1).bit_length()


########################################################################
def zeros_padding(courbe, n=0, pow2=True, left=False):
    """
    Fonction to do zeros padding on a curve before doing the fft
    
    :param courbe: input curve of the FFT
    :type courbe: pycurves.Courbe

    :param n: number of zeros that are added
    :type n: integer

    :param pow2: if pow2 == True then n is majored to obtain a final length of the curve equal to a power of 2
    :type pow2: boolean

    :param left: by default the zeros are added to the right side of the signal. If left == True the zeros are added to the left side.
    :type left: boolean    

    :return: the curve with additional zeros.
    :rtype: pycurves.Courbe 
   
    """

    if pow2 == True:
        n = pow2up(n + courbe.x.size) - courbe.x.size
    dt = courbe.x[1] - courbe.x[0]

    if left == False:
        zeros = pyc.Courbe()
        zeros.x = np.linspace(courbe.x[-1] + dt, courbe.x[-1] + dt * (n), n)
        zeros.y = zeros.x * 0

        cb = copy.deepcopy(courbe)
        cb.append(zeros)
    elif left == True:
        zeros = Courbe()
        zeros.x = np.linspace(courbe.x[0] - dt * n, courbe.x[0] - dt, n)
        zeros.y = zeros.x * 0

        cb = copy.deepcopy(courbe)
        zeros.append(cb)
        cb = zeros
    return cb


########################################################################
def mirror(courbe):
    """
    Mirror a curve, useful in some case before doing a fft or filtering
    
    :param courbe: curve to mirror
    :type courbe: pycurves.Courbe
    
    :return: the symetry of the curves.
    :rtype: pycurves.Courbe 

    """
    dt = courbe.x[1] - courbe.x[0]
    c_original = copy.deepcopy(courbe)
    c_mirror = pyc.Courbe()
    c_mirror.y = np.fliplr([c_original.y])[0]
    c_mirror.x = c_original.x + c_original.x[-1] + dt - c_original.x[0]

    c_original.append(c_mirror)
    return c_original


########################################################################
def _init_courbe_TIrevin(courbe):

    # Calls Tom Irvine method to get some signal statistics
    # and to check the sample rate

    courbe_out = copy.deepcopy(courbe)

    sr, dt, mean, sd, rms, skew, kurtosis, dur = _signal_stats(
        courbe.x, courbe.y, courbe.x.size
    )
    sr, dt = _sample_rate_check(courbe.x, courbe.y, courbe.x.size, sr, dt)

    courbe_out.mean = mean
    courbe_out.sd = sd
    courbe_out.rms = rms
    courbe_out.skew = skew
    courbe_out.kurtosis = kurtosis
    courbe_out.dur = dur
    courbe_out.sr = sr
    courbe_out.dt = dt

    return courbe_out


########################################################################


def _signal_stats(a, b, num):
    # From Tom Irvine code, modified so it doesn't talk
    # a is the time column.
    # b is the amplitude column.
    # num is the number of coordinates
    # Return
    #       sr - sample rate
    #       dt - time step
    #     mean - average
    #       sd - standard deviation
    #      rms - root mean square
    #     skew - skewness
    # kurtosis - peakedness
    #      dur - duration

    ave = np.mean(b)

    dur = a[num - 1] - a[0]

    dt = dur / float(num - 1)
    sr = 1 / dt

    rms = np.sqrt(np.var(b))
    sd = np.std(b)

    skewness = stats.skew(b)
    kurtosis = stats.kurtosis(b, fisher=False)
    return sr, dt, ave, sd, rms, skewness, kurtosis, dur


########################################################################


def _sample_rate_check(a, b, num, sr, dt):
    # From Tom Irvine code, modified so it doesn't talk
    # tells you if dtmin is too different from dtmax

    dtmin = 1e50
    dtmax = 0

    dtmin = np.min(a[1:-1] - a[0:-2])
    dtmax = np.max(a[1:-1] - a[0:-2])
    #    for i in range(1, num-1):
    #        if (a[i]-a[i-1])<dtmin:
    #            dtmin=a[i]-a[i-1];
    #            if (a[i]-a[i-1])>dtmax:
    #                dtmax=a[i]-a[i-1];

    srmax = float(1 / dtmin)
    srmin = float(1 / dtmax)

    if (srmax - srmin) > 0.01 * sr:
        print(" ")
        print(" Warning: sample rate difference ")
        print("  dtmin = %8.4g sec" % dtmin)
        print("     dt = %8.4g sec" % dt)
        print("  dtmax = %8.4g sec \n" % dtmax)
    return sr, dt


########################################################################
