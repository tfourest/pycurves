# -*- coding: utf-8 -*-
"""
Classic frequency filter functions to easily filter curves.

"""

from numpy import zeros

from pycurves.frequency_tools.useful_tools import _init_courbe_TIrevin
from math import pi, cos, sin, tan
from scipy.signal import lfilter

from pycurves.frequency_tools.fft_courbe import *

import scipy.signal as signal


def generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype):
    """
    A generic filter that is called by the other filters for simpler maintainability
    
    :param courbe: the courbe on which you want to apply the filter
    :type courbe: pycurves.Courbe

    :param freq_cut: scalar or (f1,f2) for the bandpass and bandstop filters  
              cutoff frequency must be < Nyquist frequency
    :type freq_cut: scalar or tupple

    :param filter_type: the type of the filter. The available types are : 
        'butterworth'   maximum constancy for the gain
        'cheby1'        Chebyshev type 1 filter has ripple in the gain frequency so unrecommended
        'cheby2'        Chebyshev type 1 filter has no ripple in the bandpass 
                          but it has ripple in the bandstop (not very important)
                          it start to significantly disminish the gain in the bandpass before freq_cut
        'bessel'        maximum constancy for the delay
    :type filter_type: string
    
    :param filtfilt: True if the filter must be apply forward and backward
    :type filtfilt: Boolean   
    
    :param order: the order of the filter
    :type order: integer

    :param btype: the band type of the filter 'lowpass', 'highpass', 'bandpass' and 'bandstop'
    :type btype: string
    
    :return: the filtered curve
    :rtype: pycurves.Courbe 
    

    """

    courbe_filt = _init_courbe_TIrevin(courbe)

    # First, design the filter
    Wn = freq_cut / (
        0.5 * courbe_filt.sr
    )  # cutoff frequency , Wn is normalized from 0 to 1, where 1 is the Nyquist frequency

    if filter_type == "butterworth":
        fonction = lambda x, y, z: signal.butter(x, y, z)
    elif filter_type == "cheby1":
        fonction = lambda x, y, z: signal.cheby1(x, y, z)
    elif filter_type == "cheby2":
        fonction = lambda x, y, z: signal.cheby2(x, y, z)
    elif filter_type == "bessel":
        fonction = lambda x, y, z: signal.bessel(x, y, z)

    B, A = fonction(order, Wn, btype)
    # Second, apply the filter
    if filtfilt == True:
        new_y = signal.filtfilt(B, A, courbe.y)
    else:
        new_y = signal.lfilter(B, A, courbe.y)

    courbe_filt.y = new_y
    return courbe_filt


###############################################################################
######             Butterworth filters quick calling         ##################
###############################################################################
def LP_butterworth(courbe, freq_cut, order=2):
    """
    Apply a low pass butterworth filter to a curve.
    
    :param courbe: the courbe on which you want to apply the filter
    :type courbe: pycurves.Courbe
    
    :param freq_cut: scalar cutoff frequency must be < Nyquist frequency
    :type freq_cut: scalar
    
    :param order: the order of the filter (default = 2)
    :type order: integer  
    
    :return: the filtered curve
    :rtype: pycurves.Courbe     

    """
    filter_type = "butterworth"
    filtfilt = False
    btype = "lowpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def HP_butterworth(courbe, freq_cut, order=2):
    """
    High pass butterworth filter 
    
    :param courbe: the courbe on which you want to apply the filter
    :type courbe: pycurves.Courbe
    
    :param freq_cut: scalar cutoff frequency must be < Nyquist frequency
    :type freq_cut: scalar
    
    :param order: the order of the filter (default = 2)
    :type order: integer
    
    :return: the filtered curve
    :rtype: pycurves.Courbe        
    """
    filter_type = "butterworth"
    filtfilt = False
    btype = "highpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def BP_butterworth(courbe, freq_cut, order=2):
    """
    Bandpass butterworth filter 
        
    :param courbe: the courbe on which you want to apply the filter
    :type courbe: pycurves.Courbe
    
    :param freq_cut: (f1,f2), cutoff frequency must be < Nyquist frequency
    :type freq_cut: tupple
    
    :param order: the order of the filter (default = 2)
    :type order: integer
    
    :return: the filtered curve
    :rtype: pycurves.Courbe        
    """
    filter_type = "butterworth"
    filtfilt = False
    btype = "bandpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def BS_butterworth(courbe, freq_cut, order=2):
    """
    Bandstop butterworth filter 
         
    :param courbe: the courbe on which you want to apply the filter
    :type courbe: pycurves.Courbe
    
    :param freq_cut: (f1,f2), cutoff frequency must be < Nyquist frequency
    :type freq_cut: tupple
    
    :param order: the order of the filter (default = 2)
    :type order: integer
    
    :return: the filtered curve
    :rtype: pycurves.Courbe      
    """
    filter_type = "butterworth"
    filtfilt = False
    btype = "bandstop"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def LPZ_butterworth(courbe, freq_cut, order=2):
    """
    Low pass butterworth filter with no delay, the filter is applied two times forward and backward
    
    Same inputs as LP_butterworth
    """
    filter_type = "butterworth"
    filtfilt = True
    btype = "lowpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def HPZ_butterworth(courbe, freq_cut, order=2):
    """
    High pass butterworth filter with no delay, the filter is applied two times forward and backward
    
    Same inputs as HP_butterworth
    """
    filter_type = "butterworth"
    filtfilt = True
    btype = "highpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def BPZ_butterworth(courbe, freq_cut, order=2):
    """
    Bandpass butterworth filter with no delay, the filter is applied two times forward and backward
    
    Same inputs as BP_butterworth
    """
    filter_type = "butterworth"
    filtfilt = True
    btype = "bandpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def BSZ_butterworth(courbe, freq_cut, order=2):
    """
    Bandstop butterworth filter with no delay, the filter is applied two times forward and backward
    
    Same inputs as BS_butterworth
    """
    filter_type = "butterworth"
    filtfilt = True
    btype = "bandstop"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
######                  Bessel filters quick calling         ##################
###############################################################################
def LP_bessel(courbe, freq_cut, order=2):
    """
    Low pass bessel filter 
    
    :param courbe: the courbe on which you want to apply the filter
    :type courbe: pycurves.Courbe
    
    :param freq_cut: scalar cutoff frequency must be < Nyquist frequency
    :type freq_cut: scalar
    
    :param order: the order of the filter (default = 2)
    :type order: integer
    
    :return: the filtered curve
    :rtype: pycurves.Courbe    
    
    """
    filter_type = "bessel"
    filtfilt = False
    btype = "lowpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def HP_bessel(courbe, freq_cut, order=2):
    """
    High pass bessel filter 
    
    :param courbe: the courbe on which you want to apply the filter
    :type courbe: pycurves.Courbe
    
    :param freq_cut: scalar cutoff frequency must be < Nyquist frequency
    :type freq_cut: scalar
    
    :param order: the order of the filter (default = 2)
    :type order: integer
    
    :return: the filtered curve
    :rtype: pycurves.Courbe      
    
    """
    filter_type = "bessel"
    filtfilt = False
    btype = "highpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def BP_bessel(courbe, freq_cut, order=2):
    """
    Bandpass bessel filter 
    
    :param courbe: the courbe on which you want to apply the filter
    :type courbe: pycurves.Courbe
    
    :param freq_cut: (f1,f2), cutoff frequency must be < Nyquist frequency
    :type freq_cut: tupple
    
    :param order: the order of the filter (default = 2)
    :type order: integer
    
    :return: the filtered curve
    :rtype: pycurves.Courbe      
    """
    filter_type = "bessel"
    filtfilt = False
    btype = "bandpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def BS_bessel(courbe, freq_cut, order=2):
    """
    Bandstop bessel filter
    
    :param courbe: the courbe on which you want to apply the filter
    :type courbe: pycurves.Courbe
    
    :param freq_cut: (f1,f2), cutoff frequency must be < Nyquist frequency
    :type freq_cut: tupple
    
    :param order: the order of the filter (default = 2)
    :type order: integer
    
    :return: the filtered curve
    :rtype: pycurves.Courbe      
    
    """
    filter_type = "bessel"
    filtfilt = False
    btype = "bandstop"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def LPZ_bessel(courbe, freq_cut, order=2):
    """
    Low pass bessel filter with no delay, the filter is applied two times forward and backward
    
    Same inputs as LP_bessel
    """
    filter_type = "bessel"
    filtfilt = True
    btype = "lowpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def HPZ_bessel(courbe, freq_cut, order=2):
    """
    High pass bessel filter with no delay, the filter is applied two times forward and backward
    
    Same inputs as HP_bessel
    """
    filter_type = "bessel"
    filtfilt = True
    btype = "highpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def BPZ_bessel(courbe, freq_cut, order=2):
    """
    High pass bessel filter with no delay, the filter is applied two times forward and backward
    
    Same inputs as BP_bessel

    """
    filter_type = "bessel"
    filtfilt = True
    btype = "bandpass"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)


###############################################################################
def BSZ_bessel(courbe, freq_cut, order=2):
    """
    High pass bessel filter with no delay, the filter is applied two times forward and backward
    
    Same inputs as BS_bessel   
    """
    filter_type = "bessel"
    filtfilt = True
    btype = "bandstop"
    return generic_filter(courbe, freq_cut, filter_type, filtfilt, order, btype)
