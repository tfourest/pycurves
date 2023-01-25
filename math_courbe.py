# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 11:45:00 2017

@author: tfourest

02/06/2021 Should I keep this, I haven't used this in years
"""


from pycurves.courbe import *
import numpy as np
import scipy.signal


def log(C):
    C_temp = copy.deepcopy(C)  # Courbe()
    C_temp.x = C.x
    C_temp.y = np.log(C.y)
    return C_temp


def log10(C):
    C_temp = copy.deepcopy(C)  # Courbe()
    C_temp.x = C.x
    C_temp.y = np.log10(C.y)
    return C_temp


def Cmean(lCourbe):
    """
    moyenne d'une liste de courbes
    """
    somme = 0
    for i in lCourbe:
        somme = i + somme
    somme = somme / len(lCourbe)
    return somme


def Cmax(lcourbe):
    """
    Courbe avec les maximum d'une liste de courbes
    la fonction est particulièrement compliqué car les courbes n'ont pas besoin d'avoir les mêmes X
    """
    X = lcourbe[1].x
    for i in lcourbe:
        X = np.append(X, i.x)
    X = np.unique(X)
    X = np.sort(X)
    #    print(X)
    new_lcourbe = []
    for i in lcourbe:
        new_lcourbe.append(i.project(X))
    #    print(len(new_lcourbe))
    #    print(new_lcourbe)
    Y = copy.deepcopy(new_lcourbe[1].y)
    for i in new_lcourbe:
        Y = np.append(Y, i.y, axis=0)
    Y = np.reshape(Y, (len(new_lcourbe) + 1, len(X)))
    maximum = np.amax(Y, axis=0)

    Ctemp = Courbe()
    Ctemp.x = X
    Ctemp.y = maximum
    return Ctemp


def Cmin(lcourbe):
    """
    Courbe avec les minimums d'une liste de courbes
    la fonction est particulièrement compliqué car les courbes n'ont pas besoin d'avoir les mêmes X
    """
    X = lcourbe[1].x
    for i in lcourbe:
        X = np.append(X, i.x)
    X = np.unique(X)
    X = np.sort(X)
    #    print(X)
    new_lcourbe = []
    for i in lcourbe:
        new_lcourbe.append(i.project(X))
    #    print(len(new_lcourbe))
    #    print(new_lcourbe)
    Y = copy.deepcopy(new_lcourbe[1].y)
    for i in new_lcourbe:
        Y = np.append(Y, i.y, axis=0)
    Y = np.reshape(Y, (len(new_lcourbe) + 1, len(X)))
    minimum = np.amin(Y, axis=0)

    Ctemp = Courbe()
    Ctemp.x = X
    Ctemp.y = minimum
    return Ctemp


def sin(C):
    C_temp = Courbe()
    C_temp.x = C.x
    C_temp.y = np.sin(C.y * np.pi / 180.0)
    return C_temp


def cos(C):
    C_temp = Courbe()
    C_temp.x = C.x
    C_temp.y = np.cos(C.y * np.pi / 180.0)
    return C_temp


def RMS(C1, C2):
    """
    method to compute the RMS difference of two Courbes() objects
    It is better if the two curves have the same self.x
    """
    diff = C1 - C2
    #    print (type(diff))
    return np.sqrt(np.mean(np.square(np.asarray(diff.y))))
