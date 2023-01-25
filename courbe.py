# -*- coding: utf-8 -*-
"""
The curve package is the main one in pycurves
"""


import scipy

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:  # ... no graphic display
    pass

import pycurves.os_tools as otls

# from scipy import interpolate
import scipy.signal

from scipy.interpolate import interp1d

import os.path

from scipy import stats
import copy

from math import log
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline
from pycurves.print_courbe import *


###############################################################################


class Courbe:
    """
    The Courbe object is the main object of the pycurves package.

    Attributes
    ----------
    x : numpy array
        array containing the x values of the curves
    y : numpy array
        array containing the y values of the curves
    label : string
        name to use for label when plotting the curve

    """

    def __init__(
        self, file=None, xlabel=None, ylabel=None, label=None, *args, **kwargs
    ):
        """

        """
        if file != None:
            if file[-3] + file[-2] + file[-1] == "csv":
                self.x, self.y = np.loadtxt(
                    file, unpack="true", skiprows=1, delimiter=","
                )
            elif file[-3] + file[-2] + file[-1] == "npy":
                self.load(file)
            elif file[-3] + file[-2] + file[-1] in ["txt", "asc"]:
                self.x, self.y = np.loadtxt(file, unpack="true", *args, **kwargs)
            else:
                print("pycurves doesn't know this file format")
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.show_approx = False
        if label != None:
            self.label = label
        elif file != "":
            self.label = file

    # ------------------------------------------------------------------------------
    def append(self, A):
        """
        Concatenate two Courbe() object one after anotherl'une à la suite de l'autre
        """
        if not hasattr(self, "x") and hasattr(A, "x"):
            #            print('icic')
            self.x = copy.deepcopy(A.x)
            self.y = copy.deepcopy(A.y)

        elif not hasattr(A, "x"):
            pass
        else:
            Ctemp = copy.deepcopy(self)
            Ctemp.x = np.concatenate((self.x, A.x))
            Ctemp.y = np.concatenate((self.y, A.y))
            self.x = Ctemp.x
            self.y = Ctemp.y

    # ------------------------------------------------------------------------------

    def croissanceX(self):
        """
        Function that only keeps the couples of (x,y) with growing x. The original Courbe()
        instance is modified
        """
        X = []
        Y = []
        for i in range(len(self)):
            if i == 0:
                X.append(self.x[i])
                Y.append(self.y[i])
            else:
                if self.x[i] > X[-1]:
                    #                    print('plus grand')
                    X.append(self.x[i])
                    Y.append(self.y[i])
                else:
                    pass
        #                    print('plus petit')

        self.x = np.asarray(X)
        self.y = np.asarray(Y)

    def croissanceY(self, keep_x=True):
        """
        Ensure strictly growing values by changing decreasing values by the previous ones

        :param keep_x: True if the .x values are unchanged and the .y are recopied. False if only the growing (x,y) couples are kept.
        :type keep_x: boolean

        """
        if keep_x == True:
            for i in range(1, self.x.size):
                if self.y[i] < self.y[i - 1]:
                    self.y[i] = self.y[i - 1]
        else:
            x = [self.x[0]]
            y = [self.y[0]]
            previous = self.y[0]
            for i in range(1, self.x.size):
                if self.y[i] > previous:
                    x.append(self.x[i])
                    y.append(self.y[i])
                    previous = self.y[i]
            self.x = np.array(x)
            self.y = np.array(y)

    # ------------------------------------------------------------------------------
    def cut(self, method=None, retour=False):
        """
        Method for cutting a curve

        :param method:  a dictionnary containig the chosen method keyword and the
                        corresponding parameters
        :type method: dict

        :param retour:  a dictionnary containig the chosen method keyword and the
                        corresponding parameters
        :type retour: boulean

        :return: if retour is True then a new curves is returned
        :rtype: Courbe() object

        +----------------+---------------------------------------------------------------------+
        | method key     |    Description                                                      |
        +================+=====================================================================+
        | 'inf'          |defines the method for the left hand cut. Possibilities are :        |
        |                |                                                                     |
        |                |- 'x'  cut at a x given by 'x1'                                      |
        |                |- 'ymax' cut at y corresponding to ymax                              |
        |                |- 'y' : cut at the values given with the key 'y1'                    |
        +----------------+---------------------------------------------------------------------+
        | 'sup'          |defines the method for the right hand cut, Possibilities are :       |
        |                |                                                                     |
        |                |- 'x'  cut at a x given by 'x2'                                      |
        |                |- 'ymax' cut at y corresponding to ymax                              |
        |                |- 'y' : cut at the values given with the key 'y2'                    |
        +----------------+---------------------------------------------------------------------+
        | 'keep_inf'     | defines the number of point to keep befor the left hand cut         |
        +----------------+---------------------------------------------------------------------+
        | 'keep_sup'     | defines the number of point to keep after the right hand cut        |
        +----------------+---------------------------------------------------------------------+

        .. todo:: some methods could be improuved by using np.where rather than loops


        """

        if "keep_inf" not in method:
            method["keep_inf"] = 0
        if "keep_sup" not in method:
            method["keep_sup"] = 0

        variable = self.y

        if "courbe" in method:
            borne_inf = method["courbe"].borne_cut_indice[0]
            borne_sup = method["courbe"].borne_cut_indice[1]

        else:
            borne_inf = 1
            if "inf" in method:
                if method["inf"] in ["valeur", "ymin"]:
                    while method["y1"] > self.x[borne_inf]:
                        borne_inf += 1

                elif method["inf"] == "x":
                    borne_inf = self.indice(method["x1"])
                    if borne_inf == 0:
                        borne_inf = 1
                elif method["inf"] == "y":
                    while method["y1"] > self.y[borne_inf]:
                        borne_inf += 1

                elif method["inf"] == "ymax":
                    borne_inf = variable.argmax()

            else:
                borne_inf = 1

            if "sup" in method:
                if method["sup"] == "chute":
                    borne_sup = variable.argmax()

                    max = variable.argmax()
                    while variable[borne_sup] > max * method["percent"]:
                        borne_sup += 1

                elif method["sup"] == "x":
                    borne_sup = self.indice(method["x2"]) - 1

                elif method["sup"] == "y":
                    borne_sup = borne_inf
                    while method["y2"] > self.y[borne_inf]:
                        borne_sup += 1

                elif method["sup"] == "valeur":
                    borne_sup = borne_inf + 1
                    while self.x[borne_sup] < method["xmax"]:
                        borne_sup += 1

                elif method["sup"] == "ymax":
                    borne_sup = variable.argmax() + 1

            else:
                borne_sup = self.y.size - 1

        if "courbe" not in method:
            borne_inf -= method["keep_inf"]
            borne_sup += method["keep_sup"]

            if borne_inf < 1:
                borne_inf = 1
            if borne_sup > self.x.size - 2:
                borne_sup = self.x.size - 2

        p_1 = np.zeros(borne_inf - 1)
        p_2 = np.zeros(borne_sup - borne_inf + 1) - 1
        p_3 = np.zeros(self.x.size - borne_sup)

        mask = np.concatenate((p_1, p_2, p_3))
        if retour:
            C_temp = Courbe()
            C_temp.x = self.x.compress(mask)
            C_temp.y = self.y.compress(mask)
            return C_temp

        else:
            self.x = self.x.compress(mask)
            self.y = self.y.compress(mask)
            self.borne_cut_indice = [borne_inf, borne_sup]

    # ------------------------------------------------------------------------------

    # def differentiation(self,p):
    #     """
    #     Calcul de la dérivé
    #     """
    #     CT_temp = Courbe()
    #     CT_temp.x = copy.deepcopy(self.x)
    #     CT_temp.y = scipy.signal.savgol_filter(self.y, p, 2,deriv=1,delta = self.x[2]-self.x[1])

    #     return CT_temp

    # -----------------------------------------------------------------------------

    def extrapol(self, a, b, xmin, xmax):
        """
        prolonge une droite jusqu'aux abscisses xmin et xmax choisis
        entrée : indices a et b (svt 0 et -1 pour caluler le coeff directeur de la droite) ; i est le nombre de fois auquel ds notre programme on fait appel à pente_corde avant d'utiliser droite
        """

        pas_moy = (self.x[-1] - self.x[0]) / len(self.x)  # calcul du pas moyen
        nb = int((self.x[0] - pas_moy) / pas_moy)

        curb = copy.deepcopy(self)

        if xmin != None and xmax != None:

            Ldebut = np.linspace(xmin, self.x[0] - pas_moy, num=nb).tolist()
            #            print(Ldebut)
            Lfin = np.linspace(self.x[0] + pas_moy, xmax, num=nb).tolist()

            Liste_X = np.array(Ldebut + self.x.tolist() + Lfin)
        #

        elif xmin != None:

            Ldebut = np.linspace(xmin, self.x[0] - pas_moy, num=nb).tolist()
            #            print('Ldebut',Ldebut)
            Liste_X = np.array(Ldebut + self.x.tolist())

        #
        #

        elif xmax != None:

            Lfin = np.linspace(self.x[0] + pas_moy, xmax, num=nb).tolist()
            Liste_X = np.array(self.x.tolist() + Lfin)

        n = (self.y[b] - self.y[a]) / (self.x[b] - self.x[a])

        print("n", n)

        self.offsety = self.y[a] - self.x[a] * n

        self.y = np.array([Liste_X[k] * n + self.offsety for k in range(len(Liste_X))])

        self.x = np.array(Liste_X)

    #        plt.figure()

    # ------------------------------------------------------------------------------

    #     def extract_part(self,indices   ):  #i1,i2
    #         """
    #         Methode d'extraction d'une sous partie de la courbe à partir d'une liste d'indices
    #         ===> return la sous partie de la courbe
    #         """
    # #        assert type(indices) == 'list'
    #         i1 = indices[0]
    #         i2 = indices[1]
    #         indiceslin = np.linspace(i1,i2,i2-i1+1)
    #         indices2 = indiceslin.tolist()

    #         C_temp = copy.deepcopy(self)
    # #        X = copy.deepcopy()
    # #        C_temp.plot()
    #         C_temp.x = np.take(C_temp.x,indices2)
    #         C_temp.y = np.take(C_temp.y,indices2)
    # #        C_temp.x = np.take(self.x,indices)
    # #        C_temp.y = np.take(self.y,indices)
    # #        C_temp.plot()

    #         return C_temp

    # ------------------------------------------------------------------------------------------------------------
    def indice(self, X):
        """
        Gets and index from a x value

        :param X: scallar
        :type profile: float

        :return: the first superior index to the X given
        :rtype: int

        """
        for indice in range(self.size):
            if self.x[indice] >= X:
                return indice

        print("no X with this index")

    # ------------------------------------------------------------------------------------------------------------
    #     def intersect(self, courbe, methode = 'difference'):
    #         """
    #         retourne l'indice de position de l'intersetion des deux courbes pour les deux courbes
    #         ayant la même abscisse
    #         """

    #         if methode == 'simple' :
    #             self.listey = self.y.tolist()
    #             self.listecourbe = courbe.y.tolist()
    #             H = [val for val in self.listey if val in self.listecourbe  ]
    #             self.intersectiony = self.listey.index(H[0])

    #             return self.intersectiony

    # # Rmq : si qu'une seule intersection et valeur arrondies
    #         elif methode == 'arrondie':
    #             """
    #             retourne l'indice de position de l'intersetion des deux courbes pour les deux courbes
    #             ayant la même abscisse
    #             """
    #             self.listey = self.y.tolist()
    #             listey = [round(k,3) for k in self.listey]
    #             self.listecourbe = courbe.y.tolist()
    #             listecourbe = [round(j,3) for j in self.listecourbe]
    #             H = [val for val in listey if val in listecourbe  ]
    #             J= [listey.index(i) for i in H if listey.index(i)>np.argmax(self.y)]
    #             self.intersectiony= J[0]

    #             return self.intersectiony

    #         elif methode == 'difference' :
    #             assert self.x.all() == courbe.x.all()

    #             C_temp = Courbe()
    #             C_temp.x = copy.deepcopy(courbe.x)
    #             C = courbe.__mul__(-1)
    #             C_temp.y = np.abs(self.__add__(C).y)

    #             intersect = np.argmin(C_temp.y)
    #             print(C_temp.y[intersect])
    # #            print('intersec bon = ', intersect)
    # #            print('valeur intersection =', self.y[intersect])

    #             return intersect, self.y[intersect]

    #         elif methode == 'chgt de signe' :

    #             C_temp = Courbe()
    #             C_temp.x = copy.deepcopy(courbe.x)
    #             C = courbe.__mul__(-1)
    #             C_temp.y = self.__add__(C).y
    #             signe = np.sign(C_temp.y)
    # #            print('signe',signe)
    #             chgt = ((np.roll(signe, 1) - signe) ).astype(int)
    # #            print('chgt',chgt)
    #             indices = np.nonzero(chgt)
    # #            print('indices',indices)

    #             return indices[0].tolist(), [self.y[i] for i in indices[0].tolist()]

    # ------------------------------------------------------------------------------
    def generate_x(self):
        """
        Generate fake x (0 to number of points in y)
        """
        self.x = np.linspace(0, self.y.size - 1, self.y.size)

    # ------------------------------------------------------------------------------
    def load_txt(self, filename, *args, **kwargs):
        """
        load the x and y of a curve from a two columns text file
        """
        filename = otls.replace_separator(filename)
        self.x, self.y = np.loadtxt(
            filename, unpack=True, *args, **kwargs
        )

    # ------------------------------------------------------------------------------
    def load(self, filename, *args, **kwargs):
        """
        load the x and y of a curve from a two columns npy file
        """
        filename = otls.replace_separator(filename)
        A = np.load(file= filename )
        self.x, self.y = A[:,0],A[:,1]

    # ------------------------------------------------------------------------------
    @property
    def max(self):
        """ Return the maximum of the Courbe instance """
        return max(self.y)

    # ------------------------------------------------------------------------------
    @property
    def min(self):
        """ Return the minimum of the Courbe instance  """
        return min(self.y)

    # ------------------------------------------------------------------------------
    @property
    def x_of_max(self):
        """
        :return: Return the x corresponding to ymax of the Courbe instance
        :rtype: numpy array
        """
        #        np.where(self.y==self.max)
        return self.x[np.where(self.y == self.max)]

    @property
    def index_of_max(self):
        """
        :return: Return the indexes corresponding to ymax of the Courbe instance
        :rtype: numpy array """
        return np.where(self.y == self.max)

    @property
    def x_of_min(self):
        """
        :return: Return the x corresponding to ymin of the Courbe instance
        :rtype: numpy array
        """
        #        np.where(self.y==self.max)
        return self.x[np.where(self.y == self.min)]

    @property
    def index_of_min(self):
        """
        :return: Return the indexes corresponding to ymin of the Courbe instance
        :rtype: numpy array """
        return np.where(self.y == self.min)

    @property
    def dx(self):
        """
        Return the x step (dx) of the Courbe instance
        """
        if (self.x[0:-2] - self.x[1:-1] < 1e-6).all():
            return self.x[1] - self.x[0]
        else:
            raise ValueError("dx not is not constant")

    # ------------------------------------------------------------------------------

    def offset_y(self, offset=None):
        """
        Offset the y attribute of the Courbe() instance, the Courbe instance is modified

        :param offset: the amont of the offset, if not given the courbe is offseted so that self.y[0] = 0
        :type offset: float
        """
        if offset == None:
            offset = -self.y[0]
        self.y = self.y + offset

    # ------------------------------------------------------------------------------
    def offset_x(self, offset=None, method="simple"):
        """
        Offset the x attribute of the Courbe() instance, the Courbe instance is modified

        :param offset: the amont of the offset, if not given the courbe is offseted so that self.x[0] = 0
        :type offset: float

        :param method: two methods are available *'simple'* so the self.x is modified or *'interopation'*
        so the self.x is not modified but the self.y is modified accordingly
        :type method: string

        """
        if offset == None:
            self.x = self.x - self.x[0]
        elif offset == 0.0:
            pass
        else:
            if method == "simple":
                self.x = self.x + offset

            if method == "interpolation":
                new_x = self.x + offset
                f = interp1d(new_x, self.y, bounds_error=False, fill_value=0.0)
                new_y = f(self.x)
                self.y = new_y

    # ------------------------------------------------------------------------------
    def regression(self, x1=None, x2=None, method="regression"):
        """
        Compute the slope of the Courbe between two x values

        :param x1: lower limit of the interval of interest
        :type x1: float

        :param x2:upper limit of the interval of interest
        :type x2: float

        :param method: method
        :type method: string

        Calcule la pente de la courbe dans un intervale suivant x par une méthode de corde
        les résultats sont stocké dans des listes, ce qui pemet d'appeller plusieurs fois la fonction
        axe : x (default) ou y
        method : corde (default) ou regression
        return les coefficient de la regression

        :param offset: the amont of the offset, if not given the courbe is offseted so that self.x[0] = 0
        :type offset: float


        """

        if x1 == None:
            x1 = self.x[0]
        if x2 == None:
            x2 = self.x[-1]

        if method == "corde":
            f = interp1d(self.x, self.y)
            pente = (f(x2) - f(x1)) / (x2 - x1)

        elif method == "regression":

            method = {"inf": "x", "sup": "x", "x1": x1, "x2": x2}
            options = {"retour": True}
            C_temp = self.cut(method=method, options=options)
            pente, intercept, r_value, p_value, std_err = stats.linregress(
                C_temp.x, C_temp.y
            )
        return pente

    # ------------------------------------------------------------------------------

    def plot(self, option="b", xlog=False, ylog=False, *args, **kwargs):
        """
        Call matplotlib to generate a plot of the values

        :param option: linetyle of the curve in matplotlib.pyplot format
        :type option: string

        :param xlog: True to plot with the x axis logarithmic (default = False)
        :type xlog: boolean

        :param ylog: True to plot with the y axis logarithmic (default = False)
        :type ylog: boolean

        .. note:: Other paramters understandable by matplotlib can be passed

        """
        plt.figure()
        if hasattr(self, "label"):
            plt.plot(self.x, self.y, option, label=self.label, *args, **kwargs)
        else:
            plt.plot(self.x, self.y)
        if self.xlabel != None:
            plt.xlabel(self.xlabel, fontsize=16)
        if self.ylabel != None:
            plt.ylabel(self.ylabel, fontsize=16)
        if hasattr(self, "title"):
            plt.title(self.title, fontsize=16)

        if self.show_approx == True:
            for i in range(len(self._droite)):
                plt.plot(self.x, self._droite[i], label="approx")
                plt.plot(self._points_x[i], self._points_y[i], "or", label="bornes")
        if xlog == True:
            plt.xscale("log")
        if ylog == True:
            plt.yscale("log")

        plt.legend()
        plt.show()

    # ------------------------------------------------------------------------------
    def project(self, A, kind="linear", fill_value=0):
        """
        Modify a Courbe instance so its x will match the one of a given Courbe instance

        :param A: Courbe instance with the reference x
        :type A: Courbe instance

        :param kind: kind of interpolation available default is linear for other options see
        scipy.interpolate.interp1d
        :type kind: string

        :param fill_value: value to put for x' outside of the original Courbe ones (default = 0.0)
        :type fill_value: float

        :return: new Courbe instance with A's x
        :rtype: Courbe instance
        """
        Ctemp = copy.deepcopy(self)

        f = interp1d(
            self.x, self.y, kind=kind, bounds_error=False, fill_value=fill_value
        )
        if isinstance(A, Courbe):
            Ctemp.y = f(A.x)
            Ctemp.x = A.x
        elif isinstance(A, np.ndarray):
            Ctemp.y = f(A)
            Ctemp.x = A

        return Ctemp

    # ------------------------------------------------------------------------------
    def resample(self, N):
        """
        Resample a Courbe instance with a ratio of N, modify the Courbe instance

        :param N: resampling ratio
        :type N: integer
        """
        self.x = self.x[::N]
        self.y = self.y[::N]

    # ------------------------------------------------------------------------------
    def save_txt(self, filename):
        """
        Save the Courbe instance x and y values in a two columns text file, typically a .asc
        """
        filename = otls.replace_separator(filename)
        directory = os.path.dirname(filename)

        A = np.vstack((self.x, self.y))
        B = np.transpose(A)

        np.savetxt(filename, B)

    # ------------------------------------------------------------------------------

    def save(self, filename):
        """
        Save the Courbe instance x and y values in a two columns .npy file
        """
        filename = otls.replace_separator(filename)
        directory = os.path.dirname(filename)

        A = np.vstack((self.x, self.y))
        B = np.transpose(A)

        np.save(filename, B)

    # ------------------------------------------------------------------------------
    @property
    def size(self):
        """
        Return the size of the Courbe instance x attribute
        """
        return self.x.size

    # ------------------------------------------------------------------------------------------------------------
    def suprNeg(self):
        """
        Modify the negative value of the y in the Courbe instance to 0.0
        """

        for i in range(len(self)):
            if self.y[i] < 0.0:
                self.x[i] = 0.0
                self.y[i] = 0.0

    # -----------------------------------------------------------------------------------------------------------
    def y_from_x(self, x, *args, **kwargs):
        """
        return the y values corresponding to a x value using linear interpolation

        currentlly only work if x is monotonous
        """
        np.all(np.diff(self.x) > 0)
        f = interp1d(self.x, self.y, *args, **kwargs)
        return float(f(x))

    # -----------------------------------------------------------------------------------------------------------
    def x_from_y(self, y, *args, **kwargs):
        """
        return the x values corresponding to a y value using linear interpolation

        currentlly only work if y is monotonous
        """

        np.all(np.diff(self.y) > 0)
        f = interp1d(self.y, self.x, *args, **kwargs)
        return float(f(y))

    # ------------------------------------------------------------------------------------------------------------
    def __add__(self, A):
        """
        Adding method, currently supported with interger, float and other curve with an identical self.x
        """
        if isinstance(A, int) or isinstance(A, float):
            C_temp = Courbe()
            C_temp.x = copy.deepcopy(self.x)
            C_temp.y = copy.deepcopy(self.y)
            C_temp.y = C_temp.y + A
            return C_temp

        elif isinstance(A, Courbe):
            if np.array_equal(self.x, A.x):

                #            self.x.all() == .all():
                C_temp = Courbe()
                C_temp.x = copy.deepcopy(self.x)
                C_temp.y = copy.deepcopy(self.y) + copy.deepcopy(A.y)
                return C_temp
            else:
                raise Exception("Addition between curves with different x, is not supported") 
                
                # C_temp = Courbe()
                # X = np.append(self.x, A.x)
                # X = np.unique(X)
                # X = np.sort(X)
                # new_self = self.project(X)
                # new_A = A.project(X)
                # C_temp.x = X
                # C_temp.y = new_self.y + new_A.y
                # return C_temp

        else:
            print("l'addition avec un autre type de donnée n'est pas supporté")

    # -----------------------------------------------------------------------------------------------------------
    def __radd__(self, other):
        return self.__add__(other)

    # -----------------------------------------------------------------------------------------------------------
    def __div__(self, A):
        """
        Methode de division entre courbe
        supporté pour la division avec un entier, un float et un autre courbe ayant le même self.x
        """
        if isinstance(A, int) or isinstance(A, float):
            CT_temp = copy.deepcopy(self)
            CT_temp.y = copy.deepcopy(self.y)
            CT_temp.y = CT_temp.y / A
            return CT_temp

        elif isinstance(A, Courbe):
            if np.array_equal(self.x, A.x):
                C_temp = copy.deepcopy(self)
                C_temp.y = copy.deepcopy(self.y) / copy.deepcopy(A.y)
                return C_temp
            else:
                print(
                    "la division entre deux courbe n'ayant pas le même pas de temps n'est pas supporté"
                )
        else:
            print("la division avec un autre type de donnée n'est pas supporté")

    # -----------------------------------------------------------------------------------------------------------
    def __truediv__(self, a):
        return self.__div__(a)

    # -----------------------------------------------------------------------------------------------------------
    def __rtruediv__(self, a):
        if isinstance(a, int) or isinstance(a, float):
            CT_temp = copy.deepcopy(self)
            CT_temp.y = copy.deepcopy(self.y)
            CT_temp.y = a / CT_temp.y
            return CT_temp
        elif isinstance(A, Courbe):
            if np.array_equal(self.x, A.x):
                C_temp = copy.deepcopy(self)
                C_temp.y = copy.deepcopy(A.y) / copy.deepcopy(self.y)
                return CT_temp
            else:
                print(
                    "la division entre deux courbe n'ayant pas le même pas de temps n'est pas supporté"
                )
        else:
            print("la division avec un autre type de donnée n'est pas supporté")

    # ------------------------------------------------------------------------------------------------------------
    def __mul__(self, A):
        if isinstance(A, int) or isinstance(A, float):
            CT_temp = copy.deepcopy(self)
            CT_temp.x = copy.copy(self.x)
            CT_temp.y = np.multiply(self.y, A)
            #
            return CT_temp

        elif isinstance(A, Courbe):
            C_temp = copy.deepcopy(self)
            C_temp.y = copy.deepcopy(np.multiply(self.y, A.y))
            return C_temp

        else:
            print("la multiplication avec un autre type de donnée n'est pas supporté")

    # ------------------------------------------------------------------------------------------------------------
    def __rmul__(self, other):
        return self.__mul__(other)

    # ------------------------------------------------------------------------------------------------------------
    def __neg__(self):
        C_temp = copy.deepcopy(self)
        C_temp.y = -copy.deepcopy(self.y)
        return C_temp

    # ------------------------------------------------------------------------------------------------------------
    def __pow__(self, A):
        """
        Get a Courbe instance to a certain power

        :param A: power at which to elevate the Courbe y values
        :type A: integer, float or Courbe instance with same x
        """
        if isinstance(A, int) or isinstance(A, float):
            CT_temp = copy.deepcopy(self)
            CT_temp.y = CT_temp.y ** A
            return CT_temp

        elif isinstance(A, Courbe):
            print("la puissance entre deux courbe n'est pas supporté")
        else:
            print("la puissance avec un autre type de donnée n'est pas supporté")

    def __sub__(self, A):
        """
        Sustraction method, currently supported with interger, float and other Courbe
        """
        if isinstance(A, int) or isinstance(A, float):
            C_temp = Courbe()
            C_temp.x = copy.deepcopy(self.x)
            C_temp.y = copy.deepcopy(self.y)
            C_temp.y = C_temp.y - A
            return C_temp

        elif isinstance(A, Courbe) or isinstance(A, pycurves.courbe.Courbe):
            if np.array_equal(self.x, A.x):
                C_temp = Courbe()
                C_temp.x = copy.deepcopy(self.x)
                C_temp.y = copy.deepcopy(self.y) - copy.deepcopy(A.y)
                return C_temp
            else:
                C_temp = Courbe()
                X = np.append(self.x, A.x)
                X = np.unique(X)
                X = np.sort(X)
                new_self = self.project(X)
                new_A = A.project(X)
                C_temp.x = X
                C_temp.y = new_self.y - new_A.y
                return C_temp

        else:
            print("substration with another data type is not supported")

    # -----------------------------------------------------------------------------------------------------------
    def __rsub__(self, other):
        return -1 * self.__sub__(other)

    # -----------------------------------------------------------------------------------------------------------
    def __str__(self):
        return str(self.x) + "\n" + str(self.y)

    # -----------------------------------------------------------------------------------------------------------
    def __repr__(self):
        """
        .. todo:: This function should be improved, however it seems to be rarely used
        """
        return str(self.x) + "\n" + str(self.y)

    # -----------------------------------------------------------------------------------------------------------
    def __eq__(self, A):
        return self.y.__eq__(A)
