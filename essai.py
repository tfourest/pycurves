# -*- coding: utf-8 -*-
"""
Module to provide a container for several curves that are linked, as curves from the same test.
"""

import matplotlib.pyplot as plt
import numpy as np

import pycurves as pyc
import copy


class Essai(dict):
    """
    Object to simplify dealing with numerous curves that are linked together
    in the sense that they comes from the same test or from the same simulation.

    This class inherits from dictionnary and works in the same way.

    For the creator of pycurves.essai.Essai instance :

        :param fichier: Name of the file to load. It is optional, the instance can be created empty.
        It can read '.txt' files and '.csv' files
        :type fichier: string

        :param formatf: key to determine the format to open. The possibilities are :
            None in which case the delimiter and skiprows parameters of numpy.loadtxt should be passed to the function.
            ; 'dewetron' :  delimiter = '\t'
            ; 'dewetron space' : delimiter = ' ' and skiprows = 11
            ; 'dewesoft' : delimiter = ',' and skiprows = 10
        that are to be plotted. By default all **pycurves.courbe.Courbe** are plotted.
        :type formatf: list

    """

    def __init__(self, fichier=None, formatf="dewetron", *args, **kwargs):

        dict.__init__(self)

        if fichier != None:
            if (
                fichier[-3] + fichier[-2] + fichier[-1] == "txt"
                and formatf == "dewetron"
            ):
                tableau = np.loadtxt(
                    fichier,
                    skiprows=11,
                    dtype=np.unicode_,
                    delimiter="\t",
                    ndmin=2,
                    unpack=True,
                )
            elif (
                fichier[-3] + fichier[-2] + fichier[-1] == "txt"
                and formatf == "dewesoft space"
            ):
                tableau = np.loadtxt(
                    fichier,
                    dtype=np.unicode_,
                    delimiter=" ",
                    ndmin=2,
                    unpack=True,
                    skiprows=11,
                )
            elif (
                fichier[-3] + fichier[-2] + fichier[-1] == "txt"
                and formatf == "dewesoft"
            ):
                tableau = np.loadtxt(
                    fichier,
                    dtype=np.unicode_,
                    delimiter=",",
                    ndmin=2,
                    unpack=True,
                    skiprows=10,
                )

            elif fichier[-3] + fichier[-2] + fichier[-1] == "csv":
                tableau = np.loadtxt(
                    fichier, dtype=np.unicode_, delimiter=",", ndmin=2, unpack=True
                )
            else:
                tableau = np.loadtxt(
                    fichier, dtype=np.unicode_, ndmin=2, unpack=True, *args, **kwargs
                )

            tableau2 = np.zeros(tableau.shape)

            for col in range(tableau.shape[0])[0::]:
                for line in range(tableau.shape[1])[1::]:
                    tableau[col, line] = pyc.virgule_dot(tableau[col, line])
                    tableau2[col, line] = float(tableau[col, line])

            for col in range(tableau.shape[0])[1::]:
                self[tableau[col, 0]] = pyc.Courbe()
                self[tableau[col, 0]].x = tableau2[0, 1 : tableau.shape[1]]
                self[tableau[col, 0]].y = tableau2[col, 1 : tableau.shape[1]]
                self[tableau[col, 0]].label = tableau[col, 0]

    def courbes(self):
        """
        iterator on the **pycurves.courbe.Courbe** intances contained
        in the **pycurves.essai.Essai** instance
        """
        liste = []
        for com in self:
            if isinstance(self[com], pyc.Courbe):
                liste.append(com)
        return liste

    def plot(self, c_names=None):
        """
        method to plot all  **pycurves.courbe.Courbe** intances contained in
        the **pycurves.essai.Essai** on several graphs stacked.

        :param c_names: list of the name of the keys corresponding to the
        **pycurves.courbe.Courbe** intances in the **pycurves.essai.Essai** instance
        that are to be plotted. By default all **pycurves.courbe.Courbe** are plotted.
        :type c_names: list
        """
        if c_names == None:
            nb_c = len(self.courbes())
            c_names = self.courbes()
        else:
            nb_c = len(c_names)

        i = 1
        self.plot_obj = plt.figure()
        l = []
        for com in c_names:
            l.append(self.plot_obj.add_subplot(nb_c, 1, i))
            if hasattr(self[com], "label"):
                l[-1].plot(self[com].x, self[com].y, label=self[com].label)
                l[-1].set_title(label=self[com].label)

            else:
                plt.plot(self[com].x, self[com].y)
            i += 1

        self.plot_obj.tight_layout()
        self.plot_obj.show()

    def cut(self, C, method):
        """
        Method to cut all **pycurves.courbe.Courbe** intances in the **pycurves.essai.Essai**
        instance at once.

        :param courbe: the key corresponding to the reference **pycurves.courbe.Courbe** instance in the **pycurves.essai.Essai** instance.
        :type courbe: string

        :param method: same than for the **pycurves.courbe.Courbe.cut** method.
        :type method: dictionnary
        """

        if method == None:
            method = {}

        if "courbe" in method:
            del method["courbe"]
        self[C].cut(method=method)

        method["courbe"] = self[C]

        for com in self:
            if isinstance(self[com], pyc.Courbe) and com != C:
                self[com].cut(method=method)

    def resample(self, N):
        """
        Resample all **pycurves.courbe.Courbe** intances contained
        in the **pycurves.essai.Essai** instance at once.

        :param N: resampling factor.
        :type N: integer
        """
        for com in self:
            if isinstance(self[com], pyc.Courbe):
                self[com].resample(N)


class Test(Essai):
    """Same class than Essai """

    def __init__(self, fichier=None, formatf="dewetron"):
        Essai.__init__(self, fichier=fichier, formatf=formatf)
