# -*- coding: utf-8 -*-
"""
Module to provide a container for several tests results that are linked, as tests for the same campaign done at different loading rate for instance.

"""

from pycurves.essai import *

import os
import os.path
from os.path import basename

#from misc.os_tools import *


class Campagne(dict):
    """
    Container for several Essai() object, inherit from dict
    """

    def __init__(self,):
        dict.__init__(self)

    def resize_courbes(self, nb_pt_cible=2000):
        """
        Resize all curves to get close to a number of point defined.

        """
        for essai in self:
            rapport = self[essai][self[essai].courbes()[0]].size / nb_pt_cible

            for com in self[essai].courbes():
                self[essai][com].resample(rapport)

    def generate_labels(self):
        """
        Automatically generates the label of the curves using the Campagne keys and the **pycurves.essai.Essai() instance keys.
        """
        for essai in self:
            for courbe in self[essai].courbes():
                self[essai][courbe].label = essai + "_" + courbe
