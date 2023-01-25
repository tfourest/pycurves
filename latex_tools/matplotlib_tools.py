# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from pycurves.latex_tools.matplotlib_style import init_latex_graph

import numpy as np


#%%----------------------------------------------------------------------------
def writeXcourbesMPL(
    lCourbe, directory=None, name=None, bound={}, title=None, xlabel=None, ylabel=None
):
    """
    Fonction to create pdf plot from a list of Courbe() with a LateX style
    
    ---------------------------------------------------------------------------
    :param lcourbe: list of courbe object to plot
    :type lcourbe: list
    
    :param directory: saving directory
    :type directory; str
    
    :param name: saving filename (without the extension)
    :type name; str
            
    :param bound : a dictionnary that can contain values for xmin, xmax, ymin, ymax
    :type bound: dict

    :param title: graphic title
    :type title: str
    
    :param xlabel: x axis label
    :type xlabel: str
    
    :param ylabel: y axis label
    :type ylabel: str
    ---------------------------------------------------------------------------
    :return: pdf figure
    -----------------------------------------------------------------------
    """
    leg = []

    init_latex_graph()

    plt.figure(1)

    nb_mark = 5
    markshift = 0
    nb_c = len(lCourbe)

    for c in lCourbe:

        c_size = len(c.x)
        mk_step = list(
            np.arange(c_size - 1, step=int((c_size / nb_mark)))
            + int(markshift * np.floor(c_size / (nb_mark * nb_c)))
        )

        plt.plot(c.x, c.y, markevery=mk_step)
        leg += [c.label]
        markshift += 1

    if xlabel:
        plt.xlabel(xlabel)

    if ylabel:
        plt.ylabel(ylabel)

    plt.legend(leg)

    plt.title(title)

    plt.tight_layout()

    if directory and name:
        plt.savefig(directory + name + ".pdf")

    plt.show()
