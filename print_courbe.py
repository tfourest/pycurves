# -*- coding: utf-8 -*-
"""
This module provide a simple way to visualise pycurves.Courbe instance while processing the data.
"""

from pycurves.courbe import *
import numpy as np

from pycurves.math_courbe import *
import matplotlib.pyplot as plt


import copy


def virgule_dot(word):
    new_word = word
    for i in range(0, len(word)):
        #        if ord(word[i]) != 32:
        if word[i] == ",":
            new_word = new_word.replace(word[i], ".")
    return new_word


def plotMeanMaxMin(lCourbe, show=True, save=False, Name=None):
    mean = Cmean(lCourbe)
    maximum = Cmax(lCourbe)
    minimum = Cmin(lCourbe)
    mean.label = "mean"
    maximum.label = "maximum"
    minimum.label = "minimum"

    fig, ax = plt.subplots()
    ax.fill(maximum.x, maximum.y, color="grey", alpha=0.3)
    ax.fill(
        [maximum.x[0], maximum.x[0], maximum.x[-1], maximum.x[-1]],
        [0, maximum.y[0], maximum.y[-1], 0],
        color="grey",
        alpha=0.3,
    )
    ax.plot(maximum.x, maximum.y, color="k", dashes=[30, 5, 10, 5])

    ax.plot(mean.x, mean.y, "k", label="Mean")

    #        ax.fill(mean.x,mean.y, 'b')
    ax.fill(minimum.x, minimum.y, color="white")
    ax.fill(
        [minimum.x[0], minimum.x[0], minimum.x[-1], minimum.x[-1]],
        [0, minimum.y[0], minimum.y[-1], 0],
        color="white",
    )
    ax.plot(minimum.x, minimum.y, color="k", dashes=[30, 5, 10, 5])
    #        , x, y2, 'r', alpha=0.3)
    plt.legend()
    if show == True:
        plt.show()
    if save == True:
        if Name == None:
            print("veuillez donner un nom à la courbe")
        else:
            plt.savefig(Name + ".pdf", format="pdf")
            plt.clf()


#        plotXcourbes([mean,maximum,minimum])

# ------------------------------------------------------------------------------


def plotXcourbes(lCourbe, loption=None, bounds={}, xlabel=None, ylabel=None):
    """
    Funtion to plot a matplotlib figure from a liste of curves object.
    It calls gen_plt_X_courbes to generate the matplotlib figure object.


    All inputs of **gen_plt_X_courbes** can be use, with the adition of the following:

    :param xlabel: label of the x axis of the plot
    :type xlabel: string

    :param ylabel: label of the y axis of the plot
    :type ylabel: string

    """

    fig, ax = gen_plt_X_courbes(lCourbe, loption0=loption, bounds=bounds)
    fig.canvas.draw()
    ax.legend(loc="best")
    if xlabel != None:
        plt.xlabel(xlabel)
    if ylabel != None:
        plt.ylabel(ylabel)


#    if save == True :
#        if Name == None :
#            print("veuillez donner un nom à la courbe")
#        else :
#            plt.savefig(Name+".pdf",format='pdf')
#            plt.clf()


def gen_plt_X_courbes(lCourbe, loption0=None, bounds={}):
    """
    Funtion to generate a figure object from a liste of curves.

    :param lCourbe:  list of curves to plot (pycurves.Courbe instance)
    :type lCourbe: list

    :param loption0:  list of style to be use in matplotlib.pyplot. If no list is given, ['b', 'g', 'r', 'c', 'm', 'y', 'k'] will be used.
        If there are more curves in **lCourbe** than styles in **loptions0**, the function will loop over the list of style.
    :type loption0: boulean

    :param bounds: dictionnary with boundaries to use for the plot. The keys expected are 'xmin', 'xmax', 'ymin' and 'ymax'.
        By default the boundary are automatically defined by matplotlib.pyplot.
    :type bounds: dict

    :return: the plot object
    :rtype: matplotlib.pyplot instance

    .. Note:: The label attribute of the pycurves.Courbe instance are directly use.

    .. Warning:: It is better to not use this function directly by to use **plotXcourbes** or **plotXcourbesN**.

    """

    if loption0 == None:
        loption0 = ["b", "g", "r", "c", "m", "y", "k"]
    loption = copy.deepcopy(loption0)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for C in lCourbe:
        if hasattr(C, "label"):
            if hasattr(C, "style"):
                ax.plot(C.x, C.y, C.style, label=C.label)
            else:
                ax.plot(C.x, C.y, loption[0], label=C.label)
        else:
            if hasattr(C, "style"):
                ax.plot(C.x, C.y, C.style)
            else:
                ax.plot(C.x, C.y, loption[0])
        del loption[0]
        if loption == []:
            loption = loption

    ax.legend(loc=4)

    xmind, xmaxd = ax.get_xlim()  # return the current xlim
    ymind, ymaxd = ax.get_ylim()  # return the current xlim

    if "xmin" in bounds:
        xmind = bounds["xmin"]
    if "xmax" in bounds:
        xmaxd = bounds["xmax"]
    if "ymin" in bounds:
        ymind = bounds["ymin"]
    if "ymax" in bounds:
        ymaxd = bounds["ymax"]

    ax.set_xlim(xmind, xmaxd)
    ax.set_ylim(ymind, ymaxd)

    loption = copy.deepcopy(loption0)

    return fig, ax


def plotXcourbesN(lcourbe, *arg, **kwargs):
    """
    Same as **plotXcourbes** but with all curves normalised at +1 or -1
    depending of the absolute maximum.
    All inputs are the same than **plotXcourbes**.

    """
    new_lcourbe = []
    for courbe in lcourbe:
        a = copy.deepcopy(courbe)
        a.y = a.y / np.max([np.abs(a.max), np.abs(a.min)])
        new_lcourbe.append(a)
    plotXcourbes(new_lcourbe, *arg, **kwargs)


###############################################################################


#%%
def plot_adawi(x1, y1, x2, y2, xlabel, ylabel, name, tex_dir, lim=[]):

    import pycurves.latex_tools.plot_styles_adawi as plot_styles

    xmin, xmax = 0.0, 1.0
    linewidth_ = 1.25
    markeredgewidth_ = 1
    markevery_ = 1

    pic_quality = "adawi_2"  # adawi_1 (fig width = textwidth)
    # adawi_2 (fig width set to allow for two
    #          figs side by side)
    # adawi_3 (fig size adapted to fit 3 figs
    #          above each other)

    pic_format = ["pdf"]  # format of pictures: png, pdf (list)
    Latex = True

    plt.close("all")

    # -- Read general plot style ---------------------------------------------------
    plot_style = plot_styles.start(qlty=pic_quality, frmt=pic_format, tex=Latex)

    for key in plot_style.keys():
        if key != "rcParams":
            exec(key + " = plot_style['" + key + "']")

        elif key == "rcParams":
            for key_rcP in plot_style["rcParams"].keys():
                plt.rcParams[key_rcP.replace("_", ".", 1)] = plot_style["rcParams"][
                    key_rcP
                ]

    plt.rcParams["text.latex.preamble"] = r"\usepackage{lmodern}"

    plt.figure(1)
    plt.plot(
        x1, y1, linestyle="solid", color="k", linewidth=linewidth_, label="Analytic"
    )
    plt.plot(
        x2,
        y2,
        linestyle="solid",
        color="none",
        marker="o",
        markeredgecolor="r",
        markerfacecolor="white",
        linewidth=linewidth_,
        markeredgewidth=markeredgewidth_,
        markersize=8,
        alpha=0.65,
        markevery=(1, markevery_),
        label="Grid method",
    )
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)
    plt.xlim(xmin, xmax)

    if lim:
        plt.ylim(lim[0], lim[1])
    plt.xticks(fontsize=18)  # , rotation=90
    plt.yticks(fontsize=18)  # , rotation=90
    plt.legend(ncol=1, shadow=True, fancybox=True, fontsize=18)
    plt.tight_layout()
    plt.savefig(tex_dir + name, dpi=600)

    plt.show()

    print("Done")
