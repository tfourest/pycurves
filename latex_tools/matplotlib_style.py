# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt


def init_latex_graph(*args, **kwargs):

    pic_quality = "adawi_2"  # adawi_1 (fig width = textwidth)
    # adawi_2 (fig width set to allow for two
    #          figs side by side)
    # adawi_3 (fig size adapted to fit 3 figs
    #          above each other)

    pic_format = ["pdf"]  # format of pictures: png, pdf (list)
    Latex = True

    plt.close("all")

    # -- Read general plot style ---------------------------------------------------
    plot_style = start(qlty=pic_quality, frmt=pic_format, tex=Latex)

    for key in plot_style.keys():
        if key != "rcParams":
            exec(key + " = plot_style['" + key + "']")

        elif key == "rcParams":
            for key_rcP in plot_style["rcParams"].keys():
                plt.rcParams[key_rcP.replace("_", ".", 1)] = plot_style["rcParams"][
                    key_rcP
                ]

    plt.rcParams["text.latex.preamble"] = r"\usepackage{lmodern}"

    plt.rcParams["axes.prop_cycle"] = (
        "cycler('color', ['k' , 'b', 'g', 'C44E52'])"
        + "+ cycler('linestyle', ['-', '--', '-.', '.']) "
        + "+ cycler('marker', ['o', 's', '^','D']) "
    )


def start(
    qlty="adawi_1", dpi=300, frmt=["jpg"], tex=True, xkcd=False, conf="", jour=""
):
    """
    Common definition of matplotlib plot styles used for
    
    qlty = picture quality refers to the use of the figure, for instance:
    
              # adawi_1 (fig width = textwidth)
              # adawi_2 (fig width = 0.49 textwidth to allow for two figs
                         side by side)
              # adawi_3 (fig size adapted to fit 3 figs above each other)
    
    frmt = defines a LIST of formats to be created (dflt is jpg)
           (when creating PDFs it may be nice to also create JPGs or PNGs for
            the ease of use)
    
    tex  = uses LaTex font when True (dflt)
    
    xkcd = uses XKCD style for lines etc. when True (dflt is False)
    """

    # -- Plot style ---------------------------------------------------------------
    if xkcd == True:
        tex = False  # reset LaTeX font if xkcd is ON
        plt.xkcd()

    # -- Font style (latex font or arial like) ------------------------------------
    if tex == False:  # use standard font
        font_style1 = ""
        font_style2 = ""

        plt.rcParams["text.usetex"] = False

    elif tex == True and xkcd != True:  # use Latex font (system requirements)
        font_style1 = r"$\displaystyle"
        font_style2 = r"$"

        plt.rcParams["text.latex.preamble"] = [r"\usepackage{lmodern}"]
        plt.rcParams["text.latex.preamble"] = [r"\usepackage[T1]{fontenc}"]
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["font.serif"] = "Computer Modern Roman"

        plt.rcParams["text.usetex"] = True

    # -- Picture settings ---------------------------------------------------------

    # dummy parameters (maybe later they are overwritten)
    axis_size = []

    if "adawi" in qlty:
        if qlty == "adawi_1":
            figsize_x = 179.0 / 25.4
            figsize_y = (figsize_x + 10.0 / 25.4) * 4.0 / 6.0
        elif qlty == "adawi_2":
            figsize_x = 89.0 / 25.4
            figsize_y = (figsize_x + 10.0 / 25.4) * 4.0 / 6.0
        elif qlty == "adawi_3":
            figsize_x0 = 0.38 * (179.0 + 89.0) / 25.4
            figsize_y = (figsize_x0 + 10.0 / 25.4) * 4.0 / 6.0
            figsize_x = 139.0 / 25.4

        font_size1 = 18.0
        font_size2 = font_size1 - 2
        font_size3 = font_size2 - 2

        pic_dpi = 300

        xtick_labelsize = font_size2  # font_size
        ytick_labelsize = font_size2  # font_size
        axes_linewidth = 0.5
        axes_labelsize = font_size1

        axes_grid = True
        grid_alpha = 0.5
        grid_color = "black"
        grid_linestyle = "dotted"
        grid_linewidth = 0.5

        legend_fontsize = font_size3
        legend_frameon = True
        legend_fancybox = True
        legend_numpoints = 1
        legend_borderpad = 0.25
        legend_frame_set_linewidth = 0.5
        legend_frame_set_alpha = 0.7

        lines_linestyle = "-"
        lines_linewidth = 1.25
        lines_markersize = 5
        lines_markerfacecolor = "w"
        lines_markeredgewidth = 1.0

        lines_alpha = 1.0

    # -- Assemble return parameters --------------------------------------------
    plot_style = {
        "qlty": qlty,
        "pic_dpi": pic_dpi,
        "frmt": frmt,
        "tex": tex,
        "xkcd": xkcd,
        "conf": conf,
        "jour": jour,
        "font_style1": font_size1,
        "font_style2": font_size2,
        "axis_size": axis_size,
        "figsize_x": figsize_x,
        "figsize_y": figsize_y,
        "font_size": font_size1,
        "lines_alpha": lines_alpha,
        "legend_frame_set_linewidth": legend_frame_set_linewidth,
        "legend_frame_set_alpha": legend_frame_set_alpha,
        "rcParams": {
            "xtick.labelsize": xtick_labelsize,
            "ytick.labelsize": ytick_labelsize,
            "axes.linewidth": axes_linewidth,
            "axes.labelsize": axes_labelsize,
            "axes.grid": axes_grid,
            "grid.alpha": grid_alpha,
            "grid.color": grid_color,
            "grid.linestyle": grid_linestyle,
            "grid.linewidth": grid_linewidth,
            "legend.fontsize": legend_fontsize,
            "legend.frameon": legend_frameon,
            "legend.fancybox": legend_fancybox,
            "legend.numpoints": legend_numpoints,
            "legend.borderpad": legend_borderpad,
            "lines.linestyle": lines_linestyle,
            "lines.linewidth": lines_linewidth,
            "lines.markersize": lines_markersize,
            "lines.markerfacecolor": lines_markerfacecolor,
            "lines.markeredgewidth": lines_markeredgewidth,
            "font.size": font_size1,
        },
    }

    return plot_style
