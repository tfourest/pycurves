# -*- coding: utf-8 -*-
"""
This module provide some functions to generate latex (pgfplot) graphical outputs from pycurves.Courbe instance.
"""


from pycurves.courbe import *
import numpy as np

# from math_courbe import *
import matplotlib.pyplot as plt

import os
import subprocess

from pycurves.latex_tools.pgf_styles import *


###############################################################################
def convert_letter(word):
    new_word = ""
    for i in range(0, len(word)):
        #        if ord(word[i]) != 32:
        if word[i] == "_":
            new_word += "\\_"
        else:
            new_word += word[i]
    return new_word


###############################################################################

###############################################################################


###############################################################################


def writeXcourbesPGF(
    lCourbe,
    directory=None,
    name=None,
    style=Tom_standart_BW_curves,
    bound={},
    xlabel=None,
    ylabel=None,
    title=None,
    compil=False,
):
    """
    Fonction to create pgfplot pdf file from a list of pycurves.Courbe(). 
    Writes the .tex and .bat file needed to compile the pgf figure.
    
    :param lCourbe:  list of curves to plot (pycurves.Courbe instance)
    :type lCourbe: list      

    :param directory:  name of the directory where to save the files
    :type directory: string

    :param name:  base name of the files to save
    :type directory: string

    :param style:  the style to use for latex output (pycurves.Courbe instance). Default is Tom_standart_BW_curves
    :type style:  pycurves.latex_tools.pgf_style
    
    :param bounds: dictionnary with boundaries to use for the plot. The keys expected are 'xmin', 'xmax', 'ymin' and 'ymax'. 
        By default the boundary are automatically defined by matplotlib.pyplot.
    :type bounds: dict
         
    :param xlabel: label of the x axis of the plot. Default is None. 
    :type xlabel: string
    
    :param ylabel: label of the y axis of the plot. Default is None.
    :type ylabel: string

    :param title: title of the plot. Default is None.
    :type title: string    

    :param compil: whether the latex file are compiled when generated. Default is True
    :type compil: boolean
            
       
    """

    base_directory = os.getcwd()
    os.chdir(directory)

    if not hasattr(style, "axis_type"):
        style.axis_type = "axis"

    #    fichier = open(directory + name + '.tex','w')
    fichier = open(name + ".tex", "w")

    fichier.write("\\documentclass[12pt]{article}\n")
    fichier.write("\\usepackage{tikz}\n")
    fichier.write("\\usepackage{pgfplots}\n")

    fichier.write("\\pgfrealjobname{" + name + "}\n")

    fichier.write("\\begin{document}\n")

    fichier.write(style.plotcyclelist)
    fichier.write(
        """\pgfplotscreateplotcyclelist{myonlymark}{%
only_marks\\\\%
}
"""
    )

    fichier.write("\\beginpgfgraphicnamed{" + name + "_}\n")

    fichier.write("\\begin{tikzpicture}\n")
    fichier.write("\\begin{" + style.axis_type + "}[\n")

    fichier.write(style.axis)
    if lCourbe[0].label != None:
        fichier.write("legend entries={")
        for C in lCourbe:
            fichier.write(convert_letter(str(C.label)) + ",")
        fichier.write("},\n")

    if xlabel != None:
        fichier.write("xlabel={" + xlabel + "},\n")

    if ylabel != None:
        fichier.write("ylabel={" + ylabel + "},\n")

    if title != None:
        fichier.write("title={" + title + "},\n")

    if bound == None:
        bound = {}
    if "xmin" in bound:
        fichier.write("xmin=" + str(bound["xmin"]) + ",\n")
    if "xmax" in bound:
        fichier.write("xmax=" + str(bound["xmax"]) + ",\n")
    if "ymin" in bound:
        fichier.write("ymin=" + str(bound["ymin"]) + ",\n")
    if "ymax" in bound:
        fichier.write("ymax=" + str(bound["ymax"]) + ",\n")

    fichier.write("]\n")

    j = 0
    for C in lCourbe:
        fichier.write("\\addplot+[thick,smooth,")
        if isinstance(style.addplot, list):
            fichier.write(style.addplot[j])
            j += 1
            if j >= len(style.addplot):
                j = 0
        elif isinstance(style.addplot, str):
            fichier.write(style.addplot)
        if hasattr(style, "automatic_repeat"):
            if style.repeat == False:
                fichier.write("]\n")
            else:
                fichier.write("mark repeat = " + str(int(C.x.size / 10)) + "]\n")
        else:
            fichier.write("mark repeat = " + str(int(C.x.size / 10)) + "]\n")
        fichier.write("coordinates {")

        for i in np.arange(len(C.x), step=C.x.size // 1000 + 1):
            fichier.write("(" + str(C.x[i]) + "    ,    " + str(C.y[i]) + "    )\n")
        fichier.write("};\n")
    #
    #

    fichier.write("\\end{" + style.axis_type + "}")
    fichier.write(
        """\\end{tikzpicture}
\\endpgfgraphicnamed
\\end{document}
 """
    )
    fichier.close()

    #    fichier = open(directory + name + '_compil.bat','w')
    fichier = open(name + "_compil.bat", "w")
    fichier.write("pdflatex " + name + ".tex\n")
    fichier.write("pdflatex --jobname " + name + "_ " + name + ".tex\n")
    fichier.close()
    os.chdir(base_directory)

    if compil == True:

        os.chdir(directory)
        #        os.spawnl(os.P_NOWAIT, name + '_compil.bat')
        #        print("pdflatex", directory +  name  +".tex")
        #        print("pdflatex","--jobname", name +"_ ",   name  +".tex")
        #        os.spawnl(os.P_NOWAIT, name + '_compil.bat')
        subprocess.Popen(name + "_compil.bat")

        #        subprocess.Popen(["pdflatex", name  +".tex"])
        #        subprocess.Popen(["pdflatex","--jobname", name +"_ ", name  +".tex"])

        #        os.system( name + '_compil.bat')
        #        os.chdir('..\\')
        os.chdir(base_directory)


#    os.system("pdflatex", directory +  name  +".tex")
#    nom = name+"_ "


###############################################################################


def writeXcourbesPGF_3yaxis(
    lCourbe,
    directory=None,
    name=None,
    style=[TSBWC, TSBWC, TSBWC],
    bound={},
    xlabel=None,
    ylabel=None,
    compil=False,
):
    """
    Fonction to create pgfplot pdf file from a list of Courbe(). With 3 axes. 
    This is not wery versatile so you should ask before using it.
    """

    base_directory = os.getcwd()
    os.chdir(directory)

    for stl in style:
        if not hasattr(stl, "axis_type"):
            stl.axis_type = "axis"

    #    fichier = open(directory + name + '.tex','w')
    fichier = open(name + ".tex", "w")

    fichier.write("\\documentclass[12pt]{article}\n")
    fichier.write("\\usepackage{tikz}\n")
    fichier.write("\\usepackage{pgfplots}\n")

    fichier.write("\\pgfrealjobname{" + name + "}\n")

    fichier.write("\\begin{document}\n")

    for stl in style:
        fichier.write(stl.plotcyclelist)

    fichier.write(
        """\pgfplotscreateplotcyclelist{myonlymark}{%
only_marks\\\\%
}
"""
    )

    fichier.write("\\beginpgfgraphicnamed{" + name + "_}\n")

    fichier.write("\\begin{tikzpicture}\n")

    j = 0
    for list_of_C in lCourbe:
        fichier.write("\\begin{" + style[j].axis_type + "}[\n")

        fichier.write(style[j].axis)

        if list_of_C[0].label != None:
            fichier.write("legend entries={")
            for C in list_of_C:
                fichier.write(convert_letter(str(C.label)) + ",")
            fichier.write("},\n")

        if xlabel != None:
            fichier.write("xlabel={" + xlabel + "},\n")

        if ylabel != None:
            fichier.write("ylabel={" + ylabel[j] + "},\n")

        if bound == None:
            bound = []
            for i in range(j):
                bound.append({})
        if "xmin" in bound[j]:
            fichier.write("xmin=" + str(bound["xmin"]) + ",\n")
        if "xmax" in bound[j]:
            fichier.write("xmax=" + str(bound["xmax"]) + ",\n")
        if "ymin" in bound[j]:
            fichier.write("ymin=" + str(bound["ymin"]) + ",\n")
        if "ymax" in bound[j]:
            fichier.write("ymax=" + str(bound["ymax"]) + ",\n")
        fichier.write("]\n")

        if hasattr(style[j], "plotset"):
            fichier.write(style[j].plotset)

        for C in list_of_C:
            fichier.write("\\addplot+[thick,smooth,")
            fichier.write(style[j].addplot)
            if hasattr(style[j], "repeat"):
                if style[j].repeat == False:
                    fichier.write("]\n")
                else:
                    fichier.write(
                        "mark repeat = " + str(100) + "]\n"
                    )  # int(C.x.size/10)
            else:
                fichier.write("mark repeat = " + str(100) + "]\n")  # int(C.x.size/10)
            fichier.write("coordinates {")

            for i in np.arange(len(C.x), step=C.x.size // 1000 + 1):
                fichier.write("(" + str(C.x[i]) + "    ,    " + str(C.y[i]) + "    )\n")
            fichier.write("};\n")
        fichier.write("\\end{" + style[j].axis_type + "}")
        j += 1

    fichier.write(
        """\\end{tikzpicture}
\\endpgfgraphicnamed
\\end{document}
 """
    )
    fichier.close()

    #    fichier = open(directory + name + '_compil.bat','w')
    fichier = open(name + "_compil.bat", "w")
    fichier.write("pdflatex " + name + ".tex\n")
    fichier.write("pdflatex --jobname " + name + "_ " + name + ".tex\n")
    fichier.close()
    os.chdir(base_directory)

    if compil == True:

        os.chdir(directory)
        #        os.spawnl(os.P_NOWAIT, name + '_compil.bat')
        #        print("pdflatex", directory +  name  +".tex")
        #        print("pdflatex","--jobname", name +"_ ",   name  +".tex")
        #        os.spawnl(os.P_NOWAIT, name + '_compil.bat')
        subprocess.Popen(name + "_compil.bat")

        #        subprocess.Popen(["pdflatex", name  +".tex"])
        #        subprocess.Popen(["pdflatex","--jobname", name +"_ ", name  +".tex"])

        #        os.system( name + '_compil.bat')
        #        os.chdir('..\\')
        os.chdir(base_directory)


#    if legende != None:
#        Fichier.write("%s%s%s%s%s%s%s\n" % ("legend entries={",legende[0],',',legende[1],',',legende[2],"},"))
#        Fichier.write("%s%s%s%s%s\n" % ("legend style={at={(axis cs:",str(pos[0]),",",str(pos[1]),")},anchor=south west},"))
##
#    Fichier.write("%s%s%s\n" % ("xlabel={",axes[0],"},"))
#    Fichier.write("%s%s%s\n" % ("ylabel={",axes[1],"},"))
#

#    Fichier.write("%s%s%s\n" % ("\\addplot+[thick,mark=none] table [x expr=\\thisrowno{0}, y expr=\\thisrowno{1}] {",str(liste_adresse[i]),"};"))


# 	Fichier.write("\\documentclass[12pt]{article}\n")
# 	Fichier.write("\\usepackage{tikz}\n")
# 	Fichier.write("\\usepackage{pgfplots}\n")
# 	Fichier.write("%s%s%s\n" % ("\\pgfrealjobname{",nom_courbe,"}\n"))
#
# 	Fichier.write("\\begin{document}\n")
#
# 	Fichier.write("%s%s%s\n" % ("\\beginpgfgraphicnamed{",adresse + nom_courbe,"_}" ))
#
# 	Fichier.write("\\begin{tikzpicture}\n")
# 	Fichier.write("\\begin{axis}[\n")

# 	Fichier.write("axis lines = left,\n")
# 	Fichier.write("x axis line style={-|},\n")
# 	Fichier.write("y axis line style={-|},\n")

# 	Fichier.write("width=0.40*2.2\\textwidth,\n")
# 	Fichier.write("height=0.2472*2.2\\textwidth,\n")
# 	Fichier.write("cycle list name = linestyles*,\n")
# 	if legende != []:
# 		Fichier.write("%s%s%s%s%s%s%s\n" % ("legend entries={",legende[0],',',legende[1],',',legende[2],"},"))
# 		Fichier.write("%s%s%s%s%s\n" % ("legend style={at={(axis cs:",str(pos[0]),",",str(pos[1]),")},anchor=south west},"))
#

# 	Fichier.write("%s%s%s%s%s%s%s%s%s\n" % ("xmin=",str(xmin),",xmax=",str(xmax),",ymin=",str(ymin),",ymax=",str(ymax),","))
# 	Fichier.write("compat=1.3,]\n")
#
# 	print(liste_adresse)
#
# 	for i in range(0,len(liste_adresse)) :
# 		#Fichier.write("%s%s%s\n" % ("\\addplot+[thick,mark=none] table [x expr=\\thisrowno{0}, y expr=\\thisrowno{1}] {",str(liste_adresse[i]),"};"))
# 		Fichier.write("%s%s%s\n" % ("\\addplot+[thick,mark=none] table {",str(adresse+liste_adresse[i]),"};"))
#
# 	Fichier.write("\\end{axis}\n")
# 	Fichier.write("\\end{tikzpicture}\n")
# 	Fichier.write("\\endpgfgraphicnamed\n")
# 	Fichier.write("\\end{document}\n")
# 	#os.system("complil_fig.bat")
# 	#os.system("dir")
# 	Fichier.close()
#
# 	#commandes pour faire générer les courbes par windows
# 	subprocess.Popen(["pdflatex", adresse + nom_courbe +".tex"])
# 	subprocess.Popen(["pdflatex","--jobname",adresse + nom_courbe +"_", adresse + nom_courbe +".tex"])
# 	#~ print(["pdflatex","--jobname",adresse + nom_courbe +"_", adresse + nom_courbe +".tex"])
# 	#subprocess.Popen(["AcroRd32.exe",nom_courbe +"_.pdf"])
#
# 	#{os.system("pdflatex Rayon_test.tex")
# 	#os.system("--jobname  Rayon_test_ Rayon_test.tex ")
# 	#pdflatex Figures.tex
# 	#pdflatex --jobname  Rb_KM_RP_FEM Figures.tex
# 	Fichier = open('compil_'+nom_courbe +'.bat','w')
# 	Fichier.write("%s%s%s\n" % ("pdflatex ",nom_courbe ,".tex"))
# 	Fichier.write("%s%s%s%s%s\n" % ("pdflatex --jobname ",nom_courbe ,"_ ",nom_courbe,".tex"))
# 	Fichier.close()
