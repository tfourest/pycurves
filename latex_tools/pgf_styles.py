# -*- coding: utf-8 -*-
"""
Defines standard pgfstyles for use with the writeXcourbesPGF
"""


# from pgf_tools import *


class pgf_style:
    """
    class to stock a pgfplot style that shall be use in pycurves.latex_tools.writeXcourbesPGF 
    
      
    Attributes
    ----------
    plotcyclelist : string
        Define a cycling of curve style in pgfplot language
    axis : string
            it is what will be written after 
            \begin{tikzpicture}
            \begin{axis}[ ... ]    
    addplot : string or list
        Options to add pgfplot curves in the **\addplot** command 
    automatic_repeat : boolean
        Define if the pgfplot **mark_repeat** command is use automatically. By default it is considered True.
         
        
    Predifined Styles
    --------------
    There are some predifined style defined. The best way to know what they do is to try them.
        
        Tom_standart_BW_curves
        
        Classic_pgf_color_style
        
        Classic_pgf_color_style_SIG_EPS
        
        Custom_tom_pgf_color_style
        
        Custom_tom_pgf_color_style_red
        
        Custom_tom_pgf_color_style_only_marks
      
    """

    def __init__(self):
        pass


#        self.style = s


Tom_standart_BW_curves = pgf_style()
Tom_standart_BW_curves.plotcyclelist = ""
Tom_standart_BW_curves.axis = """axis lines = left,
x axis line style={-|},
y axis line style={-|},
width=0.40*2.2\\textwidth,
height=0.2472*2.2\\textwidth,
cycle list name = linestyles*,
legend pos=north west,
legend cell align={left},
compat=1.3,
scaled y ticks={base 10:0},
"""
Tom_standart_BW_curves.addplot = ""


TSBWC = Tom_standart_BW_curves
# ==============================================================================
# Classic pgf color style
# ==============================================================================

Classic_pgf_color_style = pgf_style()
Classic_pgf_color_style.plotcyclelist = ""
Classic_pgf_color_style.axis = """grid = both,
x axis line style={-|},
y axis line style={-|},
width=0.40*2.2\\textwidth,
height=0.2472*2.2\\textwidth,
legend style={at={(0.97,0.03)},anchor= south east},
legend cell align={left},
compat=1.3,
scaled y ticks={base 10:-3},
"""
Classic_pgf_color_style.addplot = ""


Classic_pgf_color_style_SIG_EPS = pgf_style()
Classic_pgf_color_style_SIG_EPS.plotcyclelist = ""
Classic_pgf_color_style_SIG_EPS.axis = """grid = both,
x axis line style={-|},
y axis line style={-|},
width=0.40*2.2\\textwidth,
height=0.2472*2.2\\textwidth,
legend style={at={(0.97,0.03)},anchor= south east},
legend cell align={left},
compat=1.3,
scaled y ticks={base 10:-3},
"""
Classic_pgf_color_style_SIG_EPS.addplot = ""

Custom_tom_pgf_color_style = pgf_style()
Custom_tom_pgf_color_style.plotcyclelist = """\pgfplotscreateplotcyclelist{mycolorlist}{%
blue,every mark/.append style={fill=.!75!white},mark=*\\\\%
blue!80,every mark/.append style={fill=.!75!white},dashed,mark=square*\\\\%
blue!60,every mark/.append style={fill=.!75!white},dotted,mark=triangle*\\\\%
blue!40,mark=star\\\\%
black!40,every mark/.append style={fill=.!75!white},dashed,mark=diamond*\\\\%
black!60,dotted,every mark/.append style={solid,fill=.!75!white},mark=*\\\\%
black!80!black,densely dashed,every mark/.append style={
solid,fill=.!75!white},mark=square*\\\\%
black,densely dashed,every mark/.append style={solid,fill=.!75!white},mark=otimes*\\\\%
%blue,densely dashed,mark=star,every mark/.append style=solid\\\\%
%red,densely dashed,every mark/.append style={solid,fill=.!75!white},mark=diamond*\\\\%
}
"""
Custom_tom_pgf_color_style.axis = """grid = both,
x axis line style={-|},
y axis line style={-|},
width=0.40*2.2\\textwidth,
height=0.2472*2.2\\textwidth,
legend style={at={(0.97,0.03)},anchor= south east},
legend cell align={left},
cycle list name=mycolorlist,
compat=1.3,
scaled y ticks={base 10:0},
scaled x ticks={base 10:2},
"""
Custom_tom_pgf_color_style.addplot = ""

Custom_tom_pgf_color_style_red = pgf_style()
Custom_tom_pgf_color_style_red.plotcyclelist = """\pgfplotscreateplotcyclelist{mycolorlist}{%
red,every mark/.append style={fill=.!75!white},mark=*\\\\%
red!80,every mark/.append style={fill=.!75!white},dashed,mark=square*\\\\%
red!60,every mark/.append style={fill=.!75!white},dotted,mark=triangle*\\\\%
red!40,mark=star\\\\%
black!40!red,every mark/.append style={fill=.!75!white},dashed,mark=diamond*\\\\%
black!60!red,dotted,every mark/.append style={solid,fill=.!75!white},mark=*\\\\%
black!80!red,densely dashed,every mark/.append style={
solid,fill=.!75!white},mark=square*\\\\%
black,densely dashed,every mark/.append style={solid,fill=.!75!white},mark=otimes*\\\\%
%blue,densely dashed,mark=star,every mark/.append style=solid\\\\%
%red,densely dashed,every mark/.append style={solid,fill=.!75!white},mark=diamond*\\\\%
}
"""
Custom_tom_pgf_color_style_red.axis = """grid = both,
x axis line style={-|},
y axis line style={-|},
width=0.40*2.2\\textwidth,
height=0.2472*2.2\\textwidth,
legend style={at={(0.97,0.03)},anchor= south east},
legend cell align={left},
cycle list name=mycolorlist,
compat=1.3,
scaled y ticks={base 10:0},
scaled x ticks={base 10:2},
"""
Custom_tom_pgf_color_style_red.addplot = ""


Custom_tom_pgf_color_style_only_marks = pgf_style()
Custom_tom_pgf_color_style_only_marks.plotcyclelist = """\pgfplotscreateplotcyclelist{mycolorlist}{%
blue,every mark/.append style={fill=.!75!white},mark=*\\\\%
blue!80,every mark/.append style={fill=.!75!white},dashed,mark=square*\\\\%
blue!60,every mark/.append style={fill=.!75!white},dotted,mark=triangle*\\\\%
blue!40,mark=star\\\\%
black!40,every mark/.append style={fill=.!75!white},dashed,mark=diamond*\\\\%
black!60,dotted,every mark/.append style={solid,fill=.!75!white},mark=*\\\\%
black!80!black,densely dashed,every mark/.append style={
solid,fill=.!75!white},mark=square*\\\\%
black,densely dashed,every mark/.append style={solid,fill=.!75!white},mark=otimes*\\\\%
%blue,densely dashed,mark=star,every mark/.append style=solid\\\\%
%red,densely dashed,every mark/.append style={solid,fill=.!75!white},mark=diamond*\\\\%
}
"""
Custom_tom_pgf_color_style_only_marks.axis = """grid = both,
x axis line style={-|},
y axis line style={-|},
width=0.40*2.2\\textwidth,
height=0.2472*2.2\\textwidth,
legend style={at={(0.97,0.97)},anchor= north east},
legend cell align={left},
cycle list name=mycolorlist,
compat=1.3,
scaled y ticks={base 10:0},
%scaled x ticks={base 10:2},
"""
Custom_tom_pgf_color_style_only_marks.addplot = "only marks,"
