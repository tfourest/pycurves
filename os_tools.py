# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 11:46:05 2018

@author: tfourest
"""

from __future__ import print_function

import os
import fileinput

# from timeit import default_timer as timer


def _mkdir_recursive(path):
    sub_path = os.path.dirname(path)
    if not os.path.exists(sub_path):
        _mkdir_recursive(sub_path)
    if not os.path.exists(path):
        os.mkdir(path)


def make_dir(directory):
    #    if not os.path.exists(directory):    #Tu remplaces chemin par le chemin complet
    #        os.mkdir(directory)
    try:
        #        os.mkdir(directory)
        _mkdir_recursive(directory)
    except OSError:
        pass
    assert os.path.isdir(directory)


def liste_type(directory, type):
    """
    return a list of the file of type in the given directory organised by alphabetical order
    
    type : a string of the type ex: .ptw, .txt
    """
    l = []
    liste_dossier = os.listdir(directory)
    for fichier in liste_dossier:
        if type in fichier:
            l.append(directory + fichier)
    l.sort()
    return l


def clear_all():
    """Clears all the variables from the workspace of the spyder application.
    In spyder3 there is already a button to do it
    """
    gl = globals().copy()
    for var in gl:
        if var[0] == "_":
            continue
        if "func" in str(globals()[var]):
            continue
        if "module" in str(globals()[var]):
            continue

        del globals()[var]


def replace_separator(filename, separator="//"):
    """
    replace the separator given by the os separator in filenames
    this is made for easier use of the same code on linux and window
    ---------------------------------------------------------------------------
    input:
        filename :  full adress
    ---------------------------------------------------------------------------
    output:
        new_filename : filename with the separator replaced by the os separator
    """
    sep = os.sep
    new_filename = filename.replace(separator, sep)
    return new_filename


def write_line(fwrite, line, lim=0):
    """
    Write a line in a file with an automatic line break 
    (possible to add a maximum number of characters per line).
   
    	 	
    :param fwrite: file name 
    :type fwrite: str
    
    :param line: line to write in the file
    :type line: str
        
    :param lim: maximum number of characters per line (default=0)
    :type lim: int     
    """
    if lim > 0:
        for l in split_line(line, lim):
            fwrite.write(l + "\n")
    else:
        fwrite.write(line + "\n")


def split_line(line, nb):
    """
    Split a single line into lines with a maximum number of characters

    :param line: line to split
    :type line: str
        
    :param nb: maximum number of characters per line
    :type nb: int     
    
    :return: Splitted lines
    """

    return [line[start : start + nb] for start in range(0, len(line), nb)]


def replace_value(file_old, file_new, chars, values):
    """
    Replace the characters in file by the values
    ---------------------------------------------------------------------------
    :param file-old: name of the reference file
    :type file_old: str
    
    :param file-old: name of the new file
    :type file_old: str
    
    :param chars: list of characters to replace
    :type chars: array of strings
    
    :param values: list of values
    :type values: array of strings
    
    """
    oldfile = open(file_old, "r")
    path, file = os.path.split(file_new)
    make_dir(path)
    newfile = open(
        file_new, "w"
    )  # on ouvre le fichier new pour ecrire dedans, ou on le cree s'il n'existe pas
    for ligne in oldfile:
        for k in range(len(chars)):
            ligne = ligne.replace(chars[k], values[k])
        newfile.write(ligne)
    oldfile.close()
    newfile.close()


def alias(dico, dico_aliases):
    """
    create aliases in a dictionnary
    inputs :
    dico = the dictionnary in which to create aliases
    dico_aliases the dictionnary with aliases such as {'alias' : 'base_name'}
    """
    for clef in dico_aliases:
        if dico_aliases[clef] in dico:
            dico[clef] = dico[dico_aliases[clef]]


def printProgressBar(
    iteration, total, prefix="", suffix="", end="\r", decimals=1, length=50, fill="█"
):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(
        100 * ((iteration + 1) / float(total))
    )
    filledLength = int(length * (iteration + 1) // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    print("\r%s |%s| %s%% %s" % (prefix, bar, percent, suffix), end=end)
    # Print New Line on Complete
    if iteration == total:
        print()


def printime(sec):

    dm = 60
    ds = 60
    sec = int(sec)
    h = sec / (dm * ds)
    m = (sec - h * dm * ds) / ds
    s = sec - h * dm * ds - m * ds

    return "%02d hour(s) %02d minute(s) %02d second(s)" % (h, m, s)
