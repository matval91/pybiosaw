"""
Sanity checks on the input to ascot
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import a4py.classes.ReadEQDSK
from plot.plot_ascotinput import plot_bstruct
import a5py.ascot5io.ascot5 as a5                                                                                                   

def check_ascotinput(fin='/home/vallar/JT60-SA/3D/bfield/003/3Dfield_forascot_onlyTF.h5'):
    """
    Created on Thu Dec  6 10:13:44 2018
    
    Plot the magnetic field read from a file that can be used as input for ascot5
    
    @author: vallar
    """
    f=a5.Ascot(fin) 
    bstruct = f.bfield.active.read()
    plot_bstruct(bstruct)
    return bstruct