"""
Sanity checks on the input to ascot
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import a4py.classes.ReadEQDSK
from plot.plot_ascotinput import *
import a5py.ascot5io.ascot5 as a5                                                                                                   
from a5py.ascotpy.ascotpy import Ascotpy

def check_ascotinput(fin='/home/vallar/JT60-SA/3D/bfield/003/3Dfield_forascot_onlyTF.h5'):
    """
    Created on Thu Dec  6 10:13:44 2018
    
    Plot the magnetic field read from a file that can be used as input for ascot5
    
    @author: vallar
    """
    f=a5.Ascot(fin) 
    b5 = Ascotpy(fin)
    bstruct = f.bfield.active.read()
    R,z, phi = create_grid(bstruct)
    # plot_bstruct(bstruct)
    fig, ax=plot_poloidalmax(bstruct)
    b5.init(bfield=f.bfield.active.get_qid())
    b5.plotseparatrix(R,0,z,0, axes=ax)
    b5.free(bfield=True)
    fig.set_figheight(10); fig.set_figwidth(8)
    fig.tight_layout()
    return bstruct