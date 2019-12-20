"""
Sanity checks on the input to ascot
"""

from __future__ import print_function
import h5py
import numpy as np
import matplotlib.pyplot as plt
import a4py.classes.ReadEQDSK
from plot.plot_ascotinput import plot_bstruct

def check_ascotinput(fin='/home/vallar/JT60-SA/3D/bfield/003/3Dfield_forascot_onlyTF.h5'):
    """
    Created on Thu Dec  6 10:13:44 2018
    
    Plot the magnetic field read from a file that can be used as input for ascot5
    
    @author: vallar
    """
    f=h5py.File(fin)
    id =  list(f['bfield'].keys())[0]
    if '3' not in id[0:5]:
        for i, el in enumerate(list(f['bfield'].keys())):
            print(i,el)
        nn = input('Print which id number you want \n')
        id =  list(f['bfield'].keys())[int(nn)]
    print('Id chosen is ' + id)
    bstruct = f['bfield/'+id]
    plot_bstruct(bstruct)
    return bstruct