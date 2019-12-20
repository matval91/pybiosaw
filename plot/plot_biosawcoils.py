"""
Tools to plot the biosaw output

"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from mpl_toolkits.mplot3d import Axes3D

def plot_coils(fncoils='/home/vallar/JT60-SA/3D/bfield/biosaw/wholecoils.h5'):
    """
    Plot all the JT60SA coils from a biosaw h5 output
    """
    plot_TFcoils(fncoils)
    ax=plt.gca()
    plot_EFCCcoils(fncoils,ax)

def plot_TFcoils(fncoils='/home/vallar/JT60-SA/3D/bfield/biosaw/wholecoils.h5', ax=0):
    """
    Plot  the JT60SA TF coils from a biosaw h5 output
    """
    col=['k', 'g', 'r', 'b']
    f = h5py.File(fncoils)
    x={}; y={}; z={}
    xnew={}; ynew={}
    if type(ax)==int:
        fig=plt.figure(); ax = fig.add_subplot(111, projection='3d')

    for i,coil in enumerate(['TFcoil']):
        x[coil] = f[coil+'/coil_x'].value[0,:]
        y[coil] = f[coil+'/coil_y'].value[0,:]
        z[coil] = f[coil+'/coil_z'].value[0,:]
        if coil == 'TFcoil':
            for n in range(18):
                angle=360./18.*n*np.pi/180.
                xnew[coil] =  x[coil]*np.cos(angle)- y[coil]*np.sin(angle)
                ynew[coil] =  x[coil]*np.sin(angle)+ y[coil]*np.cos(angle)
                ax.plot(xnew[coil], ynew[coil],zs= z[coil], color=col[i], lw=2.5)

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')

    plt.show()

def plot_EFCCcoils(fncoils='/home/vallar/JT60-SA/3D/bfield/biosaw/wholecoils.h5', ax=0):
    """
    Plot the JT60SA EFCC coils from a biosaw h5 output    
    """
    col=['k', 'g', 'r', 'b']
    coilname=['EFCC01', 'EFCC02', 'EFCC03', 'EFCC04', 'EFCC05', 'EFCC06', \
              'EFCC07', 'EFCC08', 'EFCC09', 'EFCC10', 'EFCC11', 'EFCC12', \
              'EFCC13', 'EFCC14', 'EFCC15', 'EFCC16', 'EFCC17', 'EFCC18']
    f=h5py.File(fncoils)
    if type(ax)==int:
        fig=plt.figure(); ax = fig.add_subplot(111, projection='3d')

    for i,coil in enumerate(coilname):
        try:
            x[coil] = f[coil+'/coil_x'].value[0,:]
            y[coil] = f[coil+'/coil_y'].value[0,:]
            z[coil] = f[coil+'/coil_z'].value[0,:]
        except:
            continue
        ax.plot(x[coil], y[coil], zs= z[coil], color=col[i], lw=2.5)

    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')

    plt.show()



