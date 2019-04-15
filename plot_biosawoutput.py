"""
Tools to plot the biosaw output

"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from mpl_toolkits.mplot3d import Axes3D

def merge_biosaw(outfile='wholecoils.h5', coils=['TF','EFCC01','EFCC07','EFCC13']):
    """
    Function to merge the TF coils (h5 biosaw output) with the EFCC coils (h5 biosaw output)
    into a single h5 biosaw output
    """
    maindir = '/home/vallar/JT60-SA/3D/bfield/biosaw/'
    outfile = maindir+outfile
    outfcoil = h5py.File(outfile)
    if 'TF' in coils:
        fcoil = h5py.File(maindir+'TF/TFcoilout.h5',"r+")
    else:
        fcoil = h5py.File(maindir+'efcc/EFCC01/EFCC01coilout.h5',"r+")

    for i in ['/Rgrid', '/phigrid', '/zgrid']:
        outfcoil.create_dataset(i, data=fcoil[i])
    fcoil.close()

    for c in coils:
        if c[0]!='E':
            fcoil = h5py.File(maindir+'TF/TFcoilout.h5',"r+")
            outfcoil.create_group('/'+c+'coil')
            for field_id  in ['BR', 'Bphi', 'Bz', 'coil_x', 'coil_y', 'coil_z']:
                field=fcoil['/'+c+'coil/'+field_id]
                outfcoil.create_dataset('/'+c+'coil/'+field_id, data=field)
        else:
            try:
                fcoil = h5py.File(maindir+'efcc/'+c+'/'+c+'coilout.h5',"r+")
            except:
                print("No coil "+c+" available")
                continue
            outfcoil.create_group('/'+c)
            for field_id  in ['BR', 'Bphi', 'Bz', 'coil_x', 'coil_y', 'coil_z']:
                field=fcoil['/'+c+'coil/'+field_id]
                outfcoil.create_dataset('/'+c+'/'+field_id, data=field)

    fcoil.close()
    outfcoil.close()

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
    coilname=['EFCC01', 'EFCC02', 'EFCC03', 'EFCC04', 'EFCC05', 'EFCC06',/
              'EFCC07', 'EFCC08', 'EFCC09', 'EFCC10', 'EFCC11', 'EFCC12',/
              'EFCC13', 'EFCC14', 'EFCC15', 'EFCC16', 'EFCC17', 'EFCC18']
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



