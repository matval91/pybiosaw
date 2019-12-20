"""
Tool to merge the biosaw output of TF and EFCC in the same file
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