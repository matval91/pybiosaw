"""
Combine magnetic single-coil magnetic perturbation from TF coils and EFCCs to
create 3D magnetic field input for ASCOT5.

This script assume the field data has 720 grid points representing [0, 360] grid
in phi. One grid interval is therefore 0.5 deg.

TODO Rotation could be put into its own function (that accepts scaling factor)
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d


fncoils = '/home/vallar/JT60-SA/3D/bfield/biosaw/wholecoils.h5'
# Read magnetic field grid, this is same for all coils
coils = h5py.File(fncoils, "r")

data = {} #This is the structure going to ascot
print('Data read')
fig=plt.figure(); ax=fig.add_subplot(111)
# Place coils on correct phi coordinates, and assign correct current for each
# UPPER AND LOWER (evenly spaced)
for coil in ["EFCC01", "EFCC13"]:
    # Lower and Upper coils can be rotated like TF coils
    x,y = coils[coil]["coil_x"].value, coils[coil]["coil_y"].value
    
    for i in range(0, 6, 1):
        rotation = np.mod(60*i+191,360)
        rotation = rotation*np.pi/180.
        newx = x*np.cos(rotation)+y*np.sin(rotation)
        newy = -x*np.sin(rotation)+y*np.cos(rotation)
        ax.plot(newx[0],newy[0], 'k')
        ax.text(newx[0, int(len(newx[0])/2.)]*1.1, newy[0, int(len(newy)/2.)]*1.1, str(i))

for i in range(18):
    ax.plot([0,4.5*np.cos((179-i*20)*np.pi/180.)], [0, 4.5*np.sin((179-i*20)*np.pi/180.)], 'k')


ax.axis('equal')
ax.grid('on')
ax.set_title('from coil rotation - UL')

fig=plt.figure(); ax=fig.add_subplot(111)
coilname = ['EFCC01',  'EFCC02',  'EFCC03',  'EFCC04', \
 'EFCC05',  'EFCC06']
x={}; y={}; z={}
f=h5py.File('/home/vallar/JT60-SA/3D/bfield/COILS/EFCC/efccdata.h5')
for i, coil in enumerate(coilname):
    x[coil] = f[coil+'/coil_x'].value[:,0]
    y[coil] = f[coil+'/coil_y'].value[:,0]
    z[coil] = f[coil+'/coil_z'].value[:,0]
    if i in [0,1,2,3,4,5]: col='g'
    elif i in [6,7,8,9,10,11]: col='r'
    else: col='b'

    ax.plot(x[coil], y[coil], color=col, lw=2.5)
    ax.text(x[coil][int(len(x[coil])/2.)]*1.1, y[coil][int(len(y[coil])/2.)]*1.1, coil)
for i in range(18):
    ax.plot([0,4.5*np.cos((179-i*20)*np.pi/180.)], [0, 4.5*np.sin((179-i*20)*np.pi/180.)], 'k')

ax.set_title('From data - UL')
ax.axis('equal')
ax.grid('on')

fig=plt.figure(); ax=fig.add_subplot(111)
coilname = ['EFCC01',  'EFCC02',  'EFCC03',  'EFCC04', \
 'EFCC05',  'EFCC06']
x={}; y={}; z={}
for i, coil in enumerate(coilname):
    x[coil] = f[coil+'/coil_x'].value[:,0]
    y[coil] = f[coil+'/coil_y'].value[:,0]
    ## THIS ROTATION WORKS
    rotation = (+180+11)*np.pi/180. #+ means clockwise rotation
    ##
    xtmp=np.copy(x[coil])
    ytmp=np.copy(y[coil])
    
    x[coil] = xtmp*np.cos(rotation)+ytmp*np.sin(rotation)
    y[coil] = -xtmp*np.sin(rotation)+ytmp*np.cos(rotation)

    if i in [0,1,2,3,4,5]: col='g'
    elif i in [6,7,8,9,10,11]: col='r'
    else: col='b'

    ax.plot(x[coil], y[coil], color=col, lw=2.5)
    ax.text(x[coil][int(len(x[coil])/2.)]*1.1, y[coil][int(len(y[coil])/2.)]*1.1, coil)
for i in range(18):
    ax.plot([0,4.5*np.cos((179-i*20)*np.pi/180.)], [0, 4.5*np.sin((179-i*20)*np.pi/180.)], 'k')

ax.set_title('From data, 191 deg CW rotated - UL')
ax.axis('equal')
ax.grid('on')


coil="EFCC07"
fig=plt.figure(); ax=fig.add_subplot(111)
x,y = coils[coil]["coil_x"].value, coils[coil]["coil_y"].value
#these angles are the one in the unwrapped phi angles EFCC configuration, starting from 0
angles = np.array([50, 90, 150, 210, 290, 350])
#the angles in the true ref system (AT) are 180+angles(i)-11
for i,el in enumerate(angles):
    rotation = np.mod(el-1+180, 360)
    rotation = -1.*rotation*np.pi/180.
    newx = x*np.cos(rotation)+y*np.sin(rotation)
    newy = -x*np.sin(rotation)+y*np.cos(rotation)
    ax.plot(newx[0],newy[0], 'k')
    ax.text(newx[0, int(len(newx[0])/2.)]*1.1, newy[0, int(len(newy)/2.)]*1.1, str(12-i))

for i in range(18):
    ax.plot([0,4.5*np.cos((179-i*20)*np.pi/180.)], [0, 4.5*np.sin((179-i*20)*np.pi/180.)], 'k')

ax.axis('equal')
ax.grid('on')
ax.set_title('from coil rotation - M')
coils.close()


fig=plt.figure(); ax=fig.add_subplot(111)
coilname = [ 'EFCC07',  'EFCC08',  'EFCC09', \
 'EFCC10',  'EFCC11',  'EFCC12']
x={}; y={}; z={}
for i, coil in enumerate(coilname):
    x[coil] = f[coil+'/coil_x'].value[:,0]
    y[coil] = f[coil+'/coil_y'].value[:,0]
    z[coil] = f[coil+'/coil_z'].value[:,0]
    if i in [0,1,2,3,4,5]: col='g'
    elif i in [6,7,8,9,10,11]: col='r'
    else: col='b'

    ax.plot(x[coil], y[coil], color=col, lw=2.5)
    ax.text(x[coil][int(len(x[coil])/2.)]*1.1, y[coil][int(len(y[coil])/2.)]*1.1, coil)

for i in range(18):
    ax.plot([0,4.5*np.cos((179-i*20)*np.pi/180.)], [0, 4.5*np.sin((179-i*20)*np.pi/180.)], 'k')
ax.set_title('From data - M')
ax.axis('equal')
ax.grid('on')

fig=plt.figure(); ax=fig.add_subplot(111)
x={}; y={}; z={}
for i, coil in enumerate(coilname):
    x[coil] = f[coil+'/coil_x'].value[:,0]
    y[coil] = f[coil+'/coil_y'].value[:,0]
    ## THIS ROTATION WORKS
    rotation = (+180+11)*np.pi/180. #+ means clockwise rotation
    ##
    xtmp=np.copy(x[coil])
    ytmp=np.copy(y[coil])
    
    x[coil] = xtmp*np.cos(rotation)+ytmp*np.sin(rotation)
    y[coil] = -xtmp*np.sin(rotation)+ytmp*np.cos(rotation)

    if i in [0,1,2,3,4,5]: col='g'
    elif i in [6,7,8,9,10,11]: col='r'
    else: col='b'

    ax.plot(x[coil], y[coil], color=col, lw=2.5)
    ax.text(x[coil][int(len(x[coil])/2.)]*1.1, y[coil][int(len(y[coil])/2.)]*1.1, coil)
for i in range(18):
    ax.plot([0,4.5*np.cos((179-i*20)*np.pi/180.)], [0, 4.5*np.sin((179-i*20)*np.pi/180.)], 'k')

ax.set_title('From data, , 191 deg CW rotated - M')
ax.axis('equal')
ax.grid('on')


f.close()
plt.show()



