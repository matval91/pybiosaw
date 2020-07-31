import sys, os
sys.path.append('/home/matval/WORK/pythonscripts')
import utils.plot_utils as pu
import a5py.ascot5io.ascot5 as a5
import numpy as np
import matplotlib.pyplot as plt
import a5py.marker.evaluate as evaluate
from a5py.ascotpy.ascotpy import Ascotpy
from plot.plot_ascotinput import *
pu.common_style()
dir='/home/vallar/'
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
dir+='WORK/ASCOT/runs/SA_003/ripple/tang';
a=a5.Ascot(f'{dir}/ascot.h5')
run=a.active
#B field
wall=a.wall.wall_2D_3087769866
wall=wall.read()
inistate=run.inistate.read()
#Read RMPs
dir='/home/vallar/'
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
dir+='WORK/ASCOT/runs/SA_003/RMP_n3_120_UL/RMP_plasmaresponse/' 
f=a5.Ascot(f'{dir}bfield_ascot_mars_3D_iterlike.h5') 
b5 = Ascotpy(f'{dir}bfield_ascot_mars_3D_iterlike.h5')

#f=a5.Ascot(f'{dir}NORM_bfield_ascot_mars_3D_iterlike.h5') 
# b5 = Ascotpy(f'{dir}NORM_bfield_ascot_mars_3D_iterlike.h5')
bstruct = f.bfield.active.read()
R,z, phi = create_grid(bstruct)

# plot r,z
fig, ax=plot_poloidalmax(bstruct)
b5.init(bfield=f.bfield.active.get_qid())
b5.plotseparatrix(R,0,z,0, axes=ax)
ax.plot(wall['r'], wall['z'], 'y', lw=3)
ax.scatter(inistate['r'][1:200], inistate['z'][1:200], marker='x', color='c', alpha=0.4)
fig.set_figheight(10); fig.set_figwidth(7)
b5.free(bfield=True)
ax.set_xlabel(r'$R [m]$')
ax.set_ylabel(r'$z [m]$')
ax.axis('equal')
fig.tight_layout()
plt.show()