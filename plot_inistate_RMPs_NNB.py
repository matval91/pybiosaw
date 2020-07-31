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
dir+='WORK/ASCOT/runs/SA_003/2D/nnb/orbits'
a=a5.Ascot(f'{dir}/ascot.h5')
run=a.run_1593691660
#B field
wall=a.wall.wall_2D_3087769866
wall=wall.read()

#Read RMPs
RMP_fname='/home/matval/WORK/JT-60SA/RMPs/plasma_response/NORM_bfield_ascot_mars_3D_iterlike.h5'
f=a5.Ascot(RMP_fname) 
b5 = Ascotpy(RMP_fname)
bstruct = f.bfield.active.read()
R,z, phi = create_grid(bstruct)

#read ripple
dir='/home/vallar/'
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
dir+='WORK/ASCOT/runs/SA_003/ripple'
a=a5.Ascot(f'{dir}/ascot_TFripple.h5')
b5 = Ascotpy(f'{dir}/ascot_TFripple.h5')

# plot r,z
fig, ax=plot_poloidalmax(bstruct)
b5.init(bfield=f.bfield.active.get_qid())
b5.plotseparatrix(R,0,z,0, axes=ax)
ax.plot(wall['r'], wall['z'], 'k', lw=3)
ax.scatter(inistate['r'], inistate['z'], color='k')
fig.set_figheight(10); fig.set_figwidth(8)
b5.free(bfield=True)
ax.set_xlabel(r'$R [m]$')
ax.set_ylabel(r'$z [m]$')
ax.axis('equal')
fig.tight_layout()
plt.show()
