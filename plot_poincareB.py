import matplotlib.pyplot as plt
import numpy as np

import a5py.ascot5io.ascot5 as a5
from a5py.ascotpy.ascotpy import Ascotpy
import utils.plot_utils as pu
pu.common_style()
col = ['k', 'r', 'b', 'g', 'y', 'm']

def plot_poincareB(asc_fname='',run='', pncrid=0):
    if asc_fname=='': asc_fname='ascot.h5'
    f = a5.Ascot(asc_fname)
    if run=='': run=f.active

    #plot in rho and theta
    fig=plt.figure(); ax=fig.add_subplot(111) 
    run.orbit.poincare('rho', 'thetamod', pncrid, markersize=2, axes=ax)
    ax.set_xlabel(r'$\rho_{POL}$'); ax.set_ylabel(r'$\theta$ [°]')
    fig.tight_layout()

    # plot in r and z
    fig=plt.figure(figsize=(6,10)); ax=fig.add_subplot(111)
    run.orbit.poincare('R', 'z', 0, markersize=2, axes=ax)
    b5 = Ascotpy(asc_fname)
    b5.init(bfield=f.bfield.active.get_qid())
    b5.plotseparatrix(np.linspace(0,5,100), pncrid, np.linspace(-3,3,100), 0, axes=ax)
    ax.set_xlabel(r'R[m]'); ax.set_ylabel(r'Z [m]')
    fig.tight_layout()

    plt.show()
