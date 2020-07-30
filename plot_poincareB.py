import h5py
import matplotlib.pyplot as plt
import numpy as np
import a4py.classes.ReadEQDSK as ReadEQDSK
import a5py.ascot5io.ascot5 as a5
import utils.plot_utils as pu
pu.common_style()
col = ['k', 'r', 'b', 'g', 'y', 'm']

def plot_poincareB_h5(run='', asc_fname='', eqd_fname=''):
    if asc_fname=='': asc_fname='ascot.h5'
    f = a5.Ascot(asc_fname)
    if run=='':
    	run=f.active

 	#plot in rho and theta
    fig=plt.figure(); ax=fig.add_subplot(111) 
    run.orbit.poincare('rho', 'thetamod', 0, markersize=2, axes=ax)
    ax.set_xlabel(r'$\rho_{POL}$'); ax.set_ylabel(r'$\theta$ [°]')
    fig.tight_layout()
    # plot in r and z
    fig=plt.figure(); ax=fig.add_subplot(111)
    run.orbit.poincare('R', 'z', 0, markersize=2, axes=ax)
    if eqd_fname!='':
        _eq=ReadEQDSK.ReadEQDSK(eqd_fname)
        plt.plot(_eq.R, _eq.Z, 'k-')
    ax.set_xlabel(r'R[m]'); ax.set_ylabel(r'Z [m]')
    fig.tight_layout()

    plt.show()





def plot_poincareB(run='', asc_fname='', eqd_fname=''):
    if asc_fname=='': asc_fname='ascot.h5'
    f = h5py.File(asc_fname)
    mlpol=f['results/run-'+str(run)+'/orbits/mlpol0']
    r=mlpol['R'].value; z=mlpol['z'].value; idm=mlpol['id'].value

    fig=plt.figure(); ax=fig.add_subplot(111) 
    c=ax.scatter(r,z, c=mlpol['B_phi'], marker='.'); ax.axis('equal'); 
    plt.colorbar(c)
    ax.set_xlabel(r'R[m]'); ax.set_ylabel(r'Z [m]')
    if eqd_fname!='':
        _eq=ReadEQDSK.ReadEQDSK(eqd_fname)
        plt.plot(_eq.R, _eq.Z, 'k-')

    mltor=f['results/run-'+str(run)+'/orbits/mltor0']
    rho = mltor['rho']; phi=np.mod(mltor['phi'], 360)
    fig=plt.figure(); ax=fig.add_subplot(111)
    fig.suptitle('Poincare plot on the midplane')
    ax.scatter(rho, phi, marker='.');
    ax.set_xlabel(r'$\rho$'); ax.set_ylabel(r'$\phi$ [°]')
    plt.show()
