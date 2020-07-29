"""
Sanity checks on the input to ascot
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import a4py.classes.ReadEQDSK as ReadEQDSK
import a5py.ascot5io.ascot5 as ascot5
def plot_bstruct(bstruct):
    """
    Plots Bphi on midplane, BR and Bz at phi=0
    """
    #Creating grid
    nR = np.squeeze(bstruct['b_nr'][()])
    Rmin, Rmax = np.squeeze(bstruct['b_rmin'][()]), np.squeeze(bstruct['b_rmax'][()])
    R = np.linspace(Rmin, Rmax, nR)
    nz = np.squeeze(bstruct['b_nz'][()])
    zmin, zmax = np.squeeze(bstruct['b_zmin'][()]), np.squeeze(bstruct['b_zmax'][()])
    z = np.linspace(zmin, zmax, nz)
    nphi = np.squeeze(bstruct['b_nphi'][()])
    phimin, phimax = np.squeeze(bstruct['b_phimin'][()]), np.squeeze(bstruct['b_phimax'][()])
    phi = np.linspace(phimin, phimax, nphi)


    # Bphi has the shape (phi, z, R)
    Bphi = bstruct['bphi'][()]
    ind_midplane = np.argmin(z-0 <0)
    fig=plt.figure(figsize=(20,8)); 
    ax=fig.add_subplot(131); ax.set_title(r'$B_\phi$ on midplane')
    cs=ax.contour(R, phi, Bphi[:, :, ind_midplane].T, 40)
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'phi')
    fig.tight_layout()

    BR = bstruct['br'][()]
    ind_phi0 = np.argmin(phi-0 <0)
    ax=fig.add_subplot(132); ax.set_title(r'$B_R (\phi=0)$ ')
    cs=ax.contour(R, z, BR[:,ind_phi0, :].T, 40); 
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')
    ax.axis('equal')
    fig.tight_layout()

    Bz = bstruct['bz'][()]
    ind_phi0 = np.argmin(phi-0 <0)
    ax=fig.add_subplot(133); ax.set_title(r'$B_z (\phi=0)$ ')
    cs=ax.contour(R, z, Bz[:, ind_phi0, :].T, 40); 
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')
    ax.axis('equal')
    fig.tight_layout()
    
    plt.show()


def compare_3Dvs2D(fin_3d='/home/vallar/JT60-SA/3D/bfield/biosaw2ascot/ascot_TFfield_b0expr0exp.h5', \
    eqd_2d='/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/EQDSK_COCOS_02.OUT'):
    """
    Created on Thu Dec  6 10:13:44 2018

    Plots Bphi on midplane, Br, Bphi, Bz at phi=0

    Compares 
    3D field in a file that could be used as input to ascot5 
    VS
    2D field read in an eqdsk

    @author: vallar
    """
    #Read 3D file
    a5file = ascot5.Ascot(fin_3d)
    bstruct=a5file.bfield.active.read()

    #Read 2D file
    eq=ReadEQDSK.ReadEQDSK(eqd_2d)
    B2D = eq.B0EXP*eq.R0EXP/eq.R_grid


    #Creating grid
    nR = np.squeeze(bstruct['b_nr'])
    Rmin, Rmax = np.squeeze(bstruct['b_rmin']), np.squeeze(bstruct['b_rmax'])
    R = np.linspace(Rmin, Rmax, nR)
    nz = np.squeeze(bstruct['b_nz'])
    zmin, zmax = np.squeeze(bstruct['b_zmin']), np.squeeze(bstruct['b_zmax'])
    z = np.linspace(zmin, zmax, nz)
    nphi = np.squeeze(bstruct['b_nphi'])
    phimin, phimax = np.squeeze(bstruct['b_phimin']), np.squeeze(bstruct['b_phimax'])
    phi = np.linspace(phimin, phimax, nphi)


    # Bphi has the shape (phi, z, R)
    Bphi = bstruct['bphi']
    ind_midplane = np.argmin(z-0 <0)
    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_\phi$ on midplane')
    cs=ax.contour(R, phi, Bphi[:, :, ind_midplane].T, 40)
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'phi')
    fig.tight_layout()

    BR = bstruct['br']
    ind_phi0 = np.argmin(phi-0 <0)
    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_R (\phi=0)$ ')
    cs=ax.contour(R, z, BR[:,ind_phi0, :].T, 40); 
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')
    fig.tight_layout()

    Bz = bstruct['bz']
    ind_phi0 = np.argmin(phi-0 <0)
    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_z (\phi=0)$ ')
    cs=ax.contour(R, z, Bz[:,ind_phi0, :].T, 20);
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')
    fig.tight_layout()

    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_\phi (\phi=0)$ ')
    cs = ax.contour(R, z, Bphi[:,ind_phi0,:].T)
    cb=plt.colorbar(cs)
    ax.contour(R, z, Bphi[:, ind_phi0, :].T, [-1.*eq.B0EXP]);
    ax.plot([eq.R0EXP, eq.R0EXP],[-3., 3.])
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')

    plt.show()
