from __future__ import print_function
import h5py
import numpy as np
import matplotlib.pyplot as plt
import a4py.classes.ReadEQDSK

def check_ascotinput(fin='/home/vallar/JT60-SA/3D/bfield/003/3Dfield_forascot_onlyTF.h5'):
    """
    Created on Thu Dec  6 10:13:44 2018
    
    Plot the magnetic field read from a file that can be used as input for ascot5
    
    @author: vallar
    """
    f=h5py.File(fin)
    id =  list(f['bfield'].keys())[0]
    print('Id chosen is ' + id)
    bstruct = f['bfield/'+id]

    #Creating grid
    nR = bstruct['n_R'].value
    Rmin, Rmax = bstruct['R_min'].value, bstruct['R_max'].value
    R = np.linspace(Rmin, Rmax, nR)
    nz = bstruct['n_z'].value
    zmin, zmax = bstruct['z_min'].value, bstruct['z_max'].value
    z = np.linspace(zmin, zmax, nz)
    nphi = bstruct['n_phi'].value
    phimin, phimax = bstruct['phi_min'].value, bstruct['phi_max'].value
    phi = np.linspace(phimin, phimax, nphi)


    # Bphi has the shape (phi, z, R)
    Bphi = bstruct['B_phi'].value
    ind_midplane = np.argmin(z-0 <0)
    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_\phi$ on midplane')
    cs=ax.contour(R, phi, Bphi[:, ind_midplane, :], 40)
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'phi')
    fig.tight_layout()

    BR = bstruct['B_R'].value
    ind_phi0 = np.argmin(phi-0 <0)
    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_R (\phi=0)$ ')
    cs=ax.contour(R, z, BR[ind_phi0,:, :], 40); 
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')
    fig.tight_layout()

    Bz = bstruct['B_z'].value
    ind_phi0 = np.argmin(phi-0 <0)
    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_z (\phi=0)$ ')
    cs=ax.contour(R, z, Bz[ind_phi0,:, :], 40); 
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')
    fig.tight_layout()

    plt.show()

def compare_3Dvs2D(fin_3d='/home/vallar/JT60-SA/3D/bfield/biosaw2ascot/ascot_TFfield_b0expr0exp.h5', eqd_2d='/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/EQDSK_COCOS_02.OUT'):

    """
    Created on Thu Dec  6 10:13:44 2018

    Compares 
    3D field in a file that could be used as input to ascot5 
    VS
    2D field read in an eqdsk

    @author: vallar
    """
    #Read 3D file
    f=h5py.File(fin_3d)
    id =  list(f['bfield'].keys())[1]
    print('Id chosen is ' + id)
    bstruct = f['bfield/'+id]

    #Read 2D file
    eq=ReadEQDSK.ReadEQDSK(eq_2d)
    B2D = eq.B0EXP*eq.R0EXP/eq.R_grid


    #Creating grid
    nR = bstruct['n_R'].value
    Rmin, Rmax = bstruct['R_min'].value, bstruct['R_max'].value
    R = np.linspace(Rmin, Rmax, nR)
    nz = bstruct['n_z'].value
    zmin, zmax = bstruct['z_min'].value, bstruct['z_max'].value
    z = np.linspace(zmin, zmax, nz)
    nphi = bstruct['n_phi'].value
    phimin, phimax = bstruct['phi_min'].value, bstruct['phi_max'].value
    phi = np.linspace(phimin, phimax, nphi)


    # Bphi has the shape (phi, z, R)
    Bphi = bstruct['B_phi'].value
    ind_midplane = np.argmin(z-0 <0)
    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_\phi$ on midplane')
    cs=ax.contour(R, phi, Bphi[:, ind_midplane, :], 40)
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'phi')
    fig.tight_layout()

    BR = bstruct['B_R'].value
    ind_phi0 = np.argmin(phi-0 <0)
    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_R (\phi=0)$ ')
    cs=ax.contour(R, z, BR[ind_phi0,:, :], 40); 
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')
    fig.tight_layout()

    Bz = bstruct['B_z'].value
    ind_phi0 = np.argmin(phi-0 <0)
    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_z (\phi=0)$ ')
    cs=ax.contour(R, z, Bz[ind_phi0,:, :], 20);
    cb=plt.colorbar(cs)
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')
    fig.tight_layout()

    fig=plt.figure(); ax=fig.add_subplot(111); ax.set_title(r'$B_\phi (\phi=0)$ ')
    cs = ax.contour(R, z, Bphi[ind_phi0,:,:])
    cb=plt.colorbar(cs)
    ax.contour(R, z, Bphi[ind_phi0,:, :], [-1.*eq.B0EXP]);
    ax.plot([eq.R0EXP, eq.R0EXP],[-3., 3.])
    ax.set_xlabel(r'R [m]'); ax.set_ylabel(r'z [m]')

    plt.show()
