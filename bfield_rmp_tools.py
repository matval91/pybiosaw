"""
Functions to compute the field normal to the flux surfaces
"""
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import utils.plot_utils as pu
import pdb

pu.common_style()

def Bnormal(data, BR, Bz, eq, rho=1, Rpsi=None, zpsi=None, rhoRz=None):
    """
    theta, new_phi, Bnorm = Bnormal(data, BR, Bz, _eq)
    """
    print(f"Finding Bnormal at rho={rho}")
    R=np.linspace(data['Rmin'], data['Rmax'], data['nR'])
    phi=np.linspace(data['phimin'], data['phimax'], data['nphi'])
    z=np.linspace(data['zmin'], data['zmax'], data['nz'])
    BRparam = interp.RegularGridInterpolator((R,phi,z),BR)
    Bzparam = interp.RegularGridInterpolator((R,phi,z),Bz)
    try:
        new_R, new_z = rz_rho(Rpsi, zpsi, rhoRz, rho)
    except:
        print('Error, using separatrix for Bnorm')
        new_R = eq.R[0:-1]; new_z = eq.Z[0:-1] 
    
    new_phi = np.linspace(0,359,360)
    theta = np.arctan2(new_z,new_R-eq.Raxis)
    Bnorm = np.zeros((np.shape(new_R)[0], np.shape(new_phi)[0]))
    newBR=0; newBz=0
    for ii,iR in enumerate(new_R):
        for jj, jphi in enumerate(new_phi):
            newBR = BRparam((iR,jphi,new_z[ii]))
            newBz = Bzparam((iR,jphi,new_z[ii]))
            try:
                vector_normal = [-new_z[ii]+new_z[ii-1], iR-new_R[ii-1]]
            except:
                vector_normal = [-new_z[ii]+new_z[end], iR-new_R[end]]
            field_norm = np.dot([newBR, newBz], vector_normal)/np.linalg.norm(vector_normal)
            #field_norm = np.sqrt(newBR**2+newBz**2)
            Bnorm[ii,jj] = field_norm
            newBR = 0; newBz=0
    save_bnorm(new_R, eq.Raxis, new_z, theta, new_phi, Bnorm, rho)
    return theta, new_phi, Bnorm

def save_bnorm(R, Raxis, z, theta, new_phi, Bnorm, rho):
    """
    """
    print(f'save bnorm to Bnorm_{rho}.mat')
    dict={'R':R, 'Raxis':Raxis, 'z':z, 'theta':theta, 'phi':new_phi, \
    'Bnorm':Bnorm, 'rho':rho}
    io.savemat(f'Bnorm_{rho}.mat', dict)

def rz_rho(R, z, rhoRz, rho):
    """
    Finds R and Z of a surface of constant rho
    """
    cs = plt.contour(R, z, rhoRz, [rho]);
    RZsurface=cs.collections[0].get_paths()[0].vertices; plt.close('all')
    #must put the array so that the indexing starts from outer midplane
    ind=np.arange(1,np.size(RZsurface,0),6); ind=ind[::-1]
    try:
        ind_toroll = np.argwhere(abs(RZsurface[ind,1])<0.01)[0]
    except:
        ind_toroll = np.argwhere(abs(RZsurface[ind,1])<0.06)[0]
    rrr, zzz= RZsurface[ind,0], RZsurface[ind,1]
    return np.roll(rrr, -ind_toroll), np.roll(zzz, -ind_toroll)

def plot_bnorm(theta, phi, bnorm, q=0, rho_bnorm=0):
    _,_,_, my_cmap, _ = pu.define_colors()
    fig=plt.figure(); ax=fig.add_subplot(111)

    #levels = np.linspace(bnorm.min(), bnorm.max(), 12)*1e4
    theta*=180./np.pi
    new_theta = np.mod(theta-360, 360)-180.

    #f2=plt.figure(); ax2=f2.add_subplot(111); ax2.contourf(new_phi, new_theta, bnorm)
    CS=ax.pcolormesh(phi, new_theta, bnorm*1e4,cmap=my_cmap)# vmin=-100., vmax=100.); 
    CB=plt.colorbar(CS, label=r'B$_{perp}$ [Gauss]')
    ax.set_ylim([-180., 180.])
    if q!=0:
        delta = 360./q
        for i in range(6):
            ax.plot(phi, phi/q+delta*(i-3), 'k-', lw=2.3)
            ax.plot(phi, (phi+180)/q+delta*(i-3), 'k--', lw=2.3)
    ax.set_xlabel(r'$\phi$', fontsize='xx-large')
    ax.set_ylabel(r'$\theta$', fontsize='xx-large')
    ax.tick_params(labelsize='xx-large')
    #ax.grid('on')
    ax.set_title(f'$\\rho={rho_bnorm}, q={np.abs(q):.2f}$')
    fig.tight_layout()
    plt.show()