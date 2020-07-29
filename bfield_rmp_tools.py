"""
Functions to compute the field normal to the flux surfaces
"""
import numpy as np
import scipy.io as io
import matplotlib.pyplot as plt
import scipy.interpolate as interp

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
    ll=cs.collections[0].get_paths()[0].vertices
    #must put the array so that the indexing starts from outer midplane
    ind=np.arange(1,np.size(ll,0),6); ind=ind[::-1]
    try:
        ind_toroll = np.argwhere(abs(ll[ind,1])<0.01)[0]
    except:
        ind_toroll = np.argwhere(abs(ll[ind,1])<0.05)[0]
    rrr, zzz= ll[ind,0], ll[ind,1]
    return np.roll(rrr, -ind_toroll), np.roll(zzz, -ind_toroll)