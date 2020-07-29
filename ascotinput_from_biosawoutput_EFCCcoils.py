"""
Combine magnetic single-coil magnetic perturbation from TF coils and EFCCs to
create 3D magnetic field input for ASCOT5.

This script assume the field data has 720 grid points representing [0, 360] grid
in phi. One grid interval is therefore 0.5 deg.

TODO Rotation could be put into its own function (that accepts scaling factor)
"""
import numpy as np
import h5py
import a5py.ascot5io.B_3DS as B_3D
import a5py.ascot5io.B_2DS as B_2D

import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, RectBivariateSpline
import scipy.interpolate as interp
import a4py.classes.ReadEQDSK as ReadEQDSK
import utils.cocos_transform as ct
from utils.plot_utils import define_colors
import scipy.io as io

def produce_fields(fncoils='/home/vallar/JT60-SA/3D/bfield/biosaw/wholecoils.h5', 
    fnout='/home/vallar/JT60-SA/3D/bfield/ITER-like RMPs/EFCC_UML_n3_56deg0deg77deg.h5', 
    eqd_fname='/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/Equil_JT60_prova01_e_refined_COCOS7.eqdsk',\
                   nmode=3, U=1., M=1., L=1., phases=[56,0,77], bnorm_calculation=1, cocos=7):
    """
    data, Bphi, BR, Bz, theta, new_phi, newBnorm=produce_fields(fncoils, fnout, nmode, U,M,L,phases)
    """
    if fnout!='' and fnout[-2:]!='h5':
            print('Careful, output should be h5. Exiting')
            return 

    # Reading 2D magnetic field
    _eq=ReadEQDSK.ReadEQDSK(eqd_fname); eq = ct.cocos_transform(_eq, cocos,5)
    if 'prova01' in eqd_fname or eqd_fname=='/home/vallar/WORK/JT-60SA/input/005/JT-60SA_scenario5_eqdsk':
        eq.psi = np.reshape(eq.psi,(_eq.nzbox,_eq.nrbox))
        R_eqd = np.linspace(eq.rboxleft, eq.rboxleft+eq.rboxlength, _eq.nrbox)
        Z_eqd = np.linspace(-eq.zboxlength/2., eq.zboxlength/2., _eq.nzbox)
        _p=interp.interp2d(R_eqd,Z_eqd,eq.psi)
        new_nrbox=300; new_nzbox=300;
        R_eqd = np.linspace(eq.rboxleft, eq.rboxleft+eq.rboxlength, new_nrbox)
        Z_eqd = np.linspace(-eq.zboxlength/2., eq.zboxlength/2., new_nzbox) 
        eq.psi = _p(R_eqd,Z_eqd)
        eq.nrbox=new_nrbox; eq.nzbox=new_nzbox;

    pRmin=eq.rboxleft; pRmax=eq.rboxleft+eq.rboxlength; pnR=eq.nrbox
    pzmin=-1.*eq.zboxlength/2.; pzmax=1.*eq.zboxlength/2.; pnz=eq.nzbox
    psiaxis=eq.psiaxis; psisepx=eq.psiedge;psiRz=eq.psi
    axisR=eq.Raxis; axisz=eq.Zaxis;

    _z = np.linspace(pzmin, pzmax, pnz); _R = np.linspace(pRmin, pRmax, pnR)

    coils, data, BR, Bphi, Bz = read_data(fncoils)
    data['zmin']=-2.5
    data['zmax']=2.5
    data['rmin']=1.7
    # For EFCCs choose n mode (toroidal), currents and phase differences between coil rows
    BR,Bphi,Bz = _compute(coils,data, BR, Bphi, Bz,nmode, U, M, L, phases)
    for i,el in enumerate(coils["Rgrid"][0]):
        Bphi[i,:,:] += eq.B0EXP*eq.R0EXP/el

    # You have to interpolate the psi so that the Rzgrid over which psi 
    # is defined is the same as the B fields
    param_psi = interp.interp2d(_R, _z, eq.psi, kind='cubic')
    z = np.linspace(data['zmin'],data['zmax'], data['nz']); 
    R = np.linspace(data['Rmin'],data['Rmax'], data['nR'])
    psiRz = param_psi(R, z)

    if fnout!='':
        print('Writing to ', fnout)
        # Write output, same grid for fields and psi
        B_3D.write_hdf5(fnout,
                        data["Rmin"], data["Rmax"], data["nR"],
                        data["zmin"], data["zmax"], data["nz"],
                        data["phimin"], data["phimax"], data["nphi"],
                        axisR, axisz, psiRz.T, psiaxis, psisepx,
                        BR, Bphi, Bz,
                        psi_rmin=data["Rmin"], psi_rmax=data["Rmax"], psi_nr=data["nR"], \
                        psi_zmin=data["zmin"], psi_zmax=data["zmax"], psi_nz=data["nz"])

    if bnorm_calculation==1:
        theta,new_phi,newBnorm = Bnormal(data, BR, Bz, eq)
        plot(data, Bphi, theta, new_phi, newBnorm, nmode, phases,eq.q[-5])
        save_bnorm(R, axisR, z, theta, new_phi, new_Bnorm)
    else:
        theta=0; new_phi=0; newBnorm=0


    return data, Bphi, BR, Bz, theta, new_phi, newBnorm

def _compute(coils,data, BR, Bphi, Bz,nmode, U, M, L, phases):
    # Current for upper, middle, lower set of coils
    Ucurr=30e3*U
    Mcurr=45e3*M
    Lcurr=30e3*L
    currents = [Ucurr, Mcurr, Lcurr]
    angles = np.array([50, 90, 150, 210, 290, 350])

    # Place coils on correct phi coordinates, and assign correct current for each
    # UPPER AND LOWER (evenly spaced)
    for coil in ["EFCC01", "EFCC13"]:
        for comp in ["BR", "Bphi", "Bz"]:
            # Lower and Upper coils can be rotated like TF coils
            for i in range(0, 6, 1):
                Bcomp = np.transpose( coils[coil][comp][:], (1, 0, 2) )

                rotation = np.mod(float(data['nphi'])/360.*(60*i+191), data['nphi'])
                rotation -= rotation%1
                rotation = int(rotation)
                rotidx = np.append(np.arange(rotation,data["nphi"]),
                                   np.arange(0,rotation))
                Bcomp = Bcomp[:, rotidx, :]

                #phival is the correct coil location (-10. for offset)
                phival = coils["phigrid"][0][:][rotation]
                if coil == "EFCC01":
                    I = currents[0] * np.cos( (nmode*phival + phases[0])*np.pi/180 )
                if coil == "EFCC13":
                    I = currents[2] * np.cos( (nmode*phival + phases[2])*np.pi/180 )
                #Just multiply by I because the field is computed as it is a 1A current
                if comp == "BR":
                    BR   = BR + I*Bcomp 
                if comp == "Bphi":
                    Bphi = Bphi + I*Bcomp
                if comp == "Bz":
                    Bz   = Bz + I*Bcomp
                #print(coil, I)

                    
    for comp in ["BR", "Bphi", "Bz"]:
        # Equatorial coils are located on certain angles, offset is 189 degrees
        for i in angles:
            Bcomp = np.transpose( coils["EFCC07"][comp][:], (1, 0, 2) )
            #in this rotation there is an offset given by the fact that here the TFcoil 01 is at +x
            rotation = np.mod(float(data['nphi'])/360.*(60*i+191), data['nphi'])
            rotation -= rotation%1
            rotation = int(rotation)
            rotidx = np.append(np.arange(rotation,data["nphi"]), np.arange(0,rotation))
            Bcomp = Bcomp[:, rotidx, :]

            phival = coils["phigrid"][0][rotation]
            I = currents[1] * np.cos( (nmode*phival + phases[1])*np.pi/180 )
            if comp == "BR":
                BR   = BR + I*Bcomp
            if comp == "Bphi":
                Bphi = Bphi + I*Bcomp
            if comp == "Bz":
                Bz   = Bz + I*Bcomp
    return BR, Bphi, Bz


def read_data(fncoils):
    # Read magnetic field grid, this is same for all coils
    coils = h5py.File(fncoils, "r")

    data = {} #This is the structure going to ascot

    data["Rmin"] = coils["Rgrid"][0][0]
    data["Rmax"] = coils["Rgrid"][0][-1]
    data["nR"] = len(coils["Rgrid"][0][:])
    
    data["zmin"] = coils["zgrid"][0][0]
    data["zmax"] = coils["zgrid"][0][-1]
    data["nz"] = len(coils["zgrid"][0][:])

    # Note phimin=0 and phimax = 2*pi so B(phi=phimin) = B(phi=phimax) in input data
    # The last grid point where the data is given is not phimax, but it is the point just before. fucking periodicity
    data["phimin"] = coils["phigrid"][0][0]
    data["phimax"] = coils["phigrid"][0][-1]
    data["nphi"] = len(coils["phigrid"][0][:]) - 1
    
    # Empty grid for storing the magnetic field (remove last dublicate phi point)
    # The last grid point where the data is given is not phimax, but it is the point just before. fucking periodicity
    BR   = np.zeros( (data["nR"], data["nphi"], data["nz"]) )
    Bphi = np.zeros( (data["nR"], data["nphi"], data["nz"]) )
    Bz   = np.zeros( (data["nR"], data["nphi"], data["nz"]) )

    print('Data read')

    return coils,data, BR, Bphi, Bz

def Bnormal(data, BR, Bz, eq, rho=1):
    """
    theta, new_phi, Bnorm = Bnormal(data, BR, Bz, _eq)
    """

    print("Finding Bnormal at separatrix")
    R=np.linspace(data['Rmin'], data['Rmax'], data['nR'])
    phi=np.linspace(data['phimin'], data['phimax'], data['nphi'])
    z=np.linspace(data['zmin'], data['zmax'], data['nz'])
    BRparam = interp.RegularGridInterpolator((R,phi,z),BR)
    Bzparam = interp.RegularGridInterpolator((R,phi,z),Bz)
    if rho==1:
        new_R = eq.R[0:-1]; new_z = eq.Z[0:-1];

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
    save_bnorm(new_R, eq.Raxis, new_z, theta, new_phi, Bnorm)
    return theta, new_phi, Bnorm

def save_bnorm(R, Raxis, z, theta, new_phi, Bnorm):
    """
    """
    print('save bnorm to Bnorm.mat')
    dict={'R':R, 'Raxis':Raxis, 'z':z, 'theta':theta, 'phi':new_phi, \
    'Bnorm':Bnorm}
    io.savemat('Bnorm.mat', dict)

def plot(data, Bphi, theta, new_phi, bnorm, nmode, phases, q):
    # Plot
    # plt.figure()
    # R=np.linspace(data['Rmin'], data['Rmax'], data['nR'])
    # phi=np.linspace(data['phimin'], data['phimax'], data['nphi'])
    # plt.plot(phi, np.squeeze(Bphi[110,:,90]))
    # angles = np.mod(np.linspace(-181, 159, num=18), 360)
    # for i in angles:
    #     plt.plot([i,i],[-0.015, 0.015], 'k')
    # anglesM = np.mod(np.array([50, 90, 150, 210, 290, 350])+169., 360)
    # for i in anglesM:
    #     plt.plot([i,i],[-0.015, 0.015], 'r')

    # plt.grid('on')
    # plt.show()

    # plt.figure()
    # plt.contour(phi,R, Bphi[:,:,int(data["nz"]/2)])
    # for i in angles:
    #     plt.plot([i,i],[min(R), max(R)], 'k')
    # for i in anglesM:
    #     plt.plot([i,i],[min(R), max(R)], 'r')
    # plt.grid('on')
        
    _,_,_, my_cmap, _ = define_colors()
    fig=plt.figure(); ax=fig.add_subplot(111)

    #levels = np.linspace(bnorm.min(), bnorm.max(), 12)*1e4
    theta*=180./np.pi
    new_theta = np.mod(theta-360, 360)-180.

    #f2=plt.figure(); ax2=f2.add_subplot(111); ax2.contourf(new_phi, new_theta, bnorm)
    CS=ax.pcolormesh(new_phi, new_theta, bnorm*1e4,cmap=my_cmap)# vmin=-100., vmax=100.); 
    CB=plt.colorbar(CS, label=r'B$_{perp}$ [Gauss]')
    ax.set_ylim([-180., 180.])
    delta = 360./q
    for i in range(6):
        ax.plot(new_phi, new_phi/q+delta*(i-3), 'k-', lw=2.3)
        ax.plot(new_phi, (new_phi+180)/q+delta*(i-3), 'k--', lw=2.3)
    
    ax.set_xlabel(r'$\phi$', fontsize='xx-large')
    ax.set_ylabel(r'$\theta$', fontsize='xx-large')
    ax.tick_params(labelsize='xx-large')
    #ax.grid('on')
    ax.set_title('n='+str(nmode)+', phases='+str(phases))
    fig.tight_layout()
    plt.show()

# # Write output for mars
# fout = h5py.File(fnout, 'w')
# for component in ['R', 'phi', 'z']:
#     fout.create_dataset('/'+component+'min', data=data[component+"min"])
#     fout.create_dataset('/'+component+'max', data=data[component+"max"])
#     fout.create_dataset('/n'+component, data=data["n"+component])
# label=['BR', 'Bphi', 'Bz']
# for i, component in enumerate([BR, Bphi, Bz]):
#     fout.create_dataset('/'+label[i], data=BR)
# fout.close()
