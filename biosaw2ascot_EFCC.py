"""
Combine magnetic single-coil magnetic perturbation from TF coils and EFCCs to
create 3D magnetic field input for ASCOT5.

This script assume the field data has 720 grid points representing [0, 360] grid
in phi. One grid interval is therefore 0.5 deg.

TODO Rotation could be put into its own function (that accepts scaling factor)
"""
import numpy as np
import a5py.ascot5io.B_3DS as B_3D
import a5py.ascot5io.B_2DS as B_2D
import a5py.ascot5io.ascot5 as a5
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d, RectBivariateSpline
import scipy.interpolate as interp
import a5py.preprocessing.ReadEQDSK as ReadEQDSK
import utils.cocos_transform as ct
from utils.plot_utils import define_colors
from bfield_rmp_tools import *
from jt60sa_efcc_biosaw import *
import pdb

def produce_fields(fncoils='/home/vallar/WORK/JT-60SA/3D/biosaw/efcc_output/wholecoils_EFCC540phi.h5', 
    fnout='/home/vallar/JT60-SA/3D/bfield/ITER-like/EFCC_UML_n3_56deg0deg77deg.h5', 
    eqd_fname='/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/Equil_JT60_prova01_e_refined_COCOS7.eqdsk',
    fname_ripple='',
    nmode=3, U=1., M=1., L=1., phases=[56,0,77], bnorm_calculation=1, cocos=7, rho_bnorm=1.):
    """
    data, Bphi, BR, Bz, theta, new_phi, newBnorm=produce_fields(fncoils, fnout, nmode, U,M,L,phases)
    """
    if fnout!='' and fnout[-2:]!='h5':
            print('Careful, output should be h5. Exiting')
            return 
    coils_desc=''
    if U: coils_desc=[coils_desc, 'U']
    if M: coils_desc=[coils_desc, 'M']
    if L: coils_desc=[coils_desc, 'L']
    
    desc = f'nmode={nmode}, coils={coils_desc}, phases={phases}'
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

    # Read biosaw data from file with the coils
    coils, data, BR, Bphi, Bz = read_data(fncoils)
    data['zmin']=-2.5
    data['zmax']=2.5
    data['rmin']=1.7
    #Compute the magnetic field scaling
    # For EFCCs choose n mode (toroidal), currents and phase differences between coil rows
    BR,Bphi,Bz = compute(coils,data, BR, Bphi, Bz,nmode, U, M, L, phases)

    if fname_ripple!='':
        # Read TF ripple file
        # fname_TFripple_field = '/home/vallar/WORK/JT-60SA/ascot5_reference_files/scenario3/ascot_scen3_TFfield_wplasma.h5'
        f  = a5.Ascot(fname_ripple) 
        bfield_ripple = f.bfield.B_3DS_3104593462.read()
        #interpolate the fields
        output_shape = np.shape(BR)
        # for ir = 
        #     interp2d(x, y, z, kind='linear')
        # sum the fields
        BR = BR+bfield_ripple['br']
        Bz = Bz+bfield_ripple['bz']
        Bphi = Bphi+bfield_ripple['bphi']  
    else:
        # add the standard eq field
        for i,el in enumerate(coils["Rgrid"][0]):
            Bphi[i,:,:] += -1.*eq.B0EXP*eq.R0EXP/el

    # You have to interpolate the psi so that the Rzgrid over which psi 
    # is defined is the same as the B fields
    param_psi = interp.interp2d(_R, _z, eq.psi, kind='cubic')
    z = np.linspace(data['zmin'],data['zmax'], data['nz']); 
    R = np.linspace(data['Rmin'],data['Rmax'], data['nR'])
    psiRz = param_psi(R, z)


    # Write output, same grid for fields and psi
    if fnout!='':
        print('Writing to ', fnout)
        B_3D.write_hdf5(fnout,
                        data["Rmin"], data["Rmax"], data["nR"],
                        data["zmin"], data["zmax"], data["nz"],
                        data["phimin"], data["phimax"], data["nphi"],
                        axisR, axisz, psiRz.T, psiaxis, psisepx,
                        BR, Bphi, Bz,
                        psi_rmin=data["Rmin"], psi_rmax=data["Rmax"], psi_nr=data["nR"], \
                        psi_zmin=data["zmin"], psi_zmax=data["zmax"], psi_nz=data["nz"], desc=desc)

    #Calculation of the Bfield normal to the flux surfaces
    rho_2d=np.sqrt((psiRz-psiaxis)/(psisepx-psiaxis))
    if bnorm_calculation==1:
        theta,new_phi,newBnorm = Bnormal(data, BR, Bz, eq, rho_bnorm, R, z, rho_2d)
        plot(data, Bphi, theta, new_phi, newBnorm, nmode, phases, eq.q[-5], rho_bnorm)
    else:
        theta=0; new_phi=0; newBnorm=0


    return data, Bphi, BR, Bz, theta, new_phi, newBnorm, R, z, rho_2d, eq

def plot(data, Bphi, theta, new_phi, bnorm, nmode, phases, q, rho_bnorm):
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
    ax.set_title(f'n={nmode}, phases={phases}, rho={rho_bnorm}')
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
