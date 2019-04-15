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
from scipy.interpolate import interp2d
import ReadEQDSK
import a4py.utils.cocos_transform as ct

#Reading eqdsk with 2D fieldEquil_JT60_prova01_e_refined_COCOS7
#_eq=ReadEQDSK.ReadEQDSK('/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/EQDSK_COCOS_02.OUT')
#/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/Equil_JT60_prova01_e_refined_COCOS7.eqdsk')
#eq = ct.cocos_transform(_eq, 2,3)
# Path to file with the coil data and to output file
#fncoils = '/l/sarkimk1/misc/JT60SA_3Dfield.h5'
#fnout = '/l/sarkimk1/misc/ascot5.h5'
fncoils = '/home/vallar/JT60-SA/3D/bfield/JT60SA_3Dfield.h5'
fncoils = '/home/vallar/JT60-SA/3D/bfield/biosaw/wholecoils.h5'
fnout = '/home/vallar/JT60-SA/3D/bfield/EFCC_UML_n3_56deg0deg77deg.h5'
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

# Magnetic axis location and Bphi value on axis
# you should prompt these yourself
#Raxis = eq.Raxis
#zaxis = eq.Zaxis
#Bphi0 = eq.B0EXP

print('Data read')

# For EFCCs choose n mode (toroidal), currents and phase differences between coil rows
nmode  = 3
# Current for upper, middle, lower set of coils
Ucurr=30e3*1.
Mcurr=45e3*1.
Lcurr=30e3*1.
currents = [Ucurr, Mcurr, Lcurr]
phases = np.array([56, 0, 77]) # In degrees

# Place coils on correct phi coordinates, and assign correct current for each
# UPPER AND LOWER (evenly spaced)
for coil in ["EFCC01", "EFCC13"]:
    for comp in ["BR", "Bphi", "Bz"]:
        # Lower and Upper coils can be rotated like TF coils
        for i in range(0, 6, 1):
            Bcomp = np.transpose( coils[coil][comp][:], (1, 0, 2) )

            rotation = np.mod(2*(60*i+191),720)
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


angles = np.array([50, 90, 150, 210, 290, 350])
for comp in ["BR", "Bphi", "Bz"]:
    # Equatorial coils are located on certain angles, offset is 189 degrees
    for i in angles:
        Bcomp = np.transpose( coils["EFCC07"][comp][:], (1, 0, 2) )
        #in this rotation there is an offset given by the fact that here the TFcoil 01 is at +x
        rotation = np.mod(-2*(i+179), 360*2)
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


# TODO Add plasma response

# Plot
plt.figure()
R=np.linspace(data['Rmin'], data['Rmax'], data['nR'])
phi=np.linspace(data['phimin'], data['phimax'], data['nphi'])
plt.plot(phi, np.squeeze(Bphi[110,:,90]))
angles = np.mod(np.linspace(-181, 159, num=18), 360)
for i in angles:
    plt.plot([i,i],[-0.015, 0.015], 'k')
anglesM = np.mod(np.array([50, 90, 150, 210, 290, 350])+169., 360)
for i in anglesM:
    plt.plot([i,i],[-0.015, 0.015], 'r')

plt.grid('on')
plt.show()

plt.figure()
plt.contour(phi,R, Bphi[:,:,int(data["nz"]/2)])
for i in angles:
    plt.plot([i,i],[min(R), max(R)], 'k')
for i in anglesM:
    plt.plot([i,i],[min(R), max(R)], 'r')
plt.grid('on')
plt.show()




# # Reading 2D magnetic field
# pRmin=eq.rboxleft; pRmax=eq.rboxleft+eq.rboxlength; pnR=eq.nrbox
# pzmin=-1.*eq.zboxlength/2.; pzmax=1.*eq.zboxlength/2.; pnz=eq.nzbox
# psiaxis=eq.psiaxis; psisepx=eq.psiedge;  psiRz=eq.psi
# axisR=eq.Raxis; axisz=eq.Zaxis;

# Write output for mars

fout = h5py.File(fnout, 'w')
for component in ['R', 'phi', 'z']:
    fout.create_dataset('/'+component+'min', data=data[component+"min"])
    fout.create_dataset('/'+component+'max', data=data[component+"max"])
    fout.create_dataset('/n'+component, data=data["n"+component])
label=['BR', 'Bphi', 'Bz']
for i, component in enumerate([BR, Bphi, Bz]):
    fout.create_dataset('/'+label[i], data=BR)
fout.close()

# # Write output
# _R = np.linspace(data['Rmin'], data['Rmax'], data['nR'])
# _z = np.linspace(data['zmin'], data['zmax'], data['nz'])
# bsp = interp2d(_z, _R, Bphi[:,10,:])
# _R = np.linspace(pRmin, pRmax,pnR)
# _z = np.linspace(pzmin, pzmax, pnz)
# bphinew=bsp(_z,_R)
# B_2D.write_hdf5('./2Dfield_forascot_onlyTF.h5',
#                 pRmin, pRmax, pnR, pzmin, pzmax, pnz,
#                 axisR, axisz, psiRz, psiaxis, psisepx,
#                 np.zeros(psiRz.shape), bphinew.T, np.zeros(psiRz.shape))
