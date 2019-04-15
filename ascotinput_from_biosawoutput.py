"""
Combine magnetic single-coil magnetic perturbation from TF coils and EFCCs to
create 3D magnetic field input for ASCOT5.

This script assume the field data has 720 grid points representing [0, 360] grid
in phi. One grid interval is therefore 0.5 deg.

TODO Rotation could be put into its own function (that accepts scaling factor)
"""
import numpy as np
import h5py
import a5py.ascot5io.B_3D as B_3D
import a5py.ascot5io.B_2D as B_2D

import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
import ReadEQDSK
import a4py.utils.cocos_transform as ct

#Reading eqdsk with 2D fieldEquil_JT60_prova01_e_refined_COCOS7
_eq=ReadEQDSK.ReadEQDSK('/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/EQDSK_COCOS_02.OUT')
#/home/vallar/JT60-SA/003/eqdsk_fromRUI_20170715_SCENARIO3/Equil_JT60_prova01_e_refined_COCOS7.eqdsk')
eq = ct.cocos_transform(_eq, 2,3)
# Path to file with the coil data and to output file
#fncoils = '/l/sarkimk1/misc/JT60SA_3Dfield.h5'
#fnout = '/l/sarkimk1/misc/ascot5.h5'
fncoils = '/home/vallar/JT60-SA/3D/bfield/JT60SA_3Dfield.h5'
fnout = '/home/vallar/JT60-SA/3D/bfield/ascot.h5'
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
Raxis = eq.Raxis
zaxis = eq.Zaxis
Bphi0 = eq.B0EXP

print('Data read')

# Read field components from TF coils, copy it 17 times and rotate each "coil"
# so that there is 360 deg / 18 intervals between them.
for comp in ["BR", "Bphi", "Bz"]:
    print(comp)
    # We assume that we can rotate the grid toroidally by 20 degree intevals
    # otherwise we would have to interpolate
    for i in range(0, 18, 1):
        Bcomp = np.transpose( coils["TFcoil"][comp][:], (1, 0, 2) )

        # First TF coil is 9 degrees off, the rotation is in wrong direction so 
        # you should use 360-9. BUT the spacing between the coils is 20 deg, so you
        # can just rotate it of 20-9=11 degrees.
        # Do not think of coils labeling, at the end they will be in the correct place
        offset = 11*2 #2 is because the grid resolution is 0.5 deg
        rotation = 20*2*i + offset
        #this is just a "periodic shift" of the array
        # instead of being [0...i,j,...N] it is [j....N,1....i]
        rotidx = np.append(np.arange(rotation,data["nphi"]),
                           np.arange(0,rotation))
        #with rotidx you rotate the whole field from one coil
        Bcomp = Bcomp[:, rotidx, :]

        if comp == "BR":
            BR   = BR + Bcomp
        if comp == "Bphi":
            Bphi = Bphi + Bcomp
        if comp == "Bz":
            Bz   = Bz + Bcomp

# TF coil perturbation is in arbitrary units, scale to Teslas
# Find normalization factor by mathing average ripple on axis (in a.u.) to Bphi0
Bphi_avg = interp2d(coils["zgrid"][0][:], coils["Rgrid"][0][:], Bphi[:,10+9*2,:])
TF_scaling_factor = Bphi0/Bphi_avg(zaxis, Raxis)

BR   = BR*TF_scaling_factor
Bphi = Bphi*TF_scaling_factor
Bz   = Bz*TF_scaling_factor

print('TF coils field computed')

# For EFCCs choose n mode (toroidal), currents and phase differences between coil rows
nmode  = 3
# Current for upper, middle, lower set of coils
Ucurr=30e3*0.
Mcurr=45e3*0.
Lcurr=30e3*0.
currents = [Ucurr, Mcurr, Lcurr]
phases = np.array([10, 0, 20]) # In degrees

# Place coils on correct phi coordinates, and assign correct current for each
# UPPER AND LOWER (evenly spaced)
for coil in ["EFCC01", "EFCC13"]:
    for comp in ["BR", "Bphi", "Bz"]:
        # Lower and Upper coils can be rotated like TF coils
        for i in range(0, 6, 1):
            Bcomp = np.transpose( coils[coil][comp][:], (1, 0, 2) )

            rotation = np.mod(120*i-189*2,720)
            rotidx = np.append(np.arange(rotation,data["nphi"]),
                               np.arange(0,rotation))
            Bcomp = Bcomp[:, rotidx, :]

            #phival is the correct coil location (-10. for offset)
            phival = coils["phigrid"][0][:][rotation]-10.
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

for comp in ["BR", "Bphi", "Bz"]:
    # Equatorial coils are located on certain angles, offset is 189 degrees
    for i in np.mod(720-189*2-np.array([0, 60, 100, 160, 220, 300])*2,720):
        Bcomp = np.transpose( coils["EFCC07"][comp][:], (1, 0, 2) )

        rotation = i
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
plt.plot(np.squeeze(Bphi[110,:,90]))
plt.show()

# Reading 2D magnetic field
pRmin=eq.rboxleft; pRmax=eq.rboxleft+eq.rboxlength; pnR=eq.nrbox
pzmin=-1.*eq.zboxlength/2.; pzmax=1.*eq.zboxlength/2.; pnz=eq.nzbox
psiaxis=eq.psiaxis; psisepx=eq.psiedge;  psiRz=eq.psi
axisR=eq.Raxis; axisz=eq.Zaxis;

# Write output
B_3D.write_hdf5('./3Dfield_forascot_TFEFCC_UL.h5',
               data["Rmin"], data["Rmax"], data["nR"],
               data["zmin"], data["zmax"], data["nz"],
               data["phimin"], data["phimax"], data["nphi"],
               axisR, axisz, psiRz, psiaxis, psisepx,
               np.transpose(BR, (2,1,0)), np.transpose(Bphi, (2,1,0)), np.transpose(Bz,(2,1,0)),
               pRmin=pRmin, pRmax=pRmax, pnR=pnR, pzmin=pzmin, pzmax=pzmax, pnz=pnz)


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
