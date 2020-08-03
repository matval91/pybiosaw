"""
Function to read the plasma response as received from Leonardo Pigatto
"""
import h5py
import numpy as np

def read_data(fnfield, real=False):
    print(real)
    # Read magnetic field grid, this is same for all coils
    bfield_file = h5py.File(fnfield, "r")
	#['Bphi', 'Br', 'Bz', 'R', 'Z', 'phi']
	# Bphi, Br and Bz have shape (Z,R,phi), and they are complex
    data = {} #This is the structure going to ascot

    data["Rmin"] = bfield_file["R"][0]
    data["Rmax"] = bfield_file["R"][-1]
    data["nR"] = len(bfield_file["R"][:])
    
    data["zmin"] = bfield_file["Z"][0]
    data["zmax"] = bfield_file["Z"][-1]
    data["nz"] = len(bfield_file["Z"][:])

    # Note phimin=0 and phimax = 2*pi so B(phi=phimin) = B(phi=phimax) in input data
    # The last grid point where the data is given is not phimax, but it is the point just before. fucking periodicity
    data["phimin"] = bfield_file["phi"][0]*180./np.pi
    data["phimax"] = bfield_file["phi"][-1]*180./np.pi
    data["nphi"] = len(bfield_file["phi"][:])
    
    # Empty grid for storing the magnetic field (remove last dublicate phi point)
    # The last grid point where the data is given is not phimax, but it is the point just before. fucking periodicity
    BR   = np.transpose(bfield_file["Br"][:,:,:], [1,2,0])
    Bphi = np.transpose(bfield_file["Bphi"][:,:,:], [1,2,0])
    Bz   = np.transpose(bfield_file["Bz"][:,:,:], [1,2,0])
    if real:
    	BR = np.real(BR); Bphi=np.real(Bphi); Bz = np.real(Bz)
    else:
    	BR = np.abs(BR); Bphi=np.abs(Bphi); Bz = np.abs(Bz)
    print('Data read')
    return bfield_file,data, BR, Bphi, Bz