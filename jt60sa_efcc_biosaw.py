"""
This script reads the biosaw output for JT-60SA EFCC and builds the magnetic field as requested
"""
import h5py
import numpy as np

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
    # The last grid point where the data is given is not phimax, but it is the point just before. periodicity...
    data["phimin"] = coils["phigrid"][0][0]
    data["phimax"] = coils["phigrid"][0][-1]
    data["nphi"] = len(coils["phigrid"][0][:]) - 1
    
    # Empty grid for storing the magnetic field (remove last dublicate phi point)
    BR   = np.zeros( (data["nR"], data["nphi"], data["nz"]) )
    Bphi = np.zeros( (data["nR"], data["nphi"], data["nz"]) )
    Bz   = np.zeros( (data["nR"], data["nphi"], data["nz"]) )

    print('Data read')

    return coils,data, BR, Bphi, Bz


def compute(coils,data, BR, Bphi, Bz,nmode, U, M, L, phases):
    """
    Scales the field computed by biosaw with the correct current
    And puts it with the correct angles
    """
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

    if M is not 0:
	    # Equatorial coils are located on certain angles, offset is 189 degrees                    
	    for comp in ["BR", "Bphi", "Bz"]:
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