import a5py.preprocessing.markers_grid as markers_grid

fn='./ascot_markers_85keV.h5'

#YOU HAVE TO SET CORRECTLY THE PARAMETERS

# qid=None
# rmin=2.5
# rmax=4.1
# nr=10

# pitchmin=0.4
# pitchmax=1.
# npitch=10

# z=-0.4

# phi=0
# energy=85

Rpitch= markers_grid.R_pitch(fn, qid, \
                 rmin, rmax, nr , \
                 pitchmin, pitchmax, npitch,\
                 z,phi,energy)
desc=f'npart={npitch*nr}, rmin={rmin}, rmax={rmax}, nr={nr}, pitchmin={pitchmin}, pitchmax={pitchmax}, npitch={npitch}, z={z},\
    energy={energy}'

Rpitch.make_input(description=desc)