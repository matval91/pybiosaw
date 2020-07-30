import a5py.preprocessing.markers_grid as markers_grid

fn='./ascot_markers_500keV.h5'

qid=None
rmin=2.5
rmax=4
nr=10

pitchmin=0.4
pitchmax=1.
npitch=10

z=-0.4

phi=0
energy=500e3

Rpitch= markers_grid.R_pitch(fn=fn, qid=qid, \
                 rmin=rmin, rmax=rmax, nr=nr , \
                 pitchmin=pitchmin, pitchmax=pitchmax, npitch=npitch,\
                 z=z,phi=phi,energy=energy)
desc=f'npart={npitch*nr}, rmin={rmin}, rmax={rmax}, nr={nr}, pitchmin={pitchmin}, pitchmax={pitchmax}, npitch={npitch}, z={z},\
 energy={energy}'

Rpitch.make_input(description=desc)
print('markers written to')
print(f'{fn}: {desc}')