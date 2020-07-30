import a5py.preprocessing.markers_grid as markers_grid

dir='/home/vallar/'
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
dir+='WORK/ASCOT/runs/SA_003/RMP_n3_120_UL/RMP_vacuum/nnb'
fn=f'{dir}/ascot.h5'
a=a5.Ascot(fn)
qid=a.run.active.get_qid()

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

Rpitch.contourf_attr('deltaenergy')

print('markers written to')
print(f'{fn}: {desc}')