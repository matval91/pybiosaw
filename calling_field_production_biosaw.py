"""
calling the magnetic field rmpplot production
"""
import biosaw2ascot_EFCC as b2a   
import os
#on spcpc
dir='/home/vallar/'
eqd_fname='/home/vallar/WORK/JT-60SA/ascot5_reference_files/scenario3/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq'; cocos=2

#on personal pc
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
    eqd_fname='/home/matval/WORK/JT-60SA/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq';cocos=2
dir+='WORK/JT-60SA/3D/biosaw/efcc_output/'                                             
fncoils=dir+'wholecoils_EFCC540phi.h5'

fncoils='/home/vallar/WORK/JT-60SA/ascot5_reference_files/wholecoils_EFCC540phi.h5'
eqd_fname='/home/vallar/WORK/JT-60SA/ascot5_reference_files/scenario3/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq'; cocos=2
#eqd_fname='/home/vallar/WORK/JT-60SA/input/005/JT-60SA_scenario5_eqdsk'
#eqd_fname='JT-60SA_scenario2_highden_eqdsk_chease_cocos13_smoothed.geq'
#eqd_fname='JT-60SA_scenario5_eqdsk_chease_cocos02_smoothed.geq'
nmode=1
U=1
L=1
M=0
phases=[0,0,150] #Careful! putting numbers close to 60 will produce a 0 field!
lab_coils=''
if U: lab_coils+='U'
if M: lab_coils+='M'
if L: lab_coils+='L'

fnout=f'./EFCC_{lab_coils}_n{nmode}_{phases[0]}deg{phases[1]}deg{phases[2]}deg_003.h5'

bnorm_calculation=1
rho_bnorm = 1
data, Bphi, BR, Bz, theta, new_phi, newBnorm, R, z, rho_2d, eq=\
b2a.produce_fields(fncoils, fnout, eqd_fname, '',\
	nmode,U, M, L,phases, bnorm_calculation, cocos, rho_bnorm)
# b2a.plot(data, Bphi, theta, new_phi, newBnorm, nmode, phases, 4, rho_bnorm)