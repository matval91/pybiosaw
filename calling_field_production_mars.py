"""
calling the magnetic field rmpplot production
"""
import mars2ascot as m2a   
import os
#on spcpc
dir='/home/vallar/'
eqd_fname='/home/vallar/WORK/JT-60SA/input/002/Equil_JT60_prova01_e_refined_COCOS7.eqdsk';cocos=7
#on personal pc
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
    eqd_fname='/home/matval/WORK/JT-60SA/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq';cocos=2
dir+='WORK/RMPs/'                                             
fncoils=dir+'B_RESP_3D_ITERlike'
fnout='bfield_ascot_mars_3D_iterlike.h5'
#eqd_fname='/home/vallar/WORK/JT-60SA/input/005/JT-60SA_scenario5_eqdsk'
#eqd_fname='JT-60SA_scenario2_highden_eqdsk_chease_cocos13_smoothed.geq'
#eqd_fname='JT-60SA_scenario5_eqdsk_chease_cocos02_smoothed.geq'

bnorm_calculation=0
rho_bnorm = 1
data, Bphi, BR, Bz, theta, new_phi, newBnorm, R, z, rho_2d, eq=\
m2a.read_fields(fncoils, fnout, eqd_fname, \
	bnorm_calculation, cocos, rho_bnorm)
