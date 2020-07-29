"""
calling the magnetic field rmpplot production
"""
import biosaw2ascot_EFCC as b2a   
import os
#on spcpc
dir='/home/vallar/'
eqd_fname='/home/vallar/WORK/JT-60SA/input/002/Equil_JT60_prova01_e_refined_COCOS7.eqdsk';cocos=7
#on personal pc
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
    eqd_fname='/home/matval/WORK/JT-60SA/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq';cocos=2
dir+='WORK/JT-60SA/3D/biosaw/efcc_output/'                                             
fncoils=dir+'wholecoils_EFCC540phi.h5'
fnout='./EFCC_UL_n3_0deg0deg120deg_003.h5'
#eqd_fname='/home/vallar/WORK/JT-60SA/input/005/JT-60SA_scenario5_eqdsk'
#eqd_fname='JT-60SA_scenario2_highden_eqdsk_chease_cocos13_smoothed.geq'
#eqd_fname='JT-60SA_scenario5_eqdsk_chease_cocos02_smoothed.geq'
nmode=3
U=1
L=1
M=0
phases=[0,0,120] #Careful! putting numbers close to 60 will produce a 0 field!
bnorm_calculation=0
rho_bnorm = 1
data, Bphi, BR, Bz, theta, new_phi, newBnorm, R, z, rho_2d, eq=\
b2a.produce_fields(fncoils, fnout, eqd_fname, \
	nmode,U, M, L,phases, bnorm_calculation, cocos, rho_bnorm)
