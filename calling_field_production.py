"""
calling the magnetic field rmpplot production
"""
import ascotinput_from_biosawoutput_EFCCcoils as mf                                                
fncoils='/home/vallar/WORK/JT-60SA/3D/RMPs/ITER-like RMPs/wholecoils_EFCC540phi.h5'
fnout='./EFCC_UL_n3_0deg0deg120deg_005.h5'
fnout=''
eqd_fname='/home/vallar/WORK/JT-60SA/input/002/Equil_JT60_prova01_e_refined_COCOS7.eqdsk'
#eqd_fname='/home/vallar/WORK/JT-60SA/input/005/JT-60SA_scenario5_eqdsk'
#eqd_fname='JT-60SA_scenario2_highden_eqdsk_chease_cocos13_smoothed.geq'
#eqd_fname='JT-60SA_scenario5_eqdsk_chease_cocos02_smoothed.geq'
cocos=7
nmode=3
U=1
L=1
M=0
phases=[0,0,120] #Careful! putting numbers close to 60 will produce a 0 field!
bnorm_calculation=1
data, Bphi, BR, Bz, theta, new_phi, newBnorm=\
mf.produce_fields(fncoils, fnout, eqd_fname, \
	nmode,U, M, L,phases, bnorm_calculation, cocos)
