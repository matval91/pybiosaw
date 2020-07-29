"""
Create the bnorm at different radii
"""
import ascotinput_from_biosawoutput_EFCCcoils_plasmaresponse as mf   
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
fnout=''

bnorm_calculation=1
rho_bnorm=1
data, Bphi, BR, Bz, theta, new_phi, newBnorm, R_rho2d, z_rho2d, rho_2d, eq=\
mf.produce_fields(fncoils, fnout, eqd_fname, \
	bnorm_calculation, cocos, rho_bnorm)
for rho_bnorm in np.linspace(0.1, 0.9, 9):
	try:
		mf.Bnormal(data, BR, Bz, eq, rho_bnorm, R_rho2d, z_rho2d, rho_2d)
	except:
		print(f'Continuing for rho {rho_bnorm}')