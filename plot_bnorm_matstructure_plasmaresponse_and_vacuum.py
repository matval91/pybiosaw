"""
Script to plot Bnorm at different rho
"""
import os
import scipy.io as io
import numpy as np
import matplotlib.pyplot as plt
import utils.plot_utils as pu
import a4py.classes.ReadEQDSK as ReadEQDSK
pu.common_style()
dir ='/home/vallar/'
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
    dir+='WORK/RMPs/plasma_response/'

#Read Bnorm PLASMA RESPONSE
rho=1.0
dict_1=io.loadmat(f'{dir}Bnorm_plasmaresponse_1.mat')
rho=0.8
dict_0d8=io.loadmat(f'{dir}Bnorm_plasmaresponse_{rho}.mat')
rho=0.5
dict_0d5=io.loadmat(f'{dir}Bnorm_plasmaresponse_{rho}.mat')
rho=0.6
dict_0d6=io.loadmat(f'{dir}Bnorm_plasmaresponse_{rho}.mat')
rho=0.7
dict_0d7=io.loadmat(f'{dir}Bnorm_plasmaresponse_{rho}.mat')
rho=0.1
dict_0d1=io.loadmat(f'{dir}Bnorm_plasmaresponse_{rho}.mat')
rho=0.3
dict_0d3=io.loadmat(f'{dir}Bnorm_plasmaresponse_{rho}.mat')
rho=0.9
dict_0d9=io.loadmat(f'{dir}Bnorm_plasmaresponse_{rho}.mat')
# Last time, 0.2 and 0.4 went wrong
rho=0.4
dict_0d4=io.loadmat(f'{dir}Bnorm_plasmaresponse_{rho}.mat')
rho=0.2
dict_0d2=io.loadmat(f'{dir}Bnorm_plasmaresponse_{rho}.mat')
# filling array
x_rho=np.array([])
y_bnorm=np.array([])
for dd in [dict_1, dict_0d5, dict_0d8, dict_0d7, dict_0d8, dict_0d1,  dict_0d4, dict_0d9]:
	max = np.max(dd['Bnorm'])*1e4
	min = np.min(dd['Bnorm'])*1e4
	pplot=min
	if max>abs(min):
		pplot=max;
#	ax.scatter(dd['rho'], pplot, color='k')
	x_rho=np.append(x_rho, dd['rho'])
	y_bnorm = np.append(y_bnorm, pplot)
x_plasmaresponse=x_rho[np.argsort(x_rho)]
y_plasmaresponse=np.abs(y_bnorm[np.argsort(x_rho)])

dir ='/home/vallar/'
eqd_fname='/home/vallar/WORK/JT-60SA/input/002/Equil_JT60_prova01_e_refined_COCOS7.eqdsk';cocos=7
if os.uname().nodename!='spcpc182':
    dir='/home/matval/'
    dir+='WORK/RMPs/vacuum/'
    eqd_fname='/home/matval/WORK/JT-60SA/JT-60SA_scenario2_highden_eqdsk_chease_cocos02_smoothed.geq';cocos=2

#Read VACUUM
rho=1.0
dict_1=io.loadmat(f'{dir}Bnorm_1.mat')
rho=0.8
dict_0d8=io.loadmat(f'{dir}Bnorm_{rho}.mat')
rho=0.5
dict_0d5=io.loadmat(f'{dir}Bnorm_{rho}.mat')
rho=0.2
dict_0d2=io.loadmat(f'{dir}Bnorm_{rho}.mat')
rho=0.6
dict_0d6=io.loadmat(f'{dir}Bnorm_{rho}.mat')
rho=0.7
dict_0d7=io.loadmat(f'{dir}Bnorm_{rho}.mat')
rho=0.1
dict_0d1=io.loadmat(f'{dir}Bnorm_{rho}.mat')
x_rho=np.array([])
y_bnorm=np.array([])
for dd in [dict_1, dict_0d5, dict_0d8, dict_0d2, dict_0d7, dict_0d8, dict_0d1]:
	max = np.max(dd['Bnorm'])*1e4
	min = np.min(dd['Bnorm'])*1e4
	pplot=min
	if max>abs(min):
		pplot=max;
#	ax.scatter(dd['rho'], pplot, color='k')
	x_rho=np.append(x_rho, dd['rho'])
	y_bnorm = np.append(y_bnorm, pplot)
y_vacuum=np.abs(y_bnorm[np.argsort(x_rho)])
x_vacuum=x_rho[np.argsort(x_rho)]

#Read eqdsk
eq=ReadEQDSK.ReadEQDSK(eqd_fname)
# plot Bnorm
f=plt.figure();
ax=f.add_subplot(111)
ax.plot(x_vacuum, y_vacuum, 'kx-', label='Vacuum')
ax.plot(x_plasmaresponse, y_plasmaresponse, 'ko-', label='With plasma response')
ax2 = ax.twinx()
ax2.plot(eq.rhopsi, eq.q, 'r')
ax.set_xlabel(r'$\rho_{POL}$')
ax.set_ylabel(r'max(B$_{norm}$) [$10^{-4}$ T]')
ax2.set_ylabel(r'q', color='r')
ax.grid('on')
ax.legend(loc='best')
plt.show()
f.tight_layout()