import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from lmfit import Model
from vasp_main import Xdatcar

def tanhfit(z,rho_cm,zw,delta):
    
    # Force conservation of mass
    rho_off = (rho_avg - 2*rho_cm*zw)/(1 - (2*zw))
    
    tanhp = np.tanh( ( ((z-zcm) - np.rint(z-zcm)) + zw) /delta)
    tanhm = np.tanh( ( ((z-zcm) - np.rint(z-zcm)) - zw) /delta)
    rho = rho_off + ((rho_cm - rho_off)/2. * (tanhp - tanhm))
    return rho

def rho_off_err(rho_avg, rho_cm, zw, rho_cm_err, zw_err):
    dfdrho_cm = (-2*zw)/(1 - (2*zw))
    dfdzw = (2*(rho_avg - rho_cm))/( (1-(2*zw))**2. )
    
    df = np.sqrt((rho_cm_err**2. * dfdrho_cm**2.) + (zw_err**2. * dfdzw**2.))
    return df

########################################################################

xdat = Xdatcar()

atom_colors = {"Fe": "papayawhip", "Mg": "forestgreen", "O": "red", "Si": "dodgerblue"}

########################################################################

e_star_data = np.loadtxt(fname = 'est_star_rho',unpack=True)
rho1_data = np.loadtxt(fname = './dens1',unpack=True)
rho2_data = np.loadtxt(fname = './dens2',unpack=True)
lambda_rho_data = np.loadtxt(fname = './lambda_star_rho',unpack=True)

with open('atoms_rho') as f:
    lines = f.readlines()

atoms = []
for i,x in enumerate(lines):
    atoms.append(lines[i].strip('\n').strip('   '))

zref = e_star_data[0,:]
for i in range(1,np.size(e_star_data[:,1])):
    plt.scatter(zref,e_star_data[i,:],edgecolor='black',
                color=atom_colors[xdat.atypes[i-1]])

plt.legend(atoms)
plt.xlabel('z ($\AA$)',size=13)
plt.ylabel(r'$\Delta \rho$(z)',size=13)
plt.show()

int_rho = np.zeros(np.size(e_star_data[:,1])-1)
rho0 = np.zeros_like(e_star_data[:,:])

for i in range(1,np.size(e_star_data[:,1])):
    int_rho[i-1] += np.trapz(e_star_data[i,:],zref)
    rho0[i-1,:] += e_star_data[i,:] - (int_rho[i-1]/xdat.cell[2,2]) + (xdat.nelem[i-1]/xdat.volume)
    plt.scatter(zref,rho0[i-1,:],edgecolor='black',\
                color=atom_colors[xdat.atypes[i-1]])

plt.legend(atoms)
plt.xlabel('z ($\AA$)',size=13)
plt.ylabel(r'$\rho$(z)',size=13)
plt.show()

for i in range(1,np.size(e_star_data[:,1])):
    if i == 1:
        plt.plot(zref, lambda_rho_data[i,:],color="black")
    else:
        plt.plot(zref, lambda_rho_data[i,:],color=atom_colors[xdat.atypes[i-1]])

plt.xlabel('z ($\AA$)',size=13)
plt.ylabel('$\lambda^{*}$(z)',size=13)
plt.legend(atoms)
plt.show()

########################################################################

# gmodel = Model(tanhfit)
# print('paramter namese: {}' .format(gmodel.param_names))
# print('independent variables: {}' .format(gmodel.independent_vars))

# # delta=0.1

# zcm = 0.3
# rho_avg = xdat.nelem[0]/xdat.volume
# result = gmodel.fit(e_star_data[1,:]+rho0[0],z=zref/xdat.cell[2,2],rho_cm=0.1,zw=0.25,delta=0.1)
# # print(result.fit_report())
# # print("########################################################################")

# fe5500_cm = result.params['rho_cm'].value
# fe5500_cmerr = result.params['rho_cm'].stderr
# fe5500_off = ((e_star_data[1,:]+rho0[0]) - (2*fe5500_cm*result.params['zw'].value))/(1 - (2*result.params['zw'].value))
# fe5500_offerr = rho_off_err((e_star_data[1,:]+rho0[0]), fe5500_cm, result.params['zw'].value, fe5500_cmerr, result.params['zw'].stderr)

# ########################################################################

# zcm = 0.69
# rho_avg = rho0_mg
# result2 = gmodel.fit(e_star_mg+rho0_mg,z=zref_mg/L,rho_cm=16/64,zw=0.25)
# # print(result2.fit_report())
# # print("########################################################################")

# mg5500_cm = result2.params['rho_cm'].value
# mg5500_cmerr = result2.params['rho_cm'].stderr
# mg5500_off = (rho0_mg - (2*mg5500_cm*result2.params['zw'].value))/(1 - (2*result2.params['zw'].value))
# mg5500_offerr = rho_off_err(rho0_mg, mg5500_cm, result2.params['zw'].value, mg5500_cmerr, result2.params['zw'].stderr)

# ########################################################################

# zcm = 0.7
# rho_avg = rho0_o
# result3 = gmodel.fit(e_star_o+rho0_o,z=zref_o/L,rho_cm=16/64,zw=0.25)
# # print(result3.fit_report())
# # print("########################################################################")

# o5500_cm = result3.params['rho_cm'].value
# o5500_cmerr = result3.params['rho_cm'].stderr
# o5500_off = (rho0_o - (2*o5500_cm*result3.params['zw'].value))/(1 - (2*result3.params['zw'].value))
# o5500_offerr = rho_off_err(rho0_o, o5500_cm, result3.params['zw'].value, o5500_cmerr, result3.params['zw'].stderr)

# ########################################################################

# zcm = 0.88
# rho_avg = rho0_si
# result4 = gmodel.fit((e_star_si+rho0_si),z=zref_si/L,rho_cm=16/64,zw=0.25)
# # print(result4.fit_report())
# # print("########################################################################")

# si5500_cm = result4.params['rho_cm'].value
# si5500_cmerr = result4.params['rho_cm'].stderr
# si5500_off = (rho0_si - (2*si5500_cm*result4.params['zw'].value))/(1 - (2*result4.params['zw'].value))
# si5500_offerr = rho_off_err(rho0_si, si5500_cm, result4.params['zw'].value, si5500_cmerr, result4.params['zw'].stderr)

# ########################################################################

# # plt.plot(zref_fe, result.best_fit, ':', color='black', label = 'Fe best fit')
# # plt.plot(zref_mg, result2.best_fit, ':', color='forestgreen', label = 'Mg best fit')
# # plt.plot(zref_o, result3.best_fit, ':', color='red', label = 'O best fit')
# # plt.plot(zref_si, result4.best_fit, ':', color='dodgerblue', label = 'Si best fit')
