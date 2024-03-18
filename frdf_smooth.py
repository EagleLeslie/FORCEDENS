import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

data = np.loadtxt(fname = 'rad1')
data2 = np.loadtxt(fname = 'rad2')
data3 = np.loadtxt(fname = 'rad3')
lambda_data_gr = np.loadtxt(fname = 'lambda_star_gr')

with open('atoms') as f:
    lines = f.readlines()

atoms = []
for i,x in enumerate(lines):
    atoms.append(lines[i].strip('\n').strip('   '))

r = data[:,0]
for i in range(1,np.size(data[0,:])):
    plt.plot(r,data[:,i],color=cm.turbo(i/np.size(data[0,:])))

r2 = data2[:,0]
for i in range(1,np.size(data2[0,:])):
    plt.plot(r2,data2[:,i],color=cm.turbo(i/np.size(data2[0,:])),linestyle=':')

r3 = data3[:,0]
for i in range(1,np.size(data2[0,:])):
    plt.plot(r3,data3[:,i],color=cm.turbo(i/np.size(data3[0,:])),linestyle='--')

plt.legend(atoms)
plt.xlabel('Distance r, ($\AA$)',size=13)
plt.ylabel('g(r)',size=13)
plt.show()

plt.plot(r,data[:,1],label='g$_{\infty}$(r)')
plt.plot(r2,data2[:,1],label='g$_{0}$(r)')
plt.plot(r3,data3[:,1],label='g$_{\lambda}$(r)')

plt.legend()
plt.xlabel('Distance r, ($\AA$)',size=13)
plt.ylabel('g(r)',size=13)
plt.title('Mg-Mg')
plt.show()

r4 = lambda_data_gr[:,0]
for i in range(1,np.size(lambda_data_gr[0,:])):
    plt.plot(r4,lambda_data_gr[:,i],color=cm.turbo(i/np.size(lambda_data_gr[0,:])))

plt.legend(atoms)
plt.xlabel('Distance r, ($\AA$)',size=13)
plt.ylabel('$\lambda^{*}$(r)',size=13)
plt.show()

