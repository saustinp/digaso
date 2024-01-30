import numpy as np
import matplotlib.pyplot as plt

fig, ax1 = plt.subplots(1,)
fig, ax2 = plt.subplots(1,)

# 2 ns -- 1 indexed due to IC
arry = np.loadtxt('comp_1e13_f_line_000003.txt', skiprows=1)
# x y electron pos_ion phi electric_fld
x = arry[:,1]
ne = arry[:,2]
E = arry[:,5]
ax1.semilogy(x,ne,'k')
ax2.plot(x,E*3e6,'k')

arry = np.loadtxt('ns2.txt', delimiter=',')
# x y electron pos_ion phi electric_fld
x = arry[:,0]*1e-2
ne = arry[:,1]*1e12
E = arry[:,2]*3e6
ax1.semilogy(x,ne,'r--')
ax2.plot(x,E,'r--')

# 4 ns
arry = np.loadtxt('comp_1e13_f_line_000005.txt', skiprows=1)
# x y electron pos_ion phi electric_fld
x = arry[:,1]
ne = arry[:,2]
E = arry[:,5]
ax1.semilogy(x,ne,'k')
ax2.plot(x,E*3e6,'k')

arry = np.loadtxt('ns4.txt', delimiter=',')
# x y electron pos_ion phi electric_fld
x = arry[:,0]*1e-2
ne = arry[:,1]*1e12
E = arry[:,2]*3e6
ax1.semilogy(x,ne,'r--')
ax2.plot(x,E,'r--')

# 6 ns
arry = np.loadtxt('comp_1e13_f_line_000007.txt', skiprows=1)
# x y electron pos_ion phi electric_fld
x = arry[:,1]
ne = arry[:,2]
E = arry[:,5]
ax1.semilogy(x,ne,'k')
ax2.plot(x,E*3e6,'k')

arry = np.loadtxt('ns6.txt', delimiter=',')
# x y electron pos_ion phi electric_fld
x = arry[:,0]*1e-2
ne = arry[:,1]*1e12
E = arry[:,2]*3e6
ax1.semilogy(x,ne,'r--')
ax2.plot(x,E,'r--')

# 7 ns
arry = np.loadtxt('comp_1e13_f_line_000008.txt', skiprows=1)
# x y electron pos_ion phi electric_fld
x = arry[:,1]
ne = arry[:,2]
E = arry[:,5]
ax1.semilogy(x,ne,'k')
ax2.plot(x,E*3e6,'k')

arry = np.loadtxt('ns7.txt', delimiter=',')
# x y electron pos_ion phi electric_fld
x = arry[:,0]*1e-2
ne = arry[:,1]*1e12
E = arry[:,2]*3e6
ax1.semilogy(x,ne,'r--')
ax2.plot(x,E,'r--')

ax1.set_title('Electron number density, on-axis, particles/m^3')
ax2.set_title('|E|, on-axis, V/m')

plt.show()

