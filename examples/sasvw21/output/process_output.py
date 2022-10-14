import numpy as np
import netCDF4
import matplotlib.pyplot as plt

#import wall geometry to plot over
wall_file = '../setup/assets/gitr_rz.txt'
with open(wall_file, 'r') as file:
    wall = file.readlines()

R = np.zeros(len(wall))
Z = np.zeros(len(wall))
for i,line in enumerate(wall):
    point = line.split()
    R[i] = float(point[0])
    Z[i] = float(point[1])

#import W surface indices
surf_file = '../setup/assets/surf_ind.txt'
with open(surf_file, 'r') as file:
    surf = file.readlines()

surf = np.array(surf,dtype='int')
Rsurf = np.zeros(len(surf))
Zsurf = np.zeros(len(surf))
for i,v in enumerate(surf):
    Rsurf[i] = R[v]
    Zsurf[i] = Z[v]

area = 0.12*np.sqrt(np.power((Rsurf[:-1]-Rsurf[1:]),2) + np.power((Zsurf[:-1]-Zsurf[1:]),2))
dist = np.sqrt((Rsurf[:-1]-Rsurf[1:])**2 + (Zsurf[:-1]-Zsurf[1:])**2)
rmrs =  np.cumsum(dist)

def plot_history2D():
    history = netCDF4.Dataset("history.nc", "r", format="NETCDF4")

    x = history.variables['x']
    z = history.variables['z']

    plt.close()
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.axis('scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Particle Tracks')

    for i in range(0,len(x)):
        plt.plot(x[i][:],z[i][:], linewidth=0.2)

    plt.xlim(1.43, 1.56)
    plt.ylim(1.07, 1.25)
    plt.savefig('plots/history.pdf')

def plot_surf(ls,fs, pps_per_nP):
    surface = netCDF4.Dataset("surface.nc", "r", format="NETCDF4")

    grossEro = (surface.variables['grossErosion'][:])*pps_per_nP
    grossDep = (surface.variables['grossDeposition'][:])*pps_per_nP
    netEro = (grossEro-grossDep)
    
    grossEro = np.average([grossEro[:-1], grossEro[1:]],axis=0)/area
    grossDep = np.average([grossDep[:-1], grossDep[1:]],axis=0)/area
    netEro = np.average([netEro[:-1], netEro[1:]],axis=0)/area
    
    print('rmrs',len(rmrs))
    print('surf',len(grossEro))

    plt.close()
    plt.plot(rmrs,grossEro,'r', label='Gross Erosion',linewidth=ls)
    plt.plot(rmrs,grossDep,'g', label='Redeposition',linewidth=ls)
    plt.plot(rmrs,netEro,'k', label='Net Erosion',linewidth=ls)
    plt.xlabel('D-Dsep [m]',fontsize=fs)
    plt.ylabel('Flux [#/m2s]',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(loc='upper left', fontsize=fs)
    plt.title('GITR Predicted Erosion and Redeposition Profiles',fontsize=fs)
    plt.savefig('plots/surface.pdf')






if __name__ == "__main__":
    #plot_history2D()
    plot_surf(1,30,111879178639.80714)
