import numpy as np
import netCDF4
import matplotlib.pyplot as plt

#import wall geometry to plot over
wall_file = '../setup/assets/geom-SASV6/gitr_rz.txt'
with open(wall_file, 'r') as file:
    wall = file.readlines()

R = np.zeros(len(wall))
Z = np.zeros(len(wall))
for i,line in enumerate(wall):
    point = line.split()
    R[i] = float(point[0])
    Z[i] = float(point[1])

#import W surface indices
surf_file = '../setup/assets/geom-SASV6/surf_ind.txt'
with open(surf_file, 'r') as file:
    surf = file.readlines()

Rsurf = np.zeros(len(surf))
Zsurf = np.zeros(len(surf))
for i,v in enumerate(surf):
    Rsurf[i] = R[i]
    Zsurf[i] = Z[i]

dist = np.sqrt((Rsurf[:-1]-Rsurf[1:])**2 + (Zsurf[:-1]-Zsurf[1:])**2)
rmrs = np.append(np.zeros(1), np.cumsum(dist))

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

def plot_surf(ls,fs):
    surface = netCDF4.Dataset("surface.nc", "r", format="NETCDF4")

    grossEro = surface.variables['grossErosion'][:]
    grossDep = surface.variables['grossDeposition'][:]
    netEro = grossEro-grossDep

    plt.close()
    plt.plot(rmrs,grossEro,'r', label='Gross Erosion',linewidth=ls)
    plt.plot(rmrs,grossDep,'g', label='Redeposition',linewidth=ls)
    plt.plot(rmrs,netEro,'k', label='Net Erosion',linewidth=ls)
    plt.xlabel('r-rsep [m]',fontsize=fs)
    plt.ylabel('Particles per Second',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend(loc='upper left', fontsize=fs)
    plt.title('GITR Predicted Erosion and Redeposition Profiles',fontsize=fs)
    #plt.savefig('plots/surface.pdf')






if __name__ == "__main__":
    #plot_history2D()
    plot_surf(5,30)
