import numpy as np
import netCDF4
import matplotlib.pyplot as plt

#import wall geometry to plot over
wall_file = '../setup/assets/geom-SASV/gitr_rz.txt'
with open(wall_file, 'r') as file:
    wall = file.readlines()

R = np.zeros(len(wall))
Z = np.zeros(len(wall))
for i,line in enumerate(wall):
    point = line.split()
    R[i] = float(point[0])
    Z[i] = float(point[1])


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






if __name__ == "__main__":
    plot_history2D()
