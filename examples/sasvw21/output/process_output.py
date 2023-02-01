import numpy as np
import netCDF4
import matplotlib.pyplot as plt

def init():
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
    
    #set plotting style defaults
    plt.rcParams.update({'font.size':12})
    plt.rcParams.update({'lines.linewidth':0.5})
    plt.rcParams.update({'lines.markersize':1})

    return R,Z,Rsurf,Zsurf,area,rmrs
R,Z,Rsurf,Zsurf,area,rmrs = init()

def plot_history2D(basic=0, continuousChargeState=1, endChargeState=0):
    history = netCDF4.Dataset("history.nc", "r", format="NETCDF4")

    nP = len(history.dimensions['nP'])
    nT = len(history.dimensions['nT'])
    x = history.variables['x']
    z = history.variables['z']
    charge = history.variables['charge']

    plt.close()
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.axis('scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Particle Tracks')
    
    #define charge state to color mapping
    colors = {0:'black', 1:'red', 2:'orange', 3:'olive', 4:'green', 5:'cyan', \
              6:'purple', 7:'darkmagenta', 8:'pink', 9:'deep pink', 10:'gray'}
    
    # all particle source vars ordered as (nP, nT)
    if basic==1:
        for p in range(0,nP):
            plt.plot(x[p][:],z[p][:])
    
    if continuousChargeState==1:
        for p in range(0,nP):
            print(p,'out of', nP)
            for t in range(0,nT-1):
                plt.plot(x[p][t:t+2],z[p][t:t+2], colors[charge[p][t]])

    if endChargeState==1:
        for p in range(0,nP):
            plt.plot(x[p][:],z[p][:], colors[charge[p][-1]])
        
    plt.xlim(1.43, 1.56)
    plt.ylim(1.07, 1.25)
    plt.savefig('plots/history.pdf')

    return

def plot_surf_nc(ls,fs, pps_per_nP):
    surface = netCDF4.Dataset("surface_nP1e5_nT1e6.nc", "r", format="NETCDF4")
    #print(surface.variables['grossErosion'][:])

    grossEro = (surface.variables['grossErosion'][:])
    grossDep = (surface.variables['grossDeposition'][:])
    netEro = (grossEro-grossDep)
    
    print('total gross eroded comp particles',sum(grossEro))
    print('total redeposited comp particles',sum(grossDep))
    print('total net eroded comp particles',sum(netEro))
    
    grossEro = np.average([grossEro[:-1], grossEro[1:]],axis=0)*pps_per_nP/area
    grossDep = np.average([grossDep[:-1], grossDep[1:]],axis=0)*pps_per_nP/area
    netEro = np.average([netEro[:-1], netEro[1:]],axis=0)/area
    
    print('rmrs length',len(rmrs))
    print('surf length',len(grossEro))
    print('total gross eroded flux',sum(grossEro))
    print('total redeposited flux',sum(grossDep))
    print('total net eroded flux',sum(netEro))
    
    grossEro_cumsum = np.cumsum(grossEro)
    grossDep_cumsum = np.cumsum(grossDep)
    netEro_cumsum = np.cumsum(netEro)
    
    plt.close()
    plt.plot(rmrs,grossEro,'r', label='Gross Erosion')
    plt.plot(rmrs,grossDep,'g', label='Redeposition')
    plt.plot(rmrs,netEro,'k', label='Net Erosion')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [#/m2s]')
    plt.legend(loc='upper left')
    plt.title('GITR Predicted Erosion and Redeposition Profiles')
    plt.savefig('plots/surface.pdf')
    
    plt.close()
    plt.plot(rmrs,grossEro_cumsum,'r', label='Gross Erosion')
    plt.plot(rmrs,grossDep_cumsum,'g', label='Redeposition')
    plt.plot(rmrs,netEro_cumsum,'k', label='Net Erosion')
    plt.yscale('log')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [#/m2s]')
    plt.legend(loc='upper left')
    plt.title('GITR Predicted Cumulative Sum of\nErosion and Redeposition Profiles')
    plt.savefig('plots/surface_cumsum.pdf')

def plot_surf_plasma_params(plot_variables):
    profilesFile = '../input/profiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    r_mesh = profiles.variables['r'][:]
    z_mesh = profiles.variables['z'][:]
    
    r_mid = np.average(np.array([Rsurf[:-1],Rsurf[1:]]),axis=0)
    z_mid = np.average(np.array([Zsurf[:-1],Zsurf[1:]]),axis=0)
    r_indices = interpolate(r_mid,r_mesh)
    z_indices = interpolate(z_mid,z_mesh)
    
    #scoot the r_indices 1 cell to the left if the profiles.nc gridcell is too far off
    br = profiles.variables['br'][:][z_indices,r_indices]
    for i in range(len(br)):
        if br[i]==-1: r_indices[i] -= 1
    
    ne = profiles.variables['ne'][:][z_indices, r_indices]
    te = profiles.variables['te'][:][z_indices, r_indices]
    
    if plot_variables==1:            
        #electron density along surface
        plt.close()
        plt.plot(rmrs, ne)
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('ne [m-3]')
        plt.title('Electron Density along the SAS-VW Divertor')
        plt.savefig('plots/surf_ne.png')
    
        #electron temperature along surface
        plt.close()
        plt.plot(rmrs, te)
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Te [eV]')
        plt.title('Electron Temperature along the SAS-VW Divertor')
        plt.savefig('plots/surf_te.png')
        
        #electron temperature dependence of ionization rate
        plt.close()
        kte = np.arange(1,35,0.01)
        x0 = 7.86
        x1 = 16.37
        x2 = 26.0
        x3 = 38.2
        x4 = 51.6
        invtau0 = 1/(np.exp(x0/kte)/np.sqrt(kte))
        invtau1 = 1/(np.exp(x1/kte)/np.sqrt(kte))
        invtau2 = 1/(np.exp(x2/kte)/np.sqrt(kte))
        invtau3 = 1/(np.exp(x3/kte)/np.sqrt(kte))
        invtau4 = 1/(np.exp(x4/kte)/np.sqrt(kte))
        plt.close()
        plt.plot(kte, invtau0, 'red', label='0 → 1')
        plt.plot(kte, invtau1, 'darkorange', label='1 → 2')
        plt.plot(kte, invtau2, 'gold', label='2 → 3')
        plt.plot(kte, invtau3, 'green', label='3 → 4')
        plt.plot(kte, invtau4, 'blue', label='4 → 5')
        plt.xlabel('Te [eV]')
        plt.ylabel('Relative Ionization Rate [1/s]')
        plt.title('Te Dependence of Ionization Rate for W')
        plt.legend()
        plt.savefig('plots/IonizRate(Te)')
        
        #ionization rate along surface
        surf_IonizRate0 = ne/(np.exp(x0/te)/np.sqrt(te))
        surf_IonizRate1 = ne/(np.exp(x1/te)/np.sqrt(te))
        surf_IonizRate2 = ne/(np.exp(x2/te)/np.sqrt(te))
        surf_IonizRate3 = ne/(np.exp(x3/te)/np.sqrt(te))
        surf_IonizRate4 = ne/(np.exp(x4/te)/np.sqrt(te))
        plt.close()
        plt.plot(rmrs, surf_IonizRate0, 'red', label='0 → 1')
        plt.plot(rmrs, surf_IonizRate1, 'darkorange', label='1 → 2')
        plt.plot(rmrs, surf_IonizRate2, 'gold', label='2 → 3')
        plt.plot(rmrs, surf_IonizRate3, 'green', label='3 → 4')
        plt.plot(rmrs, surf_IonizRate4, 'blue', label='4 → 5')
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Relative Ionization Rate [1/s]')
        plt.title('Ionization Rate along SAS-VW')
        plt.legend()
        plt.savefig('plots/surf_IonizRate')



def interpolate(small,big):
    indices = np.zeros(len(small))
    for i in range(len(small)):
        diff = np.min(np.abs(big-small[i]))
        index_possibilities = np.array([np.where(big==small[i]+diff)[0], \
                             np.where(big==small[i]-diff)[0]], dtype=object)
        try: indices[i] = index_possibilities[1]
        except: indices[i] = index_possibilities[0]

    return indices.astype(int)

def plot(r,z,line,scatter, \
             color,legend,labelname, \
             ylabel,title,savename):
    if line == 1:
        plt.plot(r, z, label=labelname, color=color)
    else:
        return
    if scatter == 1:    
        plt.scatter(r, z,color=color)
    else:
        return
    if legend == 1:
        plt.legend()
    plt.axis('scaled')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel(ylabel)
    plt.title(title)
    savename0 = 'plots/'+savename+'.pdf'
    plt.savefig(savename0)



if __name__ == "__main__":
    plot_history2D()
    #plot_surf_nc(1,12,111879178639.80714)
    #plot_surf_plasma_params(1)
