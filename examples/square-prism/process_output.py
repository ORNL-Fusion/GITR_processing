import sys, os
import numpy as np
from scipy import special
import scipy.interpolate as scii
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import netCDF4

################################################
# setting directories and special constants
################################################

run_directory = '.'
history_file = 'output/history.nc'
surface_file = 'output/surface.nc'
positions_file = 'output/positions.nc'
spec_file = 'output/spec.nc'

x_wall = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 0.0])
z_wall = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0.0])
x1 = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0])
z1 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0])
x2 = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, 0.0])
z2 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 0.0])

length = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 8.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 8.0])
stretched_surface = np.zeros(len(length))
stretched_surface[0] = 0.5*length[0]
for i in range(len(length)-1):
    i+=1
    stretched_surface[i] = stretched_surface[i-1] + 0.5*length[i-1] + 0.5*length[i]

#surface = np.array([0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0])
surface = np.ones(len(x1))
W_indices = np.where(surface==1)

################################################
# function definitions
################################################

def plot_history2D(basic=1, continuousChargeState=0, endChargeState=0, \
                   cylsymm=0, plot_particle_source=1, markersize=0):
        
    if plot_particle_source:
        particleSource = netCDF4.Dataset(run_directory+"/input/particleSource.nc", "r", format="NETCDF4")
        x0 = particleSource.variables['x'][:]
        z0 = particleSource.variables['z'][:]
    
    history = netCDF4.Dataset(history_file, "r", format="NETCDF4")

    plt.rcParams.update({'lines.linewidth':0.3})
    plt.rcParams.update({'lines.markersize':markersize})
    plt.rcParams.update({'font.size':14})

    nP = len(history.dimensions['nP'])
    print('nP:',nP)
    nT = len(history.dimensions['nT'])
    x = history.variables['x'][:]
    y = history.variables['y'][:]
    z = history.variables['z'][:]
    charge = history.variables['charge'][:]
    if cylsymm: 
        r = np.sqrt(x**2 + y**2)
    else: 
        r = x

    plt.close()
    if plot_particle_source: plt.scatter(x0,z0,marker='o',s=10)
    plt.axis('scaled')
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
        
    #define charge state to color mapping
    colors = {0:'black', 1:'firebrick', 2:'darkorange', 3:'gold', 4:'limegreen', 5:'dodgerblue', \
              6:'mediumpurple', 7:'darkviolet', 8:'darkmagenta', 9:'deeppink', 10:'gray', 11:'gray', \
                  12:'gray', 13:'gray', 14:'gray', 15:'gray', 16:'gray', 17:'gray', 18:'gray', 19:'gray', \
                      20:'gray', 21:'gray', 22:'gray', 23:'gray', 24:'gray', 25:'gray', 26:'gray', 27:'gray'}
    
    # all particle source vars ordered as (nP, nT)
    if basic==1:
        for p in range(0,nP):
            #if z[p][0]<=Z[W_indices][3] and z[p][0]>Z[W_indices][4]:
                print(p,'out of', nP)
                plt.plot(r[p][:],z[p][:])
                #plt.scatter(r[p][:],z[p][:],marker='o',s=5,c='b')
    
    counter=0
    if continuousChargeState==1:
        #for p in [671, 860, 1043, 1275, 1366, 1479]:  # particles that leak out of the divertor space in history_H.nc from pa-fav (Case #1)
        for p in np.arange(0,nP,1): # plot all particles by default
            t=0
            counter+=1
            while t<nT-1:
                if r[p][t] != r[p][t+1]: 
                    print("particle #", p, "moved at timestep", t, "with charge", charge[p][t])
                    plt.plot(r[p][t:t+2],z[p][t:t+2], colors[np.round(charge[p][t])])
                t+=1
    print('total particles:',counter)

    if endChargeState==1:
        for p in range(0,nP):
            plt.plot(r[p][:],z[p][:], colors[charge[p][-1]])
        plt.title('Particle Tracks by End Charge State')
        
    legend_dict = {'+0':'black', '+1':'firebrick', '+2':'darkorange', '+3':'gold', '+4':'limegreen', '+5':'dodgerblue', \
              '+6':'mediumpurple', '+7':'darkviolet', '+8':'darkmagenta', '+9':'deeppink', '+10':'gray', '+11':'gray', \
                  '+12':'gray', '+13':'gray', '+14':'gray', '+15':'gray', '+16':'gray', '+17':'gray', '+18':'gray', \
                      '+19':'gray', '+20':'gray', '+21':'gray', '+22':'gray', '+23':'gray', '+24':'gray', '+25':'gray', \
                          '+26':'gray', '+27':'gray'}
    
    patchList = []
    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    #if basic==0: plt.legend(handles=patchList, fontsize=8, loc=2) #upper-left=2, lower-left=3
    
    #slightly larger than the box
    plt.xlim(-0.5, 8.5)
    plt.ylim(-0.5, 8.5)
    
    plt.title('particle tracks from history.nc', fontsize=24)
    plt.savefig(run_directory+'/output/plots/history.svg')
    plt.close()
    return

def plot_surf_nc(nP=2e2, dt=1e-8, nT=1e4, cylsymm=0, pps_per_nP=1, plot_blocker=True):
        
    #calculate area from wall
    dist = np.sqrt(np.power(x1-x2,2) + np.power(z1-z2,2))
    if cylsymm:
        area = np.pi*(x1+x2)*dist # conical frustum surface area
    else:
        #area = (x2[W_indices]-x1[W_indices])**2
        area = length
    
    surface = netCDF4.Dataset(surface_file, "r", format="NETCDF4")
    positions = netCDF4.Dataset(positions_file, "r", format="NETCDF4")
    flightAngle = positions.variables['angle'][:]
    is_prompt_redep = flightAngle < 2*np.pi
    prompt_redep_rate = np.sum(is_prompt_redep) / len(flightAngle)

    grossEro = (surface.variables['grossErosion'][:])
    print(grossEro)
    grossDep = (surface.variables['grossDeposition'][:])
    print(grossDep)
    netDep = (grossDep-grossEro)
    
    print('\n')
    print('total gross eroded comp particles:',sum(grossEro))
    print('total redeposited comp particles:',sum(grossDep))
    print('total net deposited comp particles:',sum(netDep))
    
    #grossEro = np.average([grossEro[:-1], grossEro[1:]],axis=0)*pps_per_nP/area
    #grossDep = np.average([grossDep[:-1], grossDep[1:]],axis=0)*pps_per_nP/area
    grossEro = grossEro*pps_per_nP/area
    grossDep = grossDep*pps_per_nP/area
    netDep = netDep*pps_per_nP/area
    
    print('\n')
    print('surf length',len(grossEro))
    print('\n')
    print('total gross eroded flux',sum(grossEro*area)/sum(area))
    print('total redeposited flux',sum(grossDep*area)/sum(area))
    print('total net deposited flux',sum(netDep*area)/sum(area))
    print('redeposition rate',100 * (sum(grossDep*area)/sum(area)) / (sum(grossEro*area)/sum(area)), '%')
    if positions_file != '': print('prompt redeposition rate', 100 * prompt_redep_rate, '%')
    
    #take gross erosion of a slice in the filterscope range    
    V_grossEro = grossEro[3:5]
    
    print('\n')
    print('gross erosion in synthetic View:', V_grossEro)
    
    plt.rcParams.update({'font.size':24})
    plt.rcParams.update({'lines.linewidth':5}) 
    
    #plot main surface plot with 3 views
    plt.close()
    #plt.axvspan(3, 5, color='lightsalmon', alpha=0.5)
    #plt.axvspan(8, 16, color='gray', alpha=0.5)
    #plt.axvspan(24, 32, color='gray', alpha=0.5)
    
    plt.plot(stretched_surface, grossEro,'r', label='Gross Erosion')
    plt.plot(stretched_surface, grossDep,'g', label='Gross Deposition')
    #plt.plot(stretched_surface, netDep,'k', label='Net Deposition')
    
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('\u0393$_{W}$ [m$^{-2}$s$^{-1}$]')
    plt.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
    plt.legend(loc='upper left')
    plt.title('1D Erosion and Deposition', fontsize=30)
    plt.savefig('output/plots/surface.png')
    
    return grossEro

def PEC_interpolation(te, ne):
    PEC_file = '../../../SAS-VW-Data/sxb/martin_pec_4009_only_5d_ion.dat'
    PEC_data = np.loadtxt(PEC_file, dtype='double')

    te_grid = PEC_data[0][1:]
    ne_grid = PEC_data.transpose()[0][1:]
    PEC_grid = PEC_data[1:,1:] # cm3 / s
    
    PEC_interpolated = scii.interpn((ne_grid,te_grid),PEC_grid,(ne,te))[0]
    
    return PEC_interpolated

def spec_line_integration(pps_per_nP=1, num_points=100, dt=1e-8):
    
    # set start and end points for the line along which we are integrating
    # all (r,z) points are in [cm] for this calculation
    r_source = 3.5
    z_source = 8
    r_targ = 3.5
    z_targ = 0
            
    # define the points of interest along the line from (r_source,z_source) to (r_targ,z_targ)
    r_line = np.linspace(r_source, r_targ, num_points+1)
    z_line = np.linspace(z_source, z_targ, num_points+1)
    
    # import neutral W density results from spec.nc
    spec = netCDF4.Dataset(spec_file, "r", format="NETCDF4")
    
    gridz_spec_original = spec.variables['gridZ'][:]#*100 #[cm]
    gridr_spec_original = spec.variables['gridR'][:]#*100 #[cm]
    dz = gridz_spec_original[1]-gridz_spec_original[0]
    dr = gridr_spec_original[1]-gridr_spec_original[0]
    gridz_spec = np.append(gridz_spec_original, gridz_spec_original[-1]+dz)
    gridr_spec = np.append(gridr_spec_original, gridr_spec_original[-1]+dr)
    #gridz_spec = gridz_spec_original
    #gridr_spec = gridr_spec_original
    
    rr, zz = np.meshgrid(gridz_spec,gridr_spec)
    rr1 = rr[:-1,:-1]
    rr2 = rr[1:,1:]
    zz1 = zz[:-1,:-1]
    zz2 = zz[1:,1:]
    
    dV = (zz2-zz1)*0.5*(rr2**2 - rr1**2) #*2*np.pi

    ni_unitless = spec.variables['n'][:] #unitless ion densities for all charge states
    n0_unitless = ni_unitless[0] #unitless neutral densities
    # gets 0 for view==1 but you can check that it's not broken because it's non-zero for ni_unitless[2]
    ni_2D = pps_per_nP * (n0_unitless/dV) * dt # W0 cm-3
    
    plt.close()
    plt.rcParams.update({'lines.linewidth':1})
    plt.plot(x_wall,z_wall,'k',label='wall')
    plt.pcolormesh(rr,zz,ni_2D)
    
    ni_line = scii.interpn((gridz_spec_original,gridr_spec_original),ni_2D,(z_line,r_line))
    
    # interpolate plasma parameters (te,ne,ni) at each point along the line from profiles.nc
    te = 10
    ne = 1e19/1e6

    # interpolate PEC as a function of local ne and te
    PEC_line = PEC_interpolation(te, ne) #[ph cm3 s-1]
    
    # perform the integration per equation 2 of Bogen 1984
    dr = (r_targ-r_source)/num_points
    dz = (z_targ-z_source)/num_points
    dx = np.sqrt(dr**2+dz**2)
    
    epsilon = (1/(4*np.pi)) * ne * ni_line * PEC_line #[ph cm-3 s-1 str-1]
    intensity = np.sum(epsilon) * dx #[ph cm-2 s-1]
    print("intensity:", intensity)
    
    SXB = 30
    flux = 4*np.pi*intensity*SXB
    print("flux:", flux)
    
    plt.close()
    plt.rcParams.update({'lines.linewidth':1})
    plt.plot(x_wall,z_wall,'k',label='wall')
    plt.pcolormesh(rr,zz,ni_2D)
    #plt.vlines([2.5,3.5,4.5,5.5], 0, 8, 'm',linewidth=5, label='particle source')
    #plt.plot(r_line,z_line,'g',linewidth=2,label='PEC line integration')
    
    return intensity

if __name__ == "__main__":
    #plot_history2D()
    plot_surf_nc()
    #spec_line_integration()