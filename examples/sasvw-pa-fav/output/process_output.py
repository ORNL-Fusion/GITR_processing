import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import numpy as np
from scipy import special
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import netCDF4
import solps

def init(W_indices = np.arange(11,22)):
    profilesFile = '../input/plasmaProfiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    
    #check which target is the region of interest
    r_inner_target = profiles.variables['r_inner_target'][:]
    z_inner_target = profiles.variables['z_inner_target'][:]
    rmrs = profiles.variables['rmrs_inner_target_midpoints'][W_indices]
    
    #plot r,z to check that target indices are restricted to W surfaces
    plt.close()
    plt.plot(r_inner_target, z_inner_target)
    
    #set plotting style defaults
    plt.rcParams.update({'font.size':11.5})
    plt.rcParams.update({'lines.linewidth':20})
    plt.rcParams.update({'lines.markersize':1})

    return profiles, W_indices, r_inner_target, z_inner_target, rmrs

def plot_gitr_gridspace():
    profiles, W_indices, R, Z, rmrs = init()
    
    #import vars from plasmaProfiles.nc
    gridr = profiles.variables['gridr'][:]
    gridz = profiles.variables['gridz'][:]
    values = profiles.variables['values'][:]
    te = profiles.variables['te'][:]
    
    #set off grid space to 0
    off_grid_inds = np.where(te <= 0.0)
    te[off_grid_inds] = np.nan;
    
    #reformat
    nx,ny,crx,cry,region = solps.read_b2f_geometry('../setup/assets/b2fgmtry')
    meshgridr,meshgridz = np.meshgrid(gridr,gridz)
    print(meshgridr.shape)
    
    #plot gridspace
    plt.close()
    plt.plot(R,Z, '-k', linewidth=2)
    plt.pcolormesh(gridr,gridz,te,shading='nearest')
    plt.colorbar()
    plt.xlim(1.4,1.55)
    plt.ylim(1.05,1.25)
    return

def plot_particle_source():
    profiles, W_indices, R, Z, rmrs = init()
    particleSource = netCDF4.Dataset("../input/particleSource.nc", "r", format="NETCDF4")
    
    x = particleSource.variables['x'][:]
    y = particleSource.variables['y'][:]
    z = particleSource.variables['z'][:]
    
    plt.close()
    plt.hist(z,bins=10)
    plt.xlabel('z [m]')
    plt.title('Spatial Distribution in Z \n nP='+str(len(x)))
    plt.savefig('plots/zhist.png')
    
    plt.close()
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.scatter(x,z,marker='_')
    plt.axis('scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Spatial Particle Source \n nP='+str(len(x)))
    plt.savefig('plots/particleSource.png')

def plot_history2D(history_file='history.nc', \
                   basic=1, continuousChargeState=0, endChargeState=0, \
                   plot_particle_source=0, markersize=0):
    
    if plot_particle_source:
        particleSource = netCDF4.Dataset("../input/particleSource.nc", "r", format="NETCDF4")
        x0 = particleSource.variables['x'][:]
        z0 = particleSource.variables['z'][:]
    
    profiles, W_indices, R, Z, rmrs = init()
    history = netCDF4.Dataset(history_file, "r", format="NETCDF4")

    plt.rcParams.update({'lines.linewidth':1})
    plt.rcParams.update({'lines.markersize':markersize})
    plt.rcParams.update({'font.size':16})

    nP = len(history.dimensions['nP'])
    print('nP:',nP)
    nT = len(history.dimensions['nT'])
    x = history.variables['x']
    z = history.variables['z']
    charge = history.variables['charge']

    #plt.close()
    plt.close()
    if plot_particle_source: plt.scatter(x0,z0,marker='o',s=10)
    plt.plot(R,Z,'-k',linewidth=0.7)
    plt.axis('scaled')
    plt.xlabel('x [m]')
    plt.ylabel('z [m]')
    plt.title('W Impurity Trajectories', fontsize=20)
        
    #define charge state to color mapping
    colors = {0:'black', 1:'firebrick', 2:'darkorange', 3:'gold', 4:'limegreen', 5:'dodgerblue', \
              6:'mediumpurple', 7:'darkviolet', 8:'darkmagenta', 9:'deep pink', 10:'gray'}
    
    # all particle source vars ordered as (nP, nT)
    if basic==1:
        for p in range(0,nP,100):
            plt.plot(x[p][:],z[p][:])
            #plt.scatter(x[p][:],z[p][:],marker='_',s=50,c='k')
    
    if continuousChargeState==1:
        for p in np.arange(0,nP,1):
            print(p,'out of', nP)
            for t in np.arange(0,nT,100):
                plt.plot(x[p][t:t+100],z[p][t:t+100], colors[charge[p][t]])

    if endChargeState==1:
        for p in range(0,nP):
            plt.plot(x[p][:],z[p][:], colors[charge[p][-1]])
        plt.title('Particle Tracks by End Charge State')
    
    legend_dict = {'+0':'black', '+1':'firebrick', '+2':'darkorange', '+3':'gold', '+4':'limegreen', '+5':'dodgerblue', \
              '+6':'mediumpurple', '+7':'darkviolet', '+8':'darkmagenta', '+9':'deeppink', '+10':'gray'}
    patchList = []
    for key in legend_dict:
        data_key = mpatches.Patch(color=legend_dict[key], label=key)
        patchList.append(data_key)

    if basic==0: plt.legend(handles=patchList, fontsize=12)
    
    plt.xlim(1.3, 1.6)
    plt.ylim(1, 1.3)
    plt.show(block=True)
    plt.savefig('plots/history.pdf')
    plt.close()

    return

def plot_surf_nc(pps_per_nP, nP10, nT10, \
                 Bangle_shift_indices, \
                 surface_file="surface.nc", \
                 gitr_rz='../setup/assets/gitr_rz.txt', \
                 W_fine_file='../setup/assets/W_fine.txt', \
                 rmrs_fine_file='../setup/assets/rmrs_fine.txt', \
                 plot_cumsum=0):
    
    profiles, W_indices, r_inner_target, z_inner_target, rmrs = init()
    surface = netCDF4.Dataset(surface_file, "r", format="NETCDF4")
    #print(surface.variables['grossErosion'][:])
    
    #calculate area from wall
    #import wall geometry to plot over
    with open(gitr_rz, 'r') as file:
        wall = file.readlines()
    
    #import W surface indices
    with open(W_fine_file, 'r') as file:
        W_fine = file.readlines()
    W_fine = np.array(W_fine,dtype='int')
    
    #import refined rmrs at the W surface
    with open(rmrs_fine_file, 'r') as file:
        rmrs_fine = file.readlines()   
    rmrsFine = np.array(rmrs_fine,dtype='float')
    
    R = np.zeros(len(wall))
    Z = np.zeros(len(wall))
    for i,line in enumerate(wall):
        point = line.split()
        R[i] = float(point[0])
        Z[i] = float(point[1])
    
    R = R[W_fine]
    Z = Z[W_fine]
    
    r1 = R[:-1]
    r2 = R[1:]
    z1 = Z[:-1]
    z2 = Z[1:]
    
    dist = np.sqrt(np.power(r1-r2,2) + np.power(z1-z2,2))
    area = np.pi*(r1+r2)*dist

    grossEro = (surface.variables['grossErosion'][:])
    print(grossEro)
    grossDep = (surface.variables['grossDeposition'][:])
    print(grossDep)
    netDep = (grossDep-grossEro)
    
    print('total gross eroded comp particles',sum(grossEro))
    print('total redeposited comp particles',sum(grossDep))
    print('total net deposited comp particles',sum(netDep))
    print('ghost cells',grossEro[-1],grossDep[-1],netDep[-1])
    
    #grossEro = np.average([grossEro[:-1], grossEro[1:]],axis=0)*pps_per_nP/area
    #grossDep = np.average([grossDep[:-1], grossDep[1:]],axis=0)*pps_per_nP/area
    grossEro = grossEro[:-1]*pps_per_nP/area
    grossDep = grossDep[:-1]*pps_per_nP/area
    netDep = netDep[:-1]*pps_per_nP/area
    
    print('rmrs length',len(rmrsFine))
    print('surf length',len(grossEro))
    print('total gross eroded flux',sum(grossEro))
    print('total redeposited flux',sum(grossDep))
    print('total net deposited flux',sum(netDep))
    
    fluxC = 2.22351742230795e22
    grossEro_norm = grossEro/fluxC
    grossDep_norm = grossDep/fluxC
    netDep_norm = netDep/fluxC
    
    grossEro_cumsum = np.cumsum(grossEro)
    grossDep_cumsum = np.cumsum(grossDep)
    netDep_cumsum = np.cumsum(netDep)
    
    #take gross erosion of a slice in the filterscope range
    r_sp, z_sp = 1.49829829, 1.19672716
    
    r1_start, z1_start = 1.491288, 1.21629
    r1_end, z1_end = 1.49274, 1.22149
    rmrs1_start = -1*np.sqrt((z1_start-z_sp)**2 + (r1_start-r_sp)**2)
    rmrs1_end = -1*np.sqrt((z1_end-z_sp)**2 + (r1_end-r_sp)**2)
    
    r2_start, z2_start = 1.492041, 1.19418 
    r2_end, z2_end = 1.492716, 1.18709 
    rmrs2_start = 1*np.sqrt((z2_start-z_sp)**2 + (r2_start-r_sp)**2)
    rmrs2_end = 1*np.sqrt((z2_end-z_sp)**2 + (r2_end-r_sp)**2)
    
    r3_start, z3_start = 1.493556, 1.16204
    r3_end, z3_end = 1.490536, 1.1552 
    rmrs3_start = np.sqrt((z3_start-z_sp)**2 + (r3_start-r_sp)**2)
    rmrs3_end = np.sqrt((z3_end-z_sp)**2 + (r3_end-r_sp)**2)
    
    V1_grossEro = 0
    V1_area = 0
    V2_grossEro = 0
    V2_area = 0
    V3_grossEro = 0
    V3_area = 0
    for i,v in enumerate(rmrsFine):
        if v >= rmrs1_end and v <= rmrs1_start:
            V1_grossEro += grossEro[i] * (rmrsFine[i+1] - rmrsFine[i])
            V1_area += rmrsFine[i+1] - rmrsFine[i]
        elif v >= rmrs2_start and v <= rmrs2_end:
            V2_grossEro += grossEro[i] * (rmrsFine[i+1] - rmrsFine[i])
            V2_area += rmrsFine[i+1] - rmrsFine[i]
        elif v >= rmrs3_start and v <= rmrs3_end:
            V3_grossEro += grossEro[i] * (rmrsFine[i+1] - rmrsFine[i])
            V3_area += rmrsFine[i+1] - rmrsFine[i]
    V1_grossEro = V1_grossEro / V1_area
    V2_grossEro = V2_grossEro / V2_area
    V3_grossEro = V3_grossEro / V3_area
    
    print('gross erosion in View 1:', V1_grossEro)
    print('gross erosion in View 2:', V2_grossEro)
    print('gross erosion in View 3:', V3_grossEro)
    
    plt.rcParams.update({'font.size':16})
    plt.rcParams.update({'lines.linewidth':3}) 
    
    plt.close()
    plt.axvspan(rmrs1_start, rmrs1_end, color='#f7bc00', alpha=0.5)
    plt.axvspan(rmrs2_start, rmrs2_end, color='lightsalmon', alpha=0.5)
    plt.axvspan(rmrs3_start, rmrs3_end, color='#f99301', alpha=0.5)
    
    
    plt.plot(rmrsFine,np.zeros(len(rmrsFine)),'gray')
    if Bangle_shift_indices != []:
        for i in Bangle_shift_indices:
            plt.axvline(x=rmrs[i], color='k', linestyle='dotted')#, label='\u0394\u03A8$_B$')
    
    plt.plot(rmrsFine,grossEro_norm,'r', label='Gross Erosion', linewidth=10)
    plt.plot(rmrsFine,grossDep_norm,'g', label='Redeposition')
    plt.plot(rmrsFine,netDep_norm,'k', label='Net Deposition')
    
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('\u0393$_{W,outgoing}$ / \u0393$_{C,incoming}$')
    plt.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
    plt.legend()#loc='upper left')
    plt.title('GITR Predicted Erosion and \n Redeposition Profiles, nP=1e'+str(nP10)+', nT=1e'+str(nT10))
    plt.savefig('plots/surface.png')
    
    if plot_cumsum:
        plt.close()
        plt.plot(rmrsFine,grossEro_cumsum,'r', label='Gross Erosion')
        plt.plot(rmrsFine,grossDep_cumsum,'g', label='Redeposition')
        plt.plot(rmrsFine,netDep_cumsum,'k', label='Net Erosion')
        plt.yscale('log')
        plt.xlabel('D-Dsep [m]')
        plt.ylabel('Flux [m$^{-2}$s$^{-1}$]')
        plt.legend(loc='upper left')
        plt.title('GITR Predicted Cumulative Sum of \n Erosion and Redeposition Profiles',fontsize=20)
        plt.savefig('plots/surface_cumsum.pdf')
    
def spectroscopy(pps_per_nP, View=3, \
                 specFile='spec.nc'):
    
    spec = netCDF4.Dataset(specFile, "r", format="NETCDF4")
    profiles, W_indices, R, Z, rmrs = init()
    
    gridr = spec.variables['gridR'][:]
    gridz = spec.variables['gridZ'][:]
    dr = gridr[1]-gridr[0]
    dz = gridz[1]-gridz[0]
    gridr = np.append(gridr, gridr[-1]+dr)
    gridz = np.append(gridz, gridz[-1]+dz)
    rr, zz = np.meshgrid(gridr,gridz)
    density_unitless = spec.variables['n'][:]
    density_neutrals = density_unitless[0]
    #density_neutrals = np.sum(density_unitless, axis=0)
    
    plt.pcolor(gridr,gridz,density_neutrals)
    plt.axis('Scaled')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.colorbar(label='# Computational Particles')
    plt.title('Raw W0 Spec Output form GITR')
    plt.savefig('plots/spec_compParticles.png')
    
    filterscope_diameter = 0.00702068
    tokamak_circumference = 2 * np.pi * 1.4925
    toroidal_fraction = filterscope_diameter / tokamak_circumference
    
    r_mid = dr/2 + rr[:-1,:]
    dV = 2 * np.pi * r_mid[:,:-1] * dr * dz * toroidal_fraction
    # dA = np.pi * (0.00702068/2)**2
    density_volumetric = pps_per_nP * (density_neutrals/dV) * 1e-8 * 3e3 # W0 m-2 s-1
    #neutral_flux = density_volumetric * 3e3 # W0 m-2 s-1
    #print('TESTTT',neutral_flux)
    
    plt.close()
    plt.plot(R,Z)
    plt.pcolor(gridr,gridz,density_volumetric)
    plt.axis('Scaled')
    plt.xlim(min(R),max(R))
    plt.ylim(min(Z),max(Z))
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.colorbar(label='\n Density [m$^{-3}$ s$^{-1}$]')
    plt.title('Toroidally Integrated W0 Density')
    plt.savefig('plots/spec_density.png')

    view_fraction = toroidal_fraction * np.pi/4
    density_sliced = view_fraction * density_volumetric
    
    if View==3:
        line1 = -0.97*gridr + 2.600823
        line2 = -0.906*gridr + 2.515464
        
        if True:
            plt.close()
            plt.plot(R,Z)
            plt.plot(gridr,line1,'orange')
            plt.plot(gridr,line2,'orange')
            plt.pcolor(gridr,gridz,density_sliced)
            plt.axis('Scaled')
            plt.xlim(gridr[0],gridr[-1])
            plt.ylim(gridz[0],gridz[-1])
            plt.xlabel('r [m]')
            plt.ylabel('z [m]')
            plt.colorbar(label='\n Density [m$^{-3}$ s$^{-1}$]')
            plt.title('Toroidal Slice of W0 Density')
            plt.savefig('plots/spec_density_sliced.png')
        
        rstart, rend = 39, 72
        zstart, zend = 40, 80
        gridr = gridr[rstart:rend]
        gridz = gridz[zstart:zend]
        fscope = density_sliced[zstart:zend, rstart:rend]
        
        fscope[np.where(fscope==0)] = 'NaN'
        #fscope[-5:,23:] = 'NaN'
    
        line1 = -0.97*gridr + 2.600823
        line2 = -0.906*gridr + 2.515464
        
        if True:
            plt.close()
            plt.plot(R,Z)
            line1 = -0.97*gridr + 2.600823
            line2 = -0.906*gridr + 2.515464    
            plt.plot(gridr,line1,'orange')
            plt.plot(gridr,line2,'orange')
            plt.pcolor(gridr,gridz,fscope,shading='auto')
            plt.axis('Scaled')
            plt.xlim(gridr[0],gridr[-1])
            plt.ylim(gridz[0],gridz[-1])
            plt.xlabel('r [m]')
            plt.ylabel('z [m]')
            plt.colorbar(label='\n Density [m$^{-3}$ s$^{-1}$]')
    
    fscope[np.isnan(fscope)] = 0
    print('W0 Neutral Flux:', np.sum(fscope)) # in units of m-3 s-1
    plt.title('Toroidal Slice of W0 Density \n W0: %10e' %np.sum(fscope))
    plt.savefig('plots/spec_filterscope.png')
    
def analyze_leakage():
    return    

def analyze_forces(varString, component, rzlim=True, colorbarLimits=[]):
    #import wall geometry to plot over
    gitr_rz='../setup/assets/gitr_rz.txt'
    with open(gitr_rz, 'r') as file:
        wall = file.readlines()
        
    r_wall = np.zeros(len(wall))
    z_wall = np.zeros(len(wall))
    for i,line in enumerate(wall):
        point = line.split()
        r_wall[i] = float(point[0])
        z_wall[i] = float(point[1])
    
    profilesFile = '../input/plasmaProfiles.nc'
    profiles = netCDF4.Dataset(profilesFile)
    gridr = profiles.variables['gridr'][:]
    gridz = profiles.variables['gridz'][:]
    
    Br = profiles.variables['Br'][:]
    Bt = profiles.variables['Bt'][:]
    Bz = profiles.variables['Bz'][:]
    Er = profiles.variables['Er'][:]
    Et = profiles.variables['Et'][:]
    Ez = profiles.variables['Ez'][:]
    charge = profiles.variables['charge'][:]
    density = profiles.variables['ni'][:]
    ti = profiles.variables['ti'][:]
    te = profiles.variables['te'][:]
    vr_background = profiles.variables['vr'][:]
    vt_background = profiles.variables['vt'][:]
    vz_background = profiles.variables['vz'][:]
    
    q = (1.602e-19)*np.ones(np.shape(Br))
    vr = -1000*np.ones(np.shape(Br))
    vt = -500*np.ones(np.shape(Br))
    vz = 700*np.ones(np.shape(Br))
    
    if varString=='q(v x B)':
        vartype = 'F'
        titleString = ' = (q(v x B))'
        
        Fr = q*(vt*Bz - vz*Bt)
        Ft = q*(vz*Br - vr*Bz)
        Fz = q*(vr*Bt - vt*Br)
    
    elif varString=='qE':
        vartype = 'F'
        titleString = ' = (qE))'
        
        Fr = q*Er
        Ft = q*Et
        Fz = q*Ez
    
    elif varString=='Lorentz':
        vartype = 'F'
        titleString = ' = (q(E + v x B))'
        
        Fr = q*(Er + vt*Bz - vz*Bt)
        Ft = q*(Et + vz*Br - vr*Bz)
        Fz = q*(Ez + vr*Bt - vt*Br)
    
    elif varString=='ExB drift':
        vartype = 'v'
        titleString = ' = ($v_{ExB}$)'
        
        dr = (np.max(gridr[1:]-gridr[:-1])-np.min(gridr[1:]-gridr[:-1]))/2
        dz = (np.max(gridz[1:]-gridz[:-1])-np.min(gridz[1:]-gridz[:-1]))/2
        Wmass_amu = 183.84
        Wmass_kg = Wmass_amu * 1.66054e-27
        B_mag2 = Br**2 + Bt**2 + Bz**2
        B_mag = np.sqrt(B_mag2)
        v_mag = np.sqrt(vr**2 + vt**2 + vz**2)
        v_parallel = np.sqrt((v_mag*Br/B_mag)**2 + (v_mag*Bt/B_mag)**2 + (v_mag*Bz/B_mag)**2)
        v_perp = np.abs(v_mag) - np.abs(v_parallel)
        r_Larmor = Wmass_kg * v_perp / (q * B_mag)
        print('r_Larmor:', np.average(r_Larmor))
        
        ExB_r = Et*Bz - Ez*Bt
        ExB_t = Ez*Br - Er*Bz
        ExB_z = Er*Bt - Et*Br
        
        del_ExB_z = np.zeros(np.shape(Br))
        del_ExB_z[0,:] = (ExB_z[1,:] - ExB_z[0,:])/dz
        del_ExB_z[-1,:] = (ExB_z[-1,:] - ExB_z[-2,:])/dz
        for i in range(1,len(gridz)-1):
            del_ExB_z[i,:] = (ExB_z[i+1,:] - ExB_z[i-1,:])/(2*dz)
        
        del2_ExB_z = np.zeros(np.shape(Br))
        del2_ExB_z[0,:] = (del_ExB_z[1,:] - del_ExB_z[0,:])/dz
        del2_ExB_z[-1,:] = (del_ExB_z[-1,:] - del_ExB_z[-2,:])/dz
        for i in range(1,len(gridz)-1):
            del2_ExB_z[i,:] = (del_ExB_z[i+1,:] - del_ExB_z[i-1,:])/(2*dz)
        
        del_ExB_r = np.zeros(np.shape(Br))
        del_ExB_r[:,0] = gridr[0] * (ExB_r[:,1] - ExB_r[:,0])/dr
        del_ExB_r[:,-1] = gridr[-1] * (ExB_r[:,-1] - ExB_r[:,-2])/dr
        for i in range(1,len(gridr)-1):
            del_ExB_r[:,i] = gridr[i] * (ExB_r[:,i+1] - ExB_r[:,i-1])/(2*dr)
        
        del2_ExB_r = np.zeros(np.shape(Br))
        del2_ExB_r[:,0] = (del_ExB_r[:,1] - del_ExB_r[:,0])/(gridr[0]*dr)
        del2_ExB_r[:,-1] = (del_ExB_r[:,-1] - del_ExB_r[:,-2])/(gridr[-1]*dr)
        for i in range(1,len(gridr)-1):
            del2_ExB_r[:,i] = (del_ExB_r[:,i+1] - del_ExB_r[:,i-1])/(2*gridr[i]*dr)
        
        Fr = ExB_r/B_mag2 + (r_Larmor**2 * del2_ExB_r)/(4 * B_mag2)
        Ft = ExB_t/B_mag2
        Fz = ExB_z/B_mag2 + (r_Larmor**2 * del2_ExB_z)/(4 * B_mag2)
    
    elif varString=='drag':
        vartype = 'F'
        titleString = ' = ($m \\nu_S U_{\parallel}$)'
        
        amu = 183.84
        amu_D = 2
        amu_C = 12.011
        Z_D = 1
        Z_C = 6
        EPS0 = 8.854187e-12 #epsilon_0 = vacuum permissivity
        Q = 1.60217662e-19
        MI = 1.6737236e-27
        
        #inverese acceleration of D and C background ions
        a_D = amu_D*MI/(2*ti*Q)
        a_C = amu_C*MI/(2*ti*Q)
        
        relativeVelocity_r = vr - vr_background
        relativeVelocity_t = vt - vt_background
        relativeVelocity_z = vz - vz_background
        velocityNorm = np.sqrt(relativeVelocity_r**2 + relativeVelocity_t**2 + relativeVelocity_z**2)
        
        #wtf is this
        lam_d_D = np.sqrt(EPS0*te/(density*(Z_D**2)*Q))
        lam_d_C = np.sqrt(EPS0*te/(density*(Z_C**2)*Q))
        lam_D = 12.0*np.pi*density*lam_d_D**3/charge
        lam_C = 12.0*np.pi*density*lam_d_C**3/charge
        gam_electron_background_D = 0.238762895*charge**2*np.log(lam_D)/(amu*amu)
        gam_electron_background_C = 0.238762895*charge**2*np.log(lam_C)/(amu*amu)
        nu_0_D = gam_electron_background_D*density/velocityNorm**3
        nu_0_C = gam_electron_background_C*density/velocityNorm**3
        
        #wtf is any of this
        xx_D = velocityNorm**2 * a_D
        xx_C = velocityNorm**2 * a_C
        psi_prime_D = 2.0*np.sqrt(xx_D/np.pi)*np.exp(-xx_D)
        psi_prime_C = 2.0*np.sqrt(xx_C/np.pi)*np.exp(-xx_C)
        psi_psiprime_D = special.erf(np.sqrt(xx_D))
        psi_psiprime_C = special.erf(np.sqrt(xx_C))
        psi_D = psi_psiprime_D - psi_prime_D
        psi_C = psi_psiprime_C - psi_prime_C
        
        #Find zero velocity moment friction frequency:
        #nu_0 for W in D (nu_friction_D) and for W in C (nu_friction_C)
        nu_friction_D = (1+amu/amu_D)*psi_D*nu_0_D
        nu_friction_C = (1+amu/amu_C)*psi_C*nu_0_C
        
        B_mag2 = Br**2 + Bt**2 + Bz**2
        B_mag = np.sqrt(B_mag2)
        v_mag = np.sqrt(vr**2 + vt**2 + vz**2)
        v_parallel_r = v_mag*Br/B_mag
        v_parallel_t = v_mag*Bt/B_mag
        v_parallel_z = v_mag*Bz/B_mag
        
        #Calulate drag force component of the Fokker-Plank collisional forces
        drag_D_r = amu * MI * nu_friction_D * v_parallel_r
        drag_D_t = amu * MI * nu_friction_D * v_parallel_t
        drag_D_z = amu * MI * nu_friction_D * v_parallel_z
        drag_C_r = amu * MI * nu_friction_C * v_parallel_r
        drag_C_t = amu * MI * nu_friction_C * v_parallel_t
        drag_C_z = amu * MI * nu_friction_C * v_parallel_z
        
        Fr = -1* (drag_D_r + drag_C_r)
        Ft = -1* (drag_D_t + drag_C_t)
        Fz = -1* (drag_D_z + drag_C_z)
        
    elif varString=='friction':
        vartype = 'F'
        titleString = ' = ($m \\nu_S U_{\parallel}$)'
        
        amu = 183.84
        amu_D = 2
        amu_C = 12.011
        Z_D = 1
        Z_C = 6
        EPS0 = 8.854187e-12 #epsilon_0 = vacuum permissivity
        Q = 1.60217662e-19
        MI = 1.6737236e-27
        
        #inverese acceleration of D and C background ions
        a_D = amu_D*MI/(2*ti*Q)
        a_C = amu_C*MI/(2*ti*Q)
        
        relativeVelocity_r = vr - vr_background
        relativeVelocity_t = vt - vt_background
        relativeVelocity_z = vz - vz_background
        velocityNorm = np.sqrt(relativeVelocity_r**2 + relativeVelocity_t**2 + relativeVelocity_z**2)
        
        #wtf is this
        lam_d_D = np.sqrt(EPS0*te/(density*(Z_D**2)*Q))
        lam_d_C = np.sqrt(EPS0*te/(density*(Z_C**2)*Q))
        lam_D = 12.0*np.pi*density*lam_d_D**3/charge
        lam_C = 12.0*np.pi*density*lam_d_C**3/charge
        gam_electron_background_D = 0.238762895*charge**2*np.log(lam_D)/(amu*amu)
        gam_electron_background_C = 0.238762895*charge**2*np.log(lam_C)/(amu*amu)
        nu_0_D = gam_electron_background_D*density/velocityNorm**3
        nu_0_C = gam_electron_background_C*density/velocityNorm**3
        
        #wtf is any of this
        xx_D = velocityNorm**2 * a_D
        xx_C = velocityNorm**2 * a_C
        psi_prime_D = 2.0*np.sqrt(xx_D/np.pi)*np.exp(-xx_D)
        psi_prime_C = 2.0*np.sqrt(xx_C/np.pi)*np.exp(-xx_C)
        psi_psiprime_D = special.erf(np.sqrt(xx_D))
        psi_psiprime_C = special.erf(np.sqrt(xx_C))
        psi_D = psi_psiprime_D - psi_prime_D
        psi_C = psi_psiprime_C - psi_prime_C
        
        #Find zero velocity moment friction frequency:
        #nu_0 for W in D (nu_friction_D) and for W in C (nu_friction_C)
        nu_friction_D = (1+amu/amu_D)*psi_D*nu_0_D
        nu_friction_C = (1+amu/amu_C)*psi_C*nu_0_C
        
        B_mag2 = Br**2 + Bt**2 + Bz**2
        B_mag = np.sqrt(B_mag2)
        v_mag = np.sqrt(vr**2 + vt**2 + vz**2)
        v_parallel_r = v_mag*Br/B_mag
        v_parallel_t = v_mag*Bt/B_mag
        v_parallel_z = v_mag*Bz/B_mag
        
        #Calulate drag force component of the Fokker-Plank collisional forces
        drag_D_r = amu * MI * nu_friction_D * v_parallel_r
        drag_D_t = amu * MI * nu_friction_D * v_parallel_t
        drag_D_z = amu * MI * nu_friction_D * v_parallel_z
        drag_C_r = amu * MI * nu_friction_C * v_parallel_r
        drag_C_t = amu * MI * nu_friction_C * v_parallel_t
        drag_C_z = amu * MI * nu_friction_C * v_parallel_z
        
        Fr_drag = -1* (drag_D_r + drag_C_r)
        Ft_drag = -1* (drag_D_t + drag_C_t)
        Fz_drag = -1* (drag_D_z + drag_C_z)
        
        nu_parallel_D = psi_D/xx_D*nu_0_D
        nu_parallel_C = psi_C/xx_C*nu_0_C
        nu_deflection_D = 2*(psi_psiprime_D - psi_D/(2*xx_D))*nu_0_D
        nu_deflection_C = 2*(psi_psiprime_C - psi_C/(2*xx_C))*nu_0_C
        nu_energy_D = 2*(amu/amu_D*psi_D - psi_prime_D)*nu_0_D
        nu_energy_C = 2*(amu/amu_C*psi_C - psi_prime_C)*nu_0_C
        
        
        
        #vx_relative = velocityRelativeNorm*(1.0-0.5*nuEdt)*((1.0 + coeff_par) * parallel_direction[0] + std::abs(n2)*(coeff_perp1 * perp_direction1[0] + coeff_perp2 * perp_direction2[0])) - velocityRelativeNorm*dt*nu_friction*parallel_direction[0];
    
    elif varString=='gradTi':
        vartype = 'F'
        titleString = ' = ($placeholder$)'
        
        gradTir = profiles.variables['gradTir'][:]
        gradTit = profiles.variables['gradTit'][:]
        gradTiz = profiles.variables['gradTiz'][:]
        gradTer = profiles.variables['gradTer'][:]
        gradTet = profiles.variables['gradTet'][:]
        gradTez = profiles.variables['gradTez'][:]
        
        mu_D = amu / (amu_D + amu)
        mu_C = amu / (amu_C + amu)
        alpha = 0.71*charge**2
        beta_D = 3 * (mu_D + 5*np.sqrt(2.0) * charge**2 * (1.1*mu_D**(5 / 2) - 0.35*mu_D**(3 / 2)) - 1) / (2.6 - 2*mu_D + 5.4*mu_D**2)
        beta_C = 3 * (mu_C + 5*np.sqrt(2.0) * charge**2 * (1.1*mu_C**(5 / 2) - 0.35*mu_C**(3 / 2)) - 1) / (2.6 - 2*mu_C + 5.4*mu_C**2)
        beta = beta_D + beta_C #assume that D and C are in thermal equilibrium
        
        Fr = alpha*gradTer + beta*gradTir
        Ft = alpha*gradTet + beta*gradTit
        Fz = alpha*gradTez + beta*gradTiz
    
    Fr[np.where(Fr==0)] = 'nan'
    Ft[np.where(Ft==0)] = 'nan'
    Fz[np.where(Fz==0)] = 'nan'
    
    gridrz = [gridr, gridz, r_wall, z_wall]
    if component == 'r': plot_forces(Fr, vartype+'r'+titleString+'$_r$', gridrz, vartype, rzlim, colorbarLimits)
    if component == 't': plot_forces(Ft, vartype+'t'+titleString+'$_t$', gridrz, vartype, rzlim, colorbarLimits)
    if component == 'z': plot_forces(Fz, vartype+'z'+titleString+'$_z$', gridrz, vartype, rzlim, colorbarLimits)
    
    return 

def plot_forces(var, titleString, gridrz, vartype='F', rzlim=True, colorbarLimits=[]):
    [gridr, gridz, r_wall, z_wall] = gridrz
    plt.rcParams.update({'pcolor.shading':'auto'})
    plt.rcParams.update({'image.cmap':'rainbow'})
    rlim = [1.375, 1.575]
    zlim = [1.05, 1.25]
    
    if rzlim:
        r_indices = np.empty(0, dtype='int')
        z_indices = np.empty(0, dtype='int')
        
        for i,v in enumerate(gridr):
            if v>=rlim[0] and v<=rlim[1]:
                r_indices = np.append(r_indices, i)
        
        for i,v in enumerate(gridz):
            if v>=zlim[0] and v<=zlim[1]:
                z_indices = np.append(z_indices, i)

        gridr = gridr[r_indices]
        gridz = gridz[z_indices]
        var = var[z_indices[0]:z_indices[-1]+1,r_indices[0]:r_indices[-1]+1]
    
    plt.plot(r_wall,z_wall,'-k')
    plot = plt.pcolor(gridr,gridz,var)
    
    if colorbarLimits != []: plot.set_clim(vmin=colorbarLimits[0],vmax=colorbarLimits[1])
    if vartype=='F': plt.colorbar(label='\n Force [N]')
    if vartype=='v': plt.colorbar(label='\n Velocity [m/s]')
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title(titleString)
    plt.axis('Scaled')
    plt.xlim(rlim)
    plt.ylim(zlim)
    
    return

if __name__ == "__main__":
    analyze_forces('friction', 'r', rzlim=True, colorbarLimits=[])
    
    #init()
    #plot_gitr_gridspace()
    #plot_particle_source()
    #plot_history2D('history.nc')
    #plot_history2D('history-alpine.nc', plot_particle_source=1, markersize=2)
    #plot_history2D("../../../../GITR/scratch/output/history.nc")
    #plot_history2D("history_nP5e2_nT1e5.nc")
    #plot_surf_nc(7.582961536113231e+17, 'surface-alpine.nc')
    #plot_surf_nc(1006929636574578.9, 5, 5, [3,8,9], "surface_p5t5.nc")
    #plot_surf_nc(1063289762078132.4, "surface.nc")
    #plot_surf_nc(37914807680566.16, "/Users/Alyssa/Dev/SAS-VW-Data/netcdf_data/nP5/surf-5-6.nc")
    #spectroscopy(3791480768056.615,specFile='specP6T6.nc')
