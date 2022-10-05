import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))
sys.path.insert(0, os.path.abspath('../../../../pyGITR/pyGITR'))

import numpy as np
import matplotlib.pyplot as plt

import netCDF4
import gitr
import solps
import Particles

def point_source(nP = int(2e2)):
    x = 1.44*np.ones(nP)
    y = np.zeros(nP)
    z = 1.15*np.ones(nP)
    vx = 500*np.ones(nP)
    vy = 5000*np.zeros(nP)
    vz = 100*np.zeros(nP)

    #########################################
    #make NetCDF Particle Source file
    #########################################

    rootgrp = netCDF4.Dataset("particleSource.nc", "w", format="NETCDF4")
    npp = rootgrp.createDimension("nP", nP)
    xxx = rootgrp.createVariable("x","f8",("nP"))
    yyy = rootgrp.createVariable("y","f8",("nP"))
    zzz = rootgrp.createVariable("z","f8",("nP"))
    vxx = rootgrp.createVariable("vx","f8",("nP"))
    vyy = rootgrp.createVariable("vy","f8",("nP"))
    vzz = rootgrp.createVariable("vz","f8",("nP"))
    xxx[:] = x
    yyy[:] = y
    zzz[:] = z
    vxx[:] = vx
    vyy[:] = vy
    vzz[:] = vz
    rootgrp.close()


def simple2D(nP, \
            geom = '../input/gitrGeometry.cfg', \
            targFile = 'assets/rightTargOutput', \
            coordsFile = 'assets/right_target_coordinates.txt', \
            profilesFile = '../input/profiles.nc', \
            ftDFile = 'assets/ftridynBackgroundD.nc', \
            configuration = 'random', \
            plot_variables = 0, \
            r_W = None, z_W = None):
    
    #import r1,r2,z1,z2 coordinates
    x1, x2, z1, z2, length, Z, slope, inDir = gitr.plot2dGeom(geom)
    coords = np.loadtxt(coordsFile, dtype='float',skiprows=1,delimiter=' ')
    gitr_inds = [int(i) for i in coords[0:-1,3]]
    r1 = x1[gitr_inds] #coords[:,4]
    r2 = x2[gitr_inds] #coords[:,5]
    z1 = z1[gitr_inds] #coords[:,6]
    z2 = z2[gitr_inds] #coords[:,7]
    slope = slope[gitr_inds] #coords[:,7]
    area = 0.12*np.sqrt(np.power(r1-r2,2) + np.power(z1-z2,2))#*np.pi*(r1+r2)

    #get indices of rcoord and r_targ that are W
    W_ind = np.empty(len(r_W),dtype=int)
    for i in range(len(r_W)):
        W_ind[i] = np.where(r1 == r_W[i])[0]

    #define angle between material wall and major radius, x
    Alpha = np.abs(np.arctan((z2-z1) / (r2-r1)))
    Beta = np.abs(np.pi/2 - Alpha)
    Alpha = np.abs(Alpha[W_ind[:-1]]) #np.abs(np.rad2deg(Alpha[W_ind[:-1]]))
    Beta = np.abs(Beta[W_ind[:-1]]) #np.abs(np.rad2deg(Beta[W_ind[:-1]]))

    #specify W r,z coordinates and slope for the outer target
    r1 = r1[W_ind]
    z1 = z1[W_ind]
    slope = slope[W_ind[:-1]]


    #########################################
    #get W/s sputtered by D, He flux to wall
    #########################################

    #get flux grid from background D, C
    r_mid, z_mid, ti, ni, flux, te, ne = solps.read_target_file(targFile)

    #split total flux into D and C ion fluxes at the W surface
    ion_flux = np.abs(flux[:,1:][W_ind[:-1]]) #only bother with W surfaces
    Dflux = np.transpose(np.abs(flux[:,1][W_ind[:-1]]))
    Cflux = np.transpose(np.abs(flux[:,3:][W_ind[:-1]]))

    #get incoming ion energy and angle estimations where the integer input is z
    DsurfE, DsurfA = get_surf_profiles(profilesFile, r1, z1, slope, 1, plot_variables)
    CsurfE, CsurfA = get_surf_profiles(profilesFile, r1, z1, slope, 6, plot_variables)
    print('D energies',DsurfE) #should not sputter anything :'(
    print('C energies',CsurfE) #should sputter things

    #get sputtering yields for D0 and D1+ on W from fractal tridyn tables
    Dspyld = get_ft_spyld(DsurfE, DsurfA, ftDFile)[1] #file input includes He, which we aren't using
    Cspyld = get_analytic_spyld(CsurfE, CsurfA)
    print('D spyld', Dspyld)
    print('C spyld', Cspyld)

    #multiply incoming ion flux by Y_s to get sputtered W flux by each species
    sputt_flux = Dspyld*Dflux + Cspyld*np.sum(Cflux,axis=0)
    print('SPUTT FLUX',sputt_flux)

    #multiply by area to get the outgoing particles per second
    pps = np.multiply(sputt_flux,area[W_ind[:-1]])
    pps_weights = nP*pps/np.sum(pps)


    #########################################
    #get x,y,z distributions for sputtered W
    #########################################

    #confirm nP stays constant
    for i in range(len(pps_weights)): pps_weights[i] = round(pps_weights[i])
    nP_diff = int(nP-np.sum(pps_weights))

    if nP_diff > 0:
        for i in range(abs(nP_diff)):
            rand_index = np.random.choice(len(pps_weights))
            pps_weights[rand_index] += 1
    elif nP_diff < 0:
        for i in range(abs(nP_diff)):
            rand_index = np.random.choice(len(pps_weights))
            pps_weights[rand_index] -= 1
    pps_weights = pps_weights.astype(int)
    print('total nP', nP)
    print('nP(r_mid):', pps_weights)
    nP_diff = nP-np.sum(pps_weights)
    print('nP_diff should be 0: ', nP_diff)
    
    #populate a,b with Alpha,Beta
    a = np.zeros(nP)
    b = np.zeros(nP)
    counter=0
    for i in range(len(pps_weights)):
        a[counter:counter+pps_weights[i]] = Alpha[i]
        b[counter:counter+pps_weights[i]] = Beta[i]
        counter += pps_weights[i]

    #define adjustment into the sheath because particles can't start exactly on the wall
    adj = 1e-7

    #populate x,y,z with r_mid,0,z_mid
    if configuration == 'random': 
        x,y,z = random(nP,pps_weights,adj,slope,Beta, r1,z1)
    elif configuration == 'uniform': 
        x,y,z = uniform(nP,pps_weights,adj,slope,Beta, r1,z1)
    elif configuration == 'midpoint': 
        x,y,z = midpoints(nP,pps_weights, adj,slope,Beta, r1,z1)
    else:
        print('(x,y,z) configuration not set')

    if plot_variables == 1:
        zmid_surf = (z1[:-1]+z1[1:])/2
        plt.close()
        plt.plot(zmid_surf, Dflux)
        plt.xlabel('z [m]')
        plt.ylabel('Flux [#/m2s]')
        plt.title('Incoming D Ion Flux')
        plt.savefig('plots/surf_Dflux.png')

        plt.close()
        plt.plot(zmid_surf, np.sum(Cflux,axis=0))
        plt.xlabel('z [m]')
        plt.ylabel('Flux [#/m2s]')
        plt.title('Incoming C Ion Flux')
        plt.savefig('plots/surf_Cflux.png')

        plt.close()
        plt.plot(zmid_surf, Dspyld)
        plt.xlabel('z [m]')
        plt.ylabel('Yield')
        plt.title('Average Sputtering Yield for D on W')
        plt.savefig('plots/Dspyld.png')

        plt.close()
        plt.plot(zmid_surf, Cspyld)
        plt.xlabel('z [m]')
        plt.ylabel('Yield')
        plt.title('Average Sputtering Yield for C on W')
        plt.savefig('plots/Cspyld.png')

        plt.close()
        plt.plot(zmid_surf, sputt_flux)
        plt.xlabel('z [m]')
        plt.ylabel('Sputtered Flux [#/m2s]')
        plt.title('Flux of W Sputtered off Wall')
        plt.savefig('plots/sputt_flux.png')

        plt.close()
        plt.plot(r_W,z_W,'-k')
        plt.scatter(x,z)
        plt.axis('Scaled')
        plt.savefig('plots/test')


    #########################################
    #get vx,vy,vz from IEADs
    #########################################

    #use PyGITR to set up vx,vy,vz,E,theta,psi distributions
    PartDist = Particles.ParticleDistribution(nP, ListAttr=['vx','vy','vz'])

    vx = np.zeros(1)
    vy = np.zeros(1)
    vz = np.zeros(1)
    
    for i in range(len(pps_weights)):
    
        weight = int(pps_weights[i])
        m = np.sign(slope[i])
        
        #get IEADs for sputtered W
        E = PartDist.Generate(weight, 'Thomson')
        PolAng = PartDist.Generate(weight, 'SinCos', x=np.linspace(0,np.pi/2,weight))
        AziAng = PartDist.Generate(weight, 'Uniform', x=np.linspace(0,2*np.pi,weight))
        
        #convert IEADs to vx,vy,vz unit vectors in particle frame of ref
        PartDist.SetAttr('vx', np.multiply(np.cos(PolAng), np.cos(AziAng)))
        PartDist.SetAttr('vy', np.multiply(np.cos(PolAng), np.sin(AziAng)))
        PartDist.SetAttr('vz', m*np.sin(PolAng))
        
        vx_prime = PartDist.Particles['vx']
        vy_prime = PartDist.Particles['vy']
        vz_prime = PartDist.Particles['vz']

        #rotate vx,vy,vz from particle frame to lab frame
        PartDist.RotateAngle('v', -m*Alpha[i],0, Degree=False)
        
        vx_lab = PartDist.Particles['vx']
        vy_lab = PartDist.Particles['vy']
        vz_lab = PartDist.Particles['vz']

        if plot_variables == 1:
            plt.close()
            plt.scatter(vx_lab,vz_lab)
            plt.axis('Scaled')
            plt.xlabel('vx')
            plt.ylabel('vz')
            plt.title('SinCos Polar Angle in the Lab Frame')

        #convert unit vectors to vx,vy,vz
        W_kg = 183.84 * 1.6605e-27 #mass of W in kg
        vtot = np.sqrt(E*1.6022e-19/W_kg) #convert eV to m/s
        
        vx = np.append(vx, vtot*vx_lab)
        vy = np.append(vy, vtot*vy_lab)
        vz = np.append(vz, vtot*vz_lab)
    
    vx = np.delete(vx,0)
    vy = np.delete(vy,0)
    vz = np.delete(vz,0)

    if plot_variables == 1:
        #plot Thomson E dist
        plt.close()
        plt.hist(E,bins=100)
        plt.xlabel('Energy Bins [eV]')
        plt.ylabel('Histogram')
        plt.title('Thomson Energy Distribution')
        plt.savefig('plots/thomson')
        
        #plot particle framed v_dist relations
        plt.close()
        plt.scatter(vx_prime,vy_prime,s=0.3)
        plt.axis('Scaled')
        plt.xlabel('vx')
        plt.ylabel('vy')
        plt.title('Uniform Azimuthal Angle in the Particle Frame')
        plt.savefig('plots/vxvy_prime')
        
        plt.close()
        plt.scatter(vx_prime,vz_prime,s=0.3)
        plt.axis('Scaled')
        plt.xlabel('vx')
        plt.ylabel('vz')
        plt.title('SinCos Polar Angle in the Particle Frame')
        plt.savefig('plots/vxvz_prime')
        
        plt.close()
        plt.scatter(vx_lab,vz_lab,s=0.3)
        plt.axis('Scaled')
        plt.xlabel('vx')
        plt.ylabel('vz')
        plt.title('SinCos Polar Angle Distribution')
        plt.savefig('plots/vxvz_lab')
        plt.close()


    #########################################
    #make NetCDF Particle Source file
    #########################################

    rootgrp = netCDF4.Dataset("particleSource.nc", "w", format="NETCDF4")
    npp = rootgrp.createDimension("nP", nP)
    xxx = rootgrp.createVariable("x","f8",("nP"))
    yyy = rootgrp.createVariable("y","f8",("nP"))
    zzz = rootgrp.createVariable("z","f8",("nP"))
    vxx = rootgrp.createVariable("vx","f8",("nP"))
    vyy = rootgrp.createVariable("vy","f8",("nP"))
    vzz = rootgrp.createVariable("vz","f8",("nP"))
    xxx[:] = x
    yyy[:] = y
    zzz[:] = z
    vxx[:] = vx
    vyy[:] = vy
    vzz[:] = vz
    rootgrp.close()


def midpoints(nP,pps_weights,adj,slope,Beta, r1,z1):
    #get midpoints of coords
    r_mid = np.zeros(len(r1)-1)
    z_mid = np.zeros(len(z1)-1)
    for i in range(len(r1)-1):
        r_mid[i] = np.average(np.array([r1[i],r1[i+1]]))
        z_mid[i] = np.average(np.array([z1[i],z1[i+1]]))

    x = np.zeros(nP)
    y = np.zeros(nP)
    z = np.zeros(nP)
    counter = 0
    for i in range(len(pps_weights)):
        x[counter:counter+pps_weights[i]] = r_mid[i] - adj*np.abs(np.cos(Beta[i]))
        z[counter:counter+pps_weights[i]] = z_mid[i] + np.sign(slope[i]) * adj*np.abs(np.sin(Beta[i]))
        counter += pps_weights[i]

    return x,y,z


def uniform(nP,pps_weights,adj,slope,Beta, r_coord,z_coord):
    r1 = r_coord[:-1]
    z1 = z_coord[:-1]
    r2 = r_coord[1:]
    z2 = z_coord[1:]
    
    x = np.zeros(nP)
    y = np.zeros(nP)
    z = np.zeros(nP)
    counter = 0
    for i in range(len(pps_weights)):
        dr = (r2[i]-r1[i])/pps_weights[i]
        dz = (z2[i]-z1[i])/pps_weights[i]
        r_segment = np.linspace(r1[i]+dr/2, r2[i]-dr/2, pps_weights[i])
        z_segment = np.linspace(z1[i]+dz/2, z2[i]-dz/2, pps_weights[i])

        x[counter:counter+pps_weights[i]] = r_segment - adj*np.abs(np.cos(Beta[i]))
        z[counter:counter+pps_weights[i]] = z_segment + np.sign(slope[i]) * adj*np.abs(np.sin(Beta[i]))
        counter += pps_weights[i]

    return x,y,z


def random(nP,pps_weights,adj,slope,Beta, r_coord,z_coord):
    r1 = r_coord[:-1]
    z1 = z_coord[:-1]
    r2 = r_coord[1:]
    z2 = z_coord[1:]
    
    x = np.zeros(nP)
    y = np.zeros(nP)
    z = np.zeros(nP)
    counter = 0
    for i in range(len(pps_weights)):
        for j in range(pps_weights[i]):
            chi = np.random.rand(1)
            x[counter+j] = r1[i]+chi*(r2[i]-r1[i]) - adj*np.abs(np.cos(Beta[i]))
            z[counter+j] = z1[i]+chi*(z2[i]-z1[i]) + np.sign(slope[i]) * adj*np.abs(np.sin(Beta[i]))
        counter += pps_weights[i]

    return x,y,z


def interpolate(small,big):
    indices = np.zeros(len(small))
    for i in range(len(small)):
        diff = np.min(np.abs(big-small[i]))
        index_possibilities = np.array([np.where(big==small[i]+diff)[0], \
                             np.where(big==small[i]-diff)[0]], dtype=object)
        try: indices[i] = index_possibilities[1]
        except: indices[i] = index_possibilities[0]

    return indices.astype(int)


def get_surf_profiles(profilesFile, r1, z1, slope, z, plot_variables):
    profiles = netCDF4.Dataset(profilesFile)
    
    #get mesh grid for the plasma profiles used in GITR
    r_mesh = profiles.variables['r'][:]
    z_mesh = profiles.variables['z'][:]
    
    #get midpoints of coords at the target (identical to midpoint func)
    r_wall = np.zeros(len(r1)-1)
    z_wall = np.zeros(len(z1)-1)
    for i in range(len(r1)-1):
        r_wall[i] = np.average(np.array([r1[i],r1[i+1]]))
        z_wall[i] = np.average(np.array([z1[i],z1[i+1]]))
    
    #figure out which indices touch the wall on the profiles mesh
    r_indices = interpolate(r_wall,r_mesh)
    z_indices = interpolate(z_wall,z_mesh)

    #scoot the r_indices 1 cell to the left if the profiles.nc gridcell is too far off
    br = profiles.variables['br'][:][z_indices,r_indices]
    for i in range(len(br)):
        if br[i]==-1: r_indices[i] -= 1

    if plot_variables == 1:
        # plot interpolated SOLPS gridpoints on the boundry
        plt.plot(r_mesh[r_indices],z_mesh[z_indices])
        plt.xlabel('r [m]')
        plt.ylabel('z [m]')
        plt.title('Interpolated SOLPS grid cell centers')
        plt.axis('Scaled')
        plt.savefig('plots/test_rz_indices')

    #extract plasma parameters at the wall indices
    te = profiles.variables['te'][:][z_indices,r_indices]
    ti = profiles.variables['ti'][:][z_indices,r_indices]
    SimpleEnergyEst = 2*ti+3*te*z
    v_para0 = profiles.variables['v_para0'][:][z_indices,r_indices]
    v_para1 = profiles.variables['v_para1'][:][z_indices,r_indices]
    ni0 = profiles.variables['ni0'][:][z_indices,r_indices]
    ni1 = profiles.variables['ni1'][:][z_indices,r_indices]
    br = profiles.variables['br'][:][z_indices,r_indices]
    bt = profiles.variables['bt'][:][z_indices,r_indices]
    bz = profiles.variables['bz'][:][z_indices,r_indices]

    #TODO: will need v, n, and fluxes for all 9 species

    #get normal directions at each wall segment
    #get bmag and bnorm at wall indices from br, bt, bz
    norm_slope = -1/slope
    bmag = np.sqrt(br**2 + bt**2 + bz**2)
    Bangle = np.zeros(len(norm_slope))
    for i in range(len(norm_slope)):
        norm_mag = np.sqrt(1+norm_slope[i]**2)
        norm = np.array([-1/norm_mag, 0, norm_slope[i]/norm_mag])

        bnorm = np.array([br[i], bt[i], bz[i]])/bmag[i]
        Bangle[i] = np.arccos(np.dot(norm,bnorm))*180/np.pi

    if plot_variables == 1:
        #plot plasma parameters along the surface
        plt.close()
        plt.plot(z_mesh[z_indices], te)
        plt.xlabel('z [m]')
        plt.ylabel('Te [eV]')
        plt.title('Electron Temperature along the SAS-V Divertor')
        plt.savefig('plots/surf_te.png')

        plt.close()
        plt.plot(z_mesh[z_indices], SimpleEnergyEst)
        plt.xlabel('z [m]')
        plt.ylabel('energyy [eV]')
        plt.title('Estimate of incoming D,C Ion Energy along the SAS-V Divertor')
        plt.savefig('plots/surf_energyest')

        plt.close()
        plt.plot(z_mesh[z_indices], Bangle)
        plt.xlabel('z [m]')
        plt.ylabel('angle [degrees]')
        plt.title('Angle of B Field with respect to the surface normal')
        plt.savefig('plots/surf_Bangle.png')
        plt.close()

    return SimpleEnergyEst, Bangle


def get_ft_spyld(surfE, surfA, ftBFile):
    #import sputtering yield tables for D0 and D1+ on W
    ftB = netCDF4.Dataset(ftBFile, "r", format="NETCDF4")
    spyld = ftB.variables['spyld'][:]
    ftE = ftB.variables['E'][:]
    ftA = ftB.variables['A'][:]

    #find which FT index contains the energies and angles at the surface
    E_indices = interpolate(surfE, ftE)
    A_indices = interpolate(surfA, ftA)

    surfY = spyld[:,E_indices,A_indices]
    return surfY


def FittingParameters_NonW(E, Eth):
    Esp = np.ones(len(E))
    f = np.zeros(len(E))
    b = np.zeros(len(E))
    c = np.zeros(len(E))
    for i,v in enumerate(E):
        if v<Eth:
            continue
        elif v>=Eth and v<49:
            f[i] = 2.9557
            b[i] = 5.8879
            c[i] = 0.9465
        elif v>=65 and v<75:
            f[i] = 0.4622
            b[i] = 2.5095
            c[i] = 1.0118
        elif v>=85 and v<95:
            f[i] = 2.1152
            b[i] = 2.6541
            c[i] = 0.9226
        elif v>=170 and v<250:
            f[i] = 2.0138
            b[i] = 1.3460
            c[i] = 1.0316
        elif v>=250 and v<400:
            f[i] = 2.2531
            b[i] = 1.2151
            c[i] = 1.0310
        else:
            print('WARNING: TABULAR SPUTTERING FITTING PARAMETERS MISSING')
    return Esp, f, b, c


def get_analytic_spyld(surfE, surfA, Z1=6, M1=12, Z2=74, M2=183.84, FitParam='Ni', Eth=45.3362, lam=0.0921, q=1.4389, mu=2.0225):
    #this entire function that is simply the Eckstein formula
    #defaults are for C on W but with fitting parameters, Eth, and Esp for N on W
    #M1, Z1 are for the projectile
    #M2, Z2 are for the target
    #Eth is the threshold energy for any sputtering to occur
    #Esp is the SBE for self-bombardment
    #lam, q, mu, f, b, c are fitting parameters from Eckstein tables

    e = -1.602176634e-19

    a_L = 0.0529177*((9*np.pi**2/128)**(1/3)) + \
        (Z1**(2/3) + Z2**(2/3))**-0.5

    epsilon_L = surfE * (M2/(M1+M2)) * (a_L/(Z1*Z2*e**2))

    omega = epsilon_L + 0.1728*np.sqrt(epsilon_L) + \
        0.008*epsilon_L**0.1504

    snKrC = 0.5*np.log(1+1.2288*epsilon_L)/omega
    
    #spyld(E, angle = 0) normal incidence
    Y_0 = q*snKrC*(surfE/Eth-1)**mu / \
        (lam/omega + (surfE/Eth-1)**mu)
    
    Y_0[np.where(surfE<Eth)[0]] = 0

    #choose basis for fitting parameters for different energies
    Esp, f, b, c = FittingParameters_NonW(surfE,Eth)

    theta_0 = np.pi - np.arccos(np.sqrt(1/(1+surfE/Esp)))

    #spyld(E,A)    
    Y = Y_0 * (np.cos((np.pi*surfA/(2*theta_0))**c))**(-f) * \
        np.exp(b*(1-1/np.cos((np.pi*surfA/(2*theta_0))**c)))

    return Y_0








