import sys, os
sys.path.insert(0, os.path.abspath('../input/'))
sys.path.insert(0, os.path.abspath('../../../python/'))
sys.path.insert(0, os.path.abspath('../../../../pyGITR/pyGITR'))

import numpy as np
import scipy.interpolate as scii
import matplotlib.pyplot as plt
import math

import gitr
import solps
import Particles
import netCDF4

def init():
    #set plotting style defaults
    plt.rcParams.update({'font.size':11.5})
    plt.rcParams.update({'lines.linewidth':1.2})
    plt.rcParams.update({'lines.markersize':1})

def line_source(nP, \
                bField = '../input/bField.nc', \
                part_type = 'helium', \
                targFile = 'assets/rightTargOutput', \
                wallFile = 'assets/gitr_rz.txt', \
                surf_coarse = 'assets/surf_coarse.txt', \
                surf_ind = 'assets/surf_ind.txt', \
                configuration = 'random', \
                plot_variables = 0, \
                r_W = None, z_W = None, rCoarse = None, zCoarse = None, addedPoints = 0):
    
    # defined here since they depend on particle type
    geom = '../'+part_type+'/input/gitrGeometry.cfg'
    profilesFile = '../'+part_type+'/input/profiles.nc'
    
    # create random points on r axis for particles
    x = np.empty(nP)
    for i in range(nP):
        chi = np.random.rand(1)
        x[i] = 1.45 + 0.04*chi
    lineheight = 1.16
    
    # Import b field info
    bProfile = netCDF4.Dataset(bField)
    br = bProfile.variables['br'][:]
    bt = bProfile.variables['bt'][:]
    bz = bProfile.variables['bz'][:]
    r = bProfile.variables['r'][:]
    z = bProfile.variables['z'][:]
    
    # Find index where z exists
    acc = 100
    zloc = np.where((np.round(z*acc))/acc==lineheight)[0][0]
    
    # Calculate B field and velocities at each point with a particle
    B = np.empty(nP)
    for i in range(nP):
        xloc = np.where((np.round(r*acc))/acc==(np.round(x[i]*acc)/acc))[0][0]
        brloc = br[zloc,xloc]
        bzloc = bz[zloc,xloc]
        B[i] = (brloc**2+bzloc**2)**0.5
    
    
    # creates y and z value
    z = lineheight*np.ones(nP)
    y = np.zeros(nP)
    vx = np.zeros(nP)
    vy = 500*np.zeros(nP)
    vz = 300*np.ones(nP)

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
    
def conical_source(nP):
    
    # create random points on r axis for particles
    x = np.empty(nP)
    y = np.empty(nP)
    z = np.empty(nP)
    
    # adjusts closeness to wall
    adj = 0.007
    
    # Creates round surface as a function of y and z position
    for i in range(nP):
        chiy = np.random.rand(1)
        chiz = np.random.rand(1)
        y[i] = 0.05*chiy
        z[i] = 1.12 + 0.05*chiz
        x[i] = math.sin(y[i]*math.pi/.05)*math.sin((z[i]-1.12)*math.pi/0.05)/25+1.45+adj
        y[i] = y[i]-0.025
        
    # creates y and z value
    vx = 200*np.ones(nP)
    vy = 50*np.ones(nP)
    vz = 50*np.ones(nP)

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
    
def point_source(nP,r,z):
    x = r*np.ones(nP)
    y = np.zeros(nP)
    z = z*np.ones(nP)
    vx = 500*np.ones(nP)
    vy = 500*np.ones(nP)
    vz = 300*np.ones(nP)

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




def get_incoming_IEADs(q, profiles, surfW, rmrsCoarse, rmrsFine):
    #extract plasma parameters at the wall indices
    teCoarse = profiles.variables['te_inner_target'][:]
    tiCoarse = profiles.variables['ti_inner_target'][:]
    
    #get temp as a function of rmrsCoarse
    fte = scii.interp1d(rmrsCoarse, teCoarse, fill_value='extrapolate')
    fti = scii.interp1d(rmrsCoarse, tiCoarse, fill_value='extrapolate')
    
    #interpolate temp at points in rmrsFine
    te = fte(rmrsFine)
    ti = fti(rmrsFine)

    SimpleEnergyEst = 2*ti+3*te*q

    if q<0:
        Esp, f, b, c, ThetaMax = FittingParameters_NonW(SimpleEnergyEst)
        AngleEst = ThetaMax
    else: AngleEst = 70*np.ones(len(SimpleEnergyEst))

    return SimpleEnergyEst, AngleEst

def midpoints(nP,pps_weights,adj,slope,Beta, r1,z1,r2,z2):
    #get midpoints of coords
    r_mid = np.zeros(len(r1)-1)
    z_mid = np.zeros(len(z1)-1)
    for i in range(len(r1)):
        r_mid[i] = np.average(np.array([r1[i],r2[i]]))
        z_mid[i] = np.average(np.array([z1[i],z2[i]]))
        
    x = np.zeros(nP)
    y = np.zeros(nP)
    z = np.zeros(nP)
    counter = 0
    for i in range(len(pps_weights)):
        x[counter:counter+pps_weights[i]] = r_mid[i] - adj*np.abs(np.cos(Beta[i]))
        z[counter:counter+pps_weights[i]] = z_mid[i] + np.sign(slope[i]) * adj*np.abs(np.sin(Beta[i]))
        counter += pps_weights[i]

    return x,y,z


def uniform(nP,pps_weights,adj,slope,Beta, r1,z1,r2,z2):
    x = np.zeros(nP)
    y = np.zeros(nP)
    z = np.zeros(nP)
    tally = 0
    for i in range(len(pps_weights)):
        if pps_weights[i]>0:
            dr = (r2[i]-r1[i])/pps_weights[i]
            dz = (z2[i]-z1[i])/pps_weights[i]
            r_segment = np.linspace(r1[i]+dr/2, r2[i]-dr/2, pps_weights[i])
            z_segment = np.linspace(z1[i]+dz/2, z2[i]-dz/2, pps_weights[i])
            
            x[tally:tally+pps_weights[i]] = r_segment - adj*np.abs(np.cos(Beta[i]))
            z[tally:tally+pps_weights[i]] = z_segment + np.sign(slope[i]) * adj*np.abs(np.sin(Beta[i]))
            tally += pps_weights[i]

    return x,y,z


def random(nP,pps_weights,adj,slope,Beta, r1,z1,r2,z2):
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


def refine_linear(a,b,coarse_rmrs, fine_rmrs, coarse_variable):
    fine_variable = np.zeros(len(fine_rmrs))
    for i in range(len(fine_rmrs)):
        #find the smallest coarse index near the fine index
        index = np.where(coarse_rmrs<=fine_rmrs[i])
        if index[0].any() > 0:
            lower_bound_index = np.max(index)
        else: lower_bound_index = 0
        
        if lower_bound_index < len(coarse_rmrs)-1:
            upper_bound_index = lower_bound_index + 1 
            lower_bound_rmrs = coarse_rmrs[lower_bound_index]
            upper_bound_rmrs = coarse_rmrs[upper_bound_index]
            lower_bound_variable = coarse_variable[lower_bound_index]
            upper_bound_variable = coarse_variable[upper_bound_index]
            
            fraction = (fine_rmrs[i]-lower_bound_rmrs) / (upper_bound_rmrs-lower_bound_rmrs)
            fine_variable[i] = fraction * (upper_bound_variable-lower_bound_variable) + lower_bound_variable
        else: fine_variable[i] = coarse_variable[-1]
        
    return fine_variable


def get_surf_profiles(profilesFile, r1, z1, r2, z2, q, plot_variables):
    profiles = netCDF4.Dataset(profilesFile)
    
    #get mesh grid for the plasma profiles used in GITR
    r_mesh = profiles.variables['r'][:]
    z_mesh = profiles.variables['z'][:]

    #get midpoints for the refined W surface        
    r_mid = np.average(np.array([r1,r2]),axis=0)
    z_mid = np.average(np.array([z1,z2]),axis=0)
        
    #figure out which indices touch the wall on the profiles mesh
    r_indices = interpolate(r_mid,r_mesh)
    z_indices = interpolate(z_mid,z_mesh)
        
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
    '''
    te = profiles.variables['te'][:]
    ti = profiles.variables['ti'][:]
    print('meshes',np.shape(np.array([z_mesh,r_mesh])))
    print('te',np.shape(te))
    print('mids',np.shape(np.transpose([r_mid,z_mid])))
    te_interp = np.zeros(len(r_mid))
    for i in range(len(r_mid)):
        te_interp[i] = scii.interpn([z_mesh,r_mesh],te,[z_mid[i],r_mid[i]])
    print('test te_interp',np.shape(te_interp))
    '''
    
    SimpleEnergyEst = 2*ti+3*te*q

    if q==7 or q==6:
        Esp, f, b, c, ThetaMax = FittingParameters_NonW(SimpleEnergyEst)
        AngleEst = ThetaMax
    else: AngleEst = (np.pi/2)*np.ones(len(SimpleEnergyEst))

    if plot_variables == 1:
        dist = np.sqrt((r1-r2)**2 + (z1-z2)**2)
        rmrs = np.append(np.zeros(1), np.cumsum(dist))
        rmrs_mid = np.average(np.array([rmrs[:-1],rmrs[1:]]),axis=0)
        
        ls, fs = 5, 14
        #plot plasma parameters along the surface
        plt.close()
        plt.plot(rmrs_mid, te, linewidth=ls)
        plt.xlabel('D-Dsep [m]',fontsize=fs)
        plt.ylabel('Te [eV]',fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.title('Electron Temperature along the SAS-V Divertor',fontsize=fs)
        plt.savefig('plots/surf_te.png')

        plt.close()
        plt.plot(rmrs_mid, AngleEst)
        plt.xlabel('z [m]')
        plt.ylabel('angle [degrees]')
        plt.title('Angle of B Field with respect to the surface normal')
        plt.savefig('plots/surf_Bangle.png')
        plt.close()

    return rmrs_mid, SimpleEnergyEst, AngleEst

# Fractal Trinyd Surface Sputtering Yeild
def get_ft_spyld(S, surfE, surfA, ftBFile):
    #import sputtering yield tables for incident ions on W
    ftB = netCDF4.Dataset(ftBFile, "r", format="NETCDF4")
    spyld = ftB.variables['spyld'][S][:]
    ftE = ftB.variables['E'][:]
    ftA = ftB.variables['A'][:]
    
    surfY = scii.interpn((ftE,ftA), spyld, (surfE,surfA))
    
    return surfY


def FittingParameters_NonW(E):
    Eth=45.3362
    Esp = np.ones(len(E))
    f = np.zeros(len(E))
    b = np.zeros(len(E))
    c = np.zeros(len(E))
    ThetaMax = np.zeros(len(E))
    for i,v in enumerate(E):
        if v<Eth:
            continue
        elif v>=Eth and v<49:
            f[i],b[i],c[i] = 2.9557, 5.8879, 0.9465
        elif v>=49 and v<51:
            f[i],b[i],c[i] = 1.7735, 4.3144, 0.9468
        elif v>=51 and v<53.5:
            f[i],b[i],c[i] = 1.2707, 3.6458, 0.8840
        elif v>=53.5 and v<57.5:
            f[i],b[i],c[i] = 1.1002, 3.3751, 0.9604
        elif v>=57.5 and v<65:
            f[i],b[i],c[i] = 0.4622, 2.5095, 1.0118
        elif v>=65 and v<75:
            f[i],b[i],c[i] = 0.4622, 2.5095, 1.0118
        elif v>=75 and v<85:
            f[i],b[i],c[i] = 2.7960, 3.4029, 0.8841
        elif v>=85 and v<95:
            f[i],b[i],c[i] = 2.1152, 2.6541, 0.9226
        elif v>=95 and v<110:
            f[i],b[i],c[i] = 1.7312, 2.1735, 0.9489
        elif v>=110 and v<130:
            f[i],b[i],c[i] = 1.6230, 1.6737, 1.0004
        elif v>=130 and v<170:
            f[i],b[i],c[i] = 1.7195, 1.5092, 1.0176
            ThetaMax[i] = 30.54
        elif v>=170 and v<250:
            f[i],b[i],c[i] = 2.0138, 1.3460, 1.0316
            ThetaMax[i] = 51.98
        elif v>=250 and v<400:
            f[i],b[i],c[i] = 2.2531, 1.2151, 1.0310
            ThetaMax[i] = 59.47
        elif v>=400 and v<750:
            f[i],b[i],c[i] = 2.4324, 1.1313, 1.0171
            ThetaMax[i] = 62.20
        elif v>=750 and v<2000:
            f[i],b[i],c[i] = 2.4383, 0.9940, 0.9936
            ThetaMax[i] = 66.00
        else:
            print('WARNING: TABULAR SPUTTERING FITTING PARAMETERS MISSING')
    return Esp, f, b, c, ThetaMax


def get_analytic_spyld(surfE, surfA, Z1=6, M1=12, Z2=74, M2=183.84, FitParam='N', Eth=45.3362, lam=0.0921, q=1.4389, mu=2.0225):
    #this entire function that is simply the Eckstein formula
    #defaults are for C on W but with fitting parameters, Eth, and Esp for N on W
    #M1, Z1 are for the projectile
    #M2, Z2 are for the target
    #Eth is the threshold energy for any sputtering to occur
    #Esp is the SBE for self-bombardment
    #lam, q, mu, f, b, c are fitting parameters from Eckstein tables

    e = 14.399651

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
    Esp, f, b, c, ThetaMax = FittingParameters_NonW(surfE)

    theta_0 = np.pi - np.arccos(np.sqrt(1/(1+surfE/Esp)))

    #spyld(E,A)    
    Y = Y_0 * (np.cos((np.pi*surfA/(2*theta_0))**c))**(-f) * \
        np.exp(b*(1-1/np.cos((np.pi*surfA/(2*theta_0))**c)))

    return Y_0

if __name__ == "__main__":
    init()
    

nP=1e5
geom = '../input/gitrGeometry.cfg'
profiles_file = '../input/plasmaProfiles.nc'
gitr_rz = 'assets/gitr_rz.txt'
rmrs_fine_file = 'assets/rmrs_fine.txt'
W_fine_file = 'assets/W_fine.txt'
ftDFile = 'assets/ftridynBackgroundD.nc'
ftCFile = 'assets/ftridynBackgroundC.nc'
configuration = 'random'
plot_variables = 1
    
#import wall geometry to plot over
with open(gitr_rz, 'r') as file:
    wall = file.readlines()

    
#import W surface indices
with open(W_fine_file, 'r') as file:
    W_fine = file.readlines()
W_fine = np.array(W_fine,dtype='int')

#import coarse rmrs at W surface
profiles = netCDF4.Dataset(profiles_file)
r_right_target = profiles.variables['r_inner_target'][:]
z_right_target = profiles.variables['z_inner_target'][:]
rmrsCoarse_in = profiles.variables['rmrs_inner_target'][:]
r_left_target = profiles.variables['r_inner_target'][:]
z_left_target = profiles.variables['z_inner_target'][:]
rmrsCoarse_out = profiles.variables['rmrs_outer_target'][:]
surfW_left = np.arange(99,99+len(rmrsCoarse_out))
surfW_right = np.arange(137,137+len(rmrsCoarse_in))

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

if plot_variables==1:
    plt.close()
    plt.plot(r_right_target, z_right_target, '-r', linewidth=0.5, zorder=0)
    # plt.plot(r_right_target[surfW], z_right_target[surfW], 'purple', label='Profiles Tungsten', linewidth=3, zorder=1)
    plt.plot(r_left_target, z_left_target, '-b', linewidth=0.5, zorder=0)
    plt.plot(R, Z, 'violet', label='Tungsten', linewidth=0.6, zorder=2)
    plt.scatter(R, Z, marker='_', color='violet', s=8, zorder=3)
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Divertors in West \n makeParticleSource')
    plt.savefig('plots/makePSGeom.png')

slope = np.zeros(len(r1))
Alpha = np.zeros(len(r1))
for i in range(len(r1)):
    if (r2[i]-r1[i])!=0:
        slope[i] = (z2[i]-z1[i])/(r2[i]-r1[i])
        #define angle between material wall and major radius, x
        Alpha[i] = np.abs(np.arctan((z2[i]-z1[i]) / (r2[i]-r1[i])))
    elif (r2[i]-r1[i])==0:
        slope[i] = 100
        Alpha[i] = 89.999*np.pi/180

Alpha = np.abs(Alpha) #np.abs(np.rad2deg(Alpha[W_fine[:-1]]))
Beta = np.abs(np.pi/2 - Alpha)

dist = np.sqrt(np.power(r1-r2,2) + np.power(z1-z2,2))
area = np.pi*(r1+r2)*dist

##############################################
#get W/s sputtered by D, C flux to wall
##############################################
#%%
#get incoming ion energy and angle estimations where the integer input is z
#using outer divertor since sputtering is greatest there
surfW_left = surfW_left[:-1]
rmrsCoarse_out = rmrsCoarse_out[:-1]
energyD, angleD = get_incoming_IEADs(1, profiles, surfW_left, rmrsCoarse_out, rmrsFine)
energyC1, angleC1 = get_incoming_IEADs(1, profiles, surfW_left, rmrsCoarse_out, rmrsFine)
energyC2, angleC2 = get_incoming_IEADs(2, profiles, surfW_left, rmrsCoarse_out, rmrsFine)
energyC3, angleC3 = get_incoming_IEADs(3, profiles, surfW_left, rmrsCoarse_out, rmrsFine)
energyC4, angleC4 = get_incoming_IEADs(4, profiles, surfW_left, rmrsCoarse_out, rmrsFine)
energyC5, angleC5 = get_incoming_IEADs(5, profiles, surfW_left, rmrsCoarse_out, rmrsFine)
energyC6, angleC6 = get_incoming_IEADs(6, profiles, surfW_left, rmrsCoarse_out, rmrsFine)

if plot_variables == 1:
    plt.close()
    plt.plot(rmrsFine, energyD, 'black', label='D1+')
    plt.plot(rmrsFine, energyC1, 'red', label='C1+')
    plt.plot(rmrsFine, energyC2, 'darkorange', label='C2+')
    plt.plot(rmrsFine, energyC3, 'gold', label='C3+')
    plt.plot(rmrsFine, energyC4, 'green', label='C4+')
    plt.plot(rmrsFine, energyC5, 'blue', label='C5+')
    plt.plot(rmrsFine, energyC6, 'purple', label='C6+')
    plt.plot(rmrsFine, 45.3362*np.ones(len(rmrsFine)), 'gray', label='W Eth')
    plt.legend()
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('energy [eV]')
    plt.title('Estimate of incident background ion energies')
    plt.savefig('plots/particle-source/incident_energy')
    
    plt.close()
    plt.plot(rmrsFine, angleD, 'black', label='D1+')
    plt.plot(rmrsFine, angleC1, 'red', label='C1+')
    plt.plot(rmrsFine, angleC2, 'darkorange', label='C2+')
    plt.plot(rmrsFine, angleC3, 'gold', label='C3+')
    plt.plot(rmrsFine, angleC4, 'green', label='C4+')
    plt.plot(rmrsFine, angleC5, 'blue', label='C5+')
    plt.plot(rmrsFine, angleC6, 'purple', label='C6+')        
    plt.legend()
    plt.xlabel('z [m]')
    plt.ylabel('angle [degrees]')
    plt.title('Used angles of incidence for incident background ions')
    plt.savefig('plots/particle-source/incident_angle.png')
    plt.close()

#get sputtering yields for D0 and D1+ on W from fractal tridyn tables
#Cspyld = get_ft_spyld(CsurfE, CsurfA, ftCFile)[0]
spyldD = get_ft_spyld(1, energyD, angleD, ftDFile) #file input includes He, which we aren't using
spyldC1 = get_ft_spyld(0, energyC1, angleC1, ftCFile)
spyldC2 = get_ft_spyld(0, energyC2, angleC2, ftCFile)
spyldC3 = get_ft_spyld(0, energyC3, angleC3, ftCFile)
spyldC4 = get_ft_spyld(0, energyC4, angleC4, ftCFile)
spyldC5 = get_ft_spyld(0, energyC5, angleC5, ftCFile)
spyldC6 = get_ft_spyld(0, energyC6, angleC6, ftCFile)

#get coarse flux profile from background D, C and refine to rmrsFine
fluxCoarseD = np.abs(profiles.variables['flux_inner_target'][1][surfW_left])
fluxCoarseC1 = np.abs(profiles.variables['flux_inner_target'][3][surfW_left])
fluxCoarseC2 = np.abs(profiles.variables['flux_inner_target'][4][surfW_left])
fluxCoarseC3 = np.abs(profiles.variables['flux_inner_target'][5][surfW_left])
fluxCoarseC4 = np.abs(profiles.variables['flux_inner_target'][6][surfW_left])
fluxCoarseC5 = np.abs(profiles.variables['flux_inner_target'][7][surfW_left])
fluxCoarseC6 = np.abs(profiles.variables['flux_inner_target'][8][surfW_left])

ffluxD = scii.interp1d(rmrsCoarse_out,fluxCoarseD,fill_value='extrapolate')
ffluxC1 = scii.interp1d(rmrsCoarse_out,fluxCoarseC1,fill_value='extrapolate')
ffluxC2 = scii.interp1d(rmrsCoarse_out,fluxCoarseC2,fill_value='extrapolate')
ffluxC3 = scii.interp1d(rmrsCoarse_out,fluxCoarseC3,fill_value='extrapolate')
ffluxC4 = scii.interp1d(rmrsCoarse_out,fluxCoarseC4,fill_value='extrapolate')
ffluxC5 = scii.interp1d(rmrsCoarse_out,fluxCoarseC5,fill_value='extrapolate')
ffluxC6 = scii.interp1d(rmrsCoarse_out,fluxCoarseC6,fill_value='extrapolate')

fluxD = ffluxD(rmrsFine)
fluxC1 = ffluxC1(rmrsFine)
fluxC2 = ffluxC2(rmrsFine)
fluxC3 = ffluxC3(rmrsFine)
fluxC4 = ffluxC4(rmrsFine)
fluxC5 = ffluxC5(rmrsFine)
fluxC6 = ffluxC6(rmrsFine)

#multiply incoming ion flux by Y_s to get sputtered W flux by each species
sputt_fluxD = spyldD*fluxD
sputt_fluxC1 = spyldC1*fluxC1
sputt_fluxC2 = spyldC2*fluxC2
sputt_fluxC3 = spyldC3*fluxC3
sputt_fluxC4 = spyldC4*fluxC4
sputt_fluxC5 = spyldC5*fluxC5
sputt_fluxC6 = spyldC6*fluxC6
sputt_flux = sputt_fluxD + sputt_fluxC1 + sputt_fluxC2 + sputt_fluxC3 + sputt_fluxC4 + sputt_fluxC5 + sputt_fluxC6
print('SPUTT FLUX',len(sputt_flux),'\n',sputt_flux)

#multiply by area to get the outgoing particles per second
pps = np.multiply(sputt_flux,area)
pps_weights = nP*pps/np.sum(pps)

if plot_variables == 1:        
    plt.close()
    plt.plot(rmrsFine, fluxD, 'black', label='D1+')
    plt.plot(rmrsFine, fluxC1, 'red', label='C1+')
    plt.plot(rmrsFine, fluxC2, 'darkorange', label='C2+')
    plt.plot(rmrsFine, fluxC3, 'gold', label='C3+')
    plt.plot(rmrsFine, fluxC4, 'green', label='C4+')
    plt.plot(rmrsFine, fluxC5, 'blue', label='C5+')
    plt.plot(rmrsFine, fluxC6, 'purple', label='C6+')
    plt.yscale('log')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [#/m2s]')
    plt.legend(loc='upper right')
    plt.title('Incident Ion Flux')
    plt.savefig('plots/particle-source/incident_flux.png')

    plt.close()
    plt.plot(rmrsFine, spyldD, 'black', label='D1+')
    plt.plot(rmrsFine, spyldC1, 'red', label='C1+')
    plt.plot(rmrsFine, spyldC2, 'darkorange', label='C2+')
    plt.plot(rmrsFine, spyldC3, 'gold', label='C3+')
    plt.plot(rmrsFine, spyldC4, 'green', label='C4+')
    plt.plot(rmrsFine, spyldC5, 'blue', label='C5+')
    plt.plot(rmrsFine, spyldC6, 'purple', label='C6+')
    plt.legend()
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Sputtering Yield')
    plt.title('W Sputtering Yield by Incident D and C')
    plt.savefig('plots/particle-source/spyld.png')
    
    plt.close()
    plt.plot(rmrsFine, sputt_fluxD, 'black', label='D1+')
    plt.plot(rmrsFine, sputt_fluxC1, 'red', label='C1+')
    plt.plot(rmrsFine, sputt_fluxC2, 'darkorange', label='C2+')
    plt.plot(rmrsFine, sputt_fluxC3, 'gold', label='C3+')
    plt.plot(rmrsFine, sputt_fluxC4, 'green', label='C4+')
    plt.plot(rmrsFine, sputt_fluxC5, 'blue', label='C5+')
    plt.plot(rmrsFine, sputt_fluxC6, 'purple', label='C6+')
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [#/m2s]')
    plt.legend(loc='upper left')
    plt.title('Flux of W Sputtered off Wall')
    plt.savefig('plots/particle-source/sputt_flux_charge_dependent.png')

    plt.close()
    plt.plot(rmrsFine, sputt_flux)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Flux [#/m2s]')
    plt.title('Flux of W Sputtered off Wall')
    plt.savefig('plots/particle-source/sputt_flux.png') 

    plt.close()
    area_cm = area*100*100
    plt.plot(rmrsFine, area_cm)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Area [cm2]')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(1,5))        
    plt.title('Toroidal surface area of each line segment')
    plt.savefig('plots/particle-source/area.png')
    
    plt.close()
    plt.plot(rmrsFine, pps_weights)
    plt.xlabel('D-Dsep [m]')
    plt.ylabel('Weighted PPS [#/s]')
    plt.title('Computational Weights for Initially Sputtered W/s \n nP = '+str(nP))
    plt.savefig('plots/particle-source/pps_weights.png')

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
    while nP_diff<0:
        rand_index = np.random.choice(len(pps_weights))
        if pps_weights[rand_index] - 1 > 0:
            pps_weights[rand_index] -= 1
            nP_diff+=1

pps_weights = np.round(pps_weights)
int_weights = np.array(pps_weights,dtype='int')
nP_diff = nP-np.sum(int_weights)

print('total nP', nP)
print('pps over nP', np.sum(pps)/nP)
print('nP(r_mid):', int_weights)
print('nP_diff should be 0: ', nP_diff)


#define adjustment into the sheath because particles can't start exactly on the wall
adj = 1e-7

#populate x,y,z with r_mid,0,z_mid
if configuration == 'random': 
    x,y,z = random(nP,int_weights,adj,slope,Beta, r1,z1,r2,z2)
elif configuration == 'uniform': 
    x,y,z = uniform(nP,int_weights,adj,slope,Beta, r1,z1,r2,z2)
elif configuration == 'midpoint': 
    x,y,z = midpoints(nP,int_weights, adj,slope,Beta, r1,z1,r2,z2)
else:
    print('(x,y,z) configuration not set')

#########################################
#get vx,vy,vz from IEADs
#########################################

#use PyGITR to set up vx,vy,vz,E,theta,psi distributions
PartDist = Particles.ParticleDistribution(nP, ListAttr=['vx','vy','vz'])

vx = np.zeros(1)
vy = np.zeros(1)
vz = np.zeros(1)

for i in range(len(int_weights)):

    weight = int(int_weights[i])

    if weight>0:
        m = np.sign(slope[i])
        #get IEADs for sputtered W
        E = PartDist.Generate(weight, 'Thomson')
        PolAng = PartDist.Generate(weight, 'SinCos', x=np.linspace(0,np.pi/2,10*weight))
        AziAng = PartDist.Generate(weight, 'Uniform', x=np.linspace(0,2*np.pi,10*weight))
        
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
    plt.hist(E,bins=50)
    plt.xlabel('Energy Bins [eV]')
    plt.ylabel('Histogram')
    plt.title('Thomson Energy Distribution \n nP='+str(int_weights[-1]))
    plt.savefig('plots/particle-source/thomsonDist.png')
    
    #plot particle framed v_dist relations
    plt.close()
    plt.scatter(vx_prime,vy_prime,s=1)
    plt.axis('Scaled')
    plt.xlabel('vx')
    plt.ylabel('vy')
    plt.title('Uniform Azimuthal Angle in the Particle Frame \n nP='+str(int_weights[-1]))
    plt.savefig('plots/particle-source/vxvy_prime.png')
    
    plt.close()
    plt.scatter(vx_prime,vz_prime,s=1)
    plt.axis('Scaled')
    plt.xlabel('vx')
    plt.ylabel('vz')
    plt.title('SinCos Polar Angle in the Particle Frame \n nP='+str(int_weights[-1]))
    plt.savefig('plots/particle-source/vxvz_prime.png')
    
    plt.close()
    plt.scatter(vx_lab,vz_lab,s=1)
    plt.axis('Scaled')
    plt.xlabel('vx')
    plt.ylabel('vz')
    plt.title('SinCos Polar Angle Distribution in \n the Lab Frame nP='+str(int_weights[-1]))
    plt.savefig('plots/particle-source/vxvz_lab.png')
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







