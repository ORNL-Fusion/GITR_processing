import sys
sys.path.insert(0, '../../../python/')
sys.path.insert(0, '../../../../pyGITR/pyGITR')

import numpy as np
import matplotlib.pyplot as plt

import netCDF4
import gitr
import solps
import Particles

def simple2D(nP = int(1e3), \
            geom = '../input/gitrGeometry.cfg', \
            targFile = 'assets/rightTargOutput', \
            coordsFile = 'assets/right_target_coordinates.txt', \
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
    area = np.pi*(r1+r2)*np.sqrt(np.power(r1-r2,2) + np.power(z1 - z2,2))

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
    #get x,y,z distributions for sputtered W
    #########################################

    #get flux from background D
    r_mid, z_mid, ti, ni, flux, te, ne = solps.read_target_file(targFile)

    #calcualte erosion flux
    ion_flux = np.abs(flux[:,1:][W_ind[:-1]]) #only bother with W surfaces
    spyld = 0.1 #assume sputtering yield of 0.1 for everything because we're lazy
    sputt_flux = spyld*ion_flux #multiply incoming ion flux by Y_s to get sputtered W flux
    sputt_flux_total = np.sum(sputt_flux,axis=1) #add together sputtered flux from 8 ion species
    pps = 0.1*np.multiply(sputt_flux_total,area[W_ind[:-1]]) #multiply by area to get the outgoing particles per second
    pps_weights = nP*pps/np.sum(pps)

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
    #x,y,z = random(nP,pps_weights,adj,slope,Beta, r1,z1)
    #x,y,z = uniform(nP,pps_weights,adj,slope,Beta, r1,z1)
    x,y,z = midpoints(nP,pps_weights, adj,slope,Beta, r1,z1)

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
        print('E',E)
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


def check_symmetry(nP = int(1e3)):
    x = 1.75*np.ones(nP)
    y = np.zeros(nP)
    z = 0.4*np.ones(nP)
    vx = np.zeros(nP)
    vy = 1500*np.ones(nP)
    vz = np.zeros(nP)

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



def old_stuff_stolen_from_west_ex(nP = int(1e3), \
            geom = 'gitrGeometry.cfg', \
            targFile = 'assets/rightTargOutput', \
            coordsFile = 'assets/right_target_coordinates.txt'):
    x1, x2, z1, z2, length, Z, slope, inDir = gitr.plot2dGeom(geom)
    coords = np.loadtxt(coordsFile, dtype='float',skiprows=1,delimiter=' ')
    gitr_inds = [int(i) for i in coords[0:-1,3]]
    print('gitr_inds',gitr_inds)
    r1 = x1[gitr_inds] #coords[:,4]
    r2 = x2[gitr_inds] #coords[:,5]
    z1 = z1[gitr_inds] #coords[:,6]
    z2 = z2[gitr_inds] #coords[:,7]
    slope = slope[gitr_inds] #coords[:,7]
    inDir = inDir[gitr_inds] #coords[:,7]
    length = length[gitr_inds] #coords[:,7]
    area = np.pi*(r1+r2)*np.sqrt(np.power(r1-r2,2) + np.power(z1 - z2,2))

    r, z, ti, ni, flux, te, ne = solps.read_target_file(targFile)

    #calcualte erosion flux
    ion_flux = flux[:,1:]
    yf = np.hstack((0.1+0*ion_flux,ion_flux))
    nSpec = int(0.5*yf.shape[1])
    print('nSpec',nSpec)

    sputt_flux = np.zeros((yf.shape[0],nSpec+1))

    for i in range(nSpec):
        sputt_flux[:,i] = np.multiply(yf[:,i],yf[:,i+nSpec])
        sputt_flux[:,-1] = sputt_flux[:,-1] + sputt_flux[:,i]

    buff = np.array([[0,0,0]])
    sputt_flux = np.vstack((buff,sputt_flux,buff))
    print('len sputt flux', len(sputt_flux[:,-1]))

    particles_per_second = np.multiply(np.abs(sputt_flux[:,-1]),area)

    print('sputt flux', sputt_flux)
    print('area',area)
    print('pps',particles_per_second)
    pps_cdf = np.cumsum(particles_per_second)
    pps_sum = pps_cdf[-1]
    pps_cdf = pps_cdf/pps_cdf[-1]
    pps_cdf = np.transpose(np.matlib.repmat(pps_cdf,nP,1))
    rand1 = np.random.rand(nP)
    rand1 = np.matlib.repmat(rand1,sputt_flux.shape[0],1)

    print('rand1 shape',rand1.shape)

    diff = pps_cdf - rand1
    diff[diff<0.0] = 100.0
    mins = np.argmin(diff,axis=0)

    print('len mins',len(mins))
    print('mins',mins)
    print('diff',diff[mins,range(nP)])
    print('scale',pps_cdf[-1]/particles_per_second[mins])

    r_sample = r1[mins] + (r2[mins] - r1[mins])*diff[mins,range(nP)]/particles_per_second[mins]*pps_sum
    z_sample = z1[mins] + (z2[mins] - z1[mins])*diff[mins,range(nP)]/particles_per_second[mins]*pps_sum
    #print('rsamp',r_sample)
    #print('zsamp',z_sample)
    plt.close()
    plt.plot(r1,z1)
    plt.scatter(r_sample,z_sample,alpha=0.1)
    plt.title('Outer divertor')
    plt.ylabel('z [m]')
    plt.xlabel('r [m]')
    plt.savefig('particleScatter.png')

    plt.close()
    plt.hist(r_sample)
    plt.savefig('hist.png')

    rPar = np.divide((r2 - r1),length)
    zPar = np.divide((z2 - z1),length)
    #print('should be ones', np.sqrt(np.multiply(rPar,rPar) + np.multiply(zPar,zPar)))
    parVec = np.zeros([len(gitr_inds),3])
    parVec[:,0] = rPar
    parVec[:,2] = zPar

    print('inDir', inDir)
    rPerp = 0*slope;
    zPerp = 0*slope;
    for i in range(0,len(slope)):
         if slope[i]==0:
             perpSlope = 1.0e12
         else:
             perpSlope = -np.sign(slope[i])/np.abs(slope[i]);

         rPerp[i] = -inDir[i]/np.sqrt(perpSlope*perpSlope+1);
         zPerp[i] = -inDir[i]*np.sign(perpSlope)*np.sqrt(1-rPerp[i]*rPerp[i]);
    #perpSlope = -np.sign(slope)/np.abs(slope);
    #rPerp = -inDir/np.sqrt(perpSlope*perpSlope+1);
    #zPerp = -inDir*np.sign(perpSlope)*np.sqrt(1-rPerp*rPerp);
    perpVec = np.zeros([len(gitr_inds),3])
    perpVec[:,0] = rPerp
    perpVec[:,2] = zPerp
    print('perpVec',perpVec)

    #moves particles 10um away from surface
    surface_buffer = 1.0e-5
    r_sample = r_sample + surface_buffer*rPerp[mins]
    z_sample = z_sample + surface_buffer*zPerp[mins]

    particle_energy = 4*np.ones(nP)
    particle_angle = np.zeros(nP)

    randDistvTheta = np.random.rand(nP)*2*np.pi
    vr = np.sqrt(2*particle_energy*1.602e-19/184/1.66e-27);
    #vx = np.multiply(vr,np.multiply(np.cos(randDistvTheta),np.sin(np.pi/180.0*particle_angle)));
    #vy = np.multiply(vr,np.multiply(np.sin(randDistvTheta),np.sin(np.pi/180.0*particle_angle)));
    #vz = np.multiply(vr,np.cos(np.pi/180.0*particle_angle));
    vx = 0*vr;
    vy = 0*vr;
    vz = vr;
    parVec2 = np.cross(parVec,perpVec)
    vxP = 0*vr;
    vyP = 0*vr;
    vzP = 0*vr;

    for i in range(nP):
        vxP[i] = parVec[mins[i],0]*vx[i] + parVec2[mins[i],0]*vy[i] + perpVec[mins[i],0]*vz[i];
        vyP[i] = parVec[mins[i],1]*vx[i] + parVec2[mins[i],1]*vy[i] + perpVec[mins[i],1]*vz[i];
        vzP[i] = parVec[mins[i],2]*vx[i] + parVec2[mins[i],2]*vy[i] + perpVec[mins[i],2]*vz[i];
        #if z_sample[i] < -3.8 :
            #print('r,z',r_sample[i],z_sample[i])
            #print('perpVec xyz',perpVec[mins[i],0],perpVec[mins[i],1],perpVec[mins[i],2])
            #print('v xyz',vxP[i],vyP[i],vzP[i])

    rootgrp = netCDF4.Dataset("particleSource.nc", "w", format="NETCDF4")
    npp = rootgrp.createDimension("nP", nP)
    xxx = rootgrp.createVariable("x","f8",("nP"))
    yyy = rootgrp.createVariable("y","f8",("nP"))
    zzz = rootgrp.createVariable("z","f8",("nP"))
    vxx = rootgrp.createVariable("vx","f8",("nP"))
    vyy = rootgrp.createVariable("vy","f8",("nP"))
    vzz = rootgrp.createVariable("vz","f8",("nP"))
    xxx[:] = r_sample
    yyy[:] = 0*r_sample
    zzz[:] = z_sample
    vxx[:] = vxP
    vyy[:] = vyP
    vzz[:] = vzP
    rootgrp.close()
    
if __name__ == "__main__":
    simple2D()