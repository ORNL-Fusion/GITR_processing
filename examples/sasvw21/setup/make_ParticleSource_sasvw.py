import sys
sys.path.insert(0, '../../../python/')
sys.path.insert(0, '../../../../pyGITR/pyGITR')

import numpy as np
import matplotlib.pyplot as plt

import netCDF4
import gitr
import solps
import Particles

def simple(nP = int(1e3), \
            geom = 'gitrGeometry.cfg', \
            targFile = 'assets/rightTargOutput', \
            coordsFile = 'assets/right_target_coordinates.txt'):
    
    PartDist = Particles.ParticleDistribution(nP, ListAttr=['x','y','z','vx','vy','vz'])
    
    #get r,z,t distributions for sputtered W
    r = PartDist.Generate(nP, 'Uniform')
    z = PartDist.Generate(nP, 'Uniform')
    t = PartDist.Generate(nP, 'Uniform')
    
    #convert r,z,t to x,y,z
    
    #get IEADs for sputtered W
    E = PartDist.Generate(nP, 'Thomson')
    PolAng = PartDist.Generate(nP, 'SinCos')
    AziAng = PartDist.Generate(nP, 'Uniform')
    
    #convert IEADs to vx,vy,vz
    
    
def old(nP = int(1e3), \
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
    simple()