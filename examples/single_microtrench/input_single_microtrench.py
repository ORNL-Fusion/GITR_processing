# -*- coding: utf-8 -*-
import electronvolt_num as units
import pyGITR
import numpy as np
from single_microtrench  import Lxbc,Lybc
ParticleFile='particleConf.nc'
GeometryFile='gitrGeom.cfg'
B0 = 2.25
thetaB = 2
phiB = 0
AxisrotB = [0,1,0]
nP=20000
mi=2.0
Elem='D'

MassElem ={'C' : 12, 'D' : 2}
Mimp = MassElem.get(Elem)
charge = 1
Zmax = 6
Te = 25
Ti = Te
V= -3*Te
#Plasma parameter

Timp = Te
vthi = np.sqrt(units.eV*Ti/units.mp/mi)
vimp = np.sqrt(units.eV*Ti/units.mp/Mimp)
cs = np.sqrt(units.eV*2*Ti/units.mp/mi)
RootPath = '~/simulations/GITR/data/single_microtrench_{}_{}_{}__{}/'.format(Elem,charge, Ti,nP)


##

p = pyGITR.ParticleDistribution()
p.SetAttr('Np', nP+1)

# Set positions of particles
p.SetAttr('x','Uniform',xmin=-Lxbc*0.99,xmax=Lxbc*0.99) #set values of y and z with uniformly distributed values between -0.05 and 0.05
p.SetAttr('y','Uniform',xmin=-Lybc*0.99,xmax=Lybc*0.99)
p.SetAttr('z',0.01) # set all values of x to -0.01
# p.SetAttr(['x','y'],0) #set values of y and z with uniformly distributed values between -0.05 and 0.05

# Set velocities of particles
p.SetAttr(['vx'],'Gaussian',sigma=1, beta=0.5)
p.SetAttr(['vy','vz'],'Gaussian',sigma=1)
p.Particles['vx'].mean()

# Rescale velocity by characteristic velocity
vpara = vimp
vperp = vimp
p.ScaleAttr(['vy','vz'],vperp)
p.ScaleAttr('vx',vpara)
p.Rotate('v', AxisrotB, thetaB) #rotate (vx,vy,vz) along AxisrotB of angle thetaB
vx = p.Particles['vx']
vy = p.Particles['vy']
vz = p.Particles['vz']
E= 1/2*Mimp*(vx**2+vy**2+vz**2)*units.mp/units.eV
#%%
# Write particle distribution in netcdf file
p.WriteParticleFile(ParticleFile)

Input = pyGITR.Input()
Input.SetBField(B0=B0, theta = thetaB, phi = phiB)
Input.SetTimeStep(dt=1e-10, nT=200000)
Input.SetGeometryFile(GeometryFile)
Input.SetParticleSource(ParticleFile, nP=nP, Zmax=Zmax, M=Mimp, Z=charge)
Input.SetBackgroundPlasmaProfiles(Voltage=V)
Input.SetSurfaces()
Input.SetDiagnostics()
Input.Input['flags'] = {
                    'USE_CUDA':0,
                    'USE_MPI':1,
                    'USE_OPENMP':1,
                    'USE_IONIZATION':1,
                    'USE_RECOMBINATION':1,
                    'USEPERPDIFFUSION':0,
                    'USEPARDIFFUSION':0,
                    'USECOULOMBCOLLISIONS':0,
                    'USEFRICTION':0,
                    'USEANGLESCATTERING':0,
                    'USEHEATING':0,
                    'USETHERMALFORCE':0,
                    'USESURFACEMODEL':0,
                    'USE_SURFACE_POTENTIAL':0,
                    'USE_SHEATHEFIELD':1,
                    'BIASED_SURFACE':3,
                    'USEPRESHEATHEFIELD':0,
                    'BFIELD_INTERP':0,
                    'LC_INTERP':0,
                    'GENERATE_LC':0,
                    'EFIELD_INTERP':0,
                    'PRESHEATH_INTERP':0,
                    'DENSITY_INTERP':0,
                    'TEMP_INTERP':0,
                    'FLOWV_INTERP':0,
                    'GRADT_INTERP':0,
                    'ODEINT':0,
                    'FIXED_SEEDS':1,
                    'PARTICLESEEDS':0,
                    'GEOM_TRACE':0,
                    'GEOM_HASH':0,
                    'GEOM_HASH_SHEATH':0,
                    'PARTICLE_TRACKS':0,
                    'PARTICLE_SOURCE_SPACE':0,
                    'PARTICLE_SOURCE_ENERGY':0,
                    'PARTICLE_SOURCE_ANGLE':0,
                    'PARTICLE_SOURCE_FILE':1,
                    'SPECTROSCOPY':0,
                    'USE3DTETGEOM':1,
                    'USECYLSYMM':0,
                    'USEFIELDALIGNEDVALUES':0,
                    'FLUX_EA':1,
                    'FORCE_EVAL':0,
                    'USE_SORT':0,
                    'CHECK_COMPATIBILITY':1
                    }

Input.WriteInputFile(Folder='input', OverWrite=True)
#%%
Run = pyGITR.Run()
Run.GITRExecPath = '/home/guterlj/boundary/GITR/build/GITR'
Run.Verbose = True

Run.SetReferenceDirectory('.')
Run.SetSimRootPath(RootPath)
Run.AddParamScan('input/gitrGeom.cfg',{'geom.lambda':[0.00005*1.2,0.0001*1.2]})
Run.AddParamScan('input/gitrGeom.cfg',{'geom.method':[k for k in range(100)]})
#Run.AddParamScan('input/gitrGeom.cfg',{'geom.lambda':[ 0.0001, 0.005]})
#Run.AddParamScan('input/gitrInput.cfg',{'backgroundPlasmaProfiles.biasPotential':[-3*Ti, 0]})
Run.ModifParam('input/gitrInput.cfg','impurityParticleSource.nP',nP)

Run.SetupScan(OverWrite=True)


Run.Clean()
Run.LaunchSbatch(p='preemptable',n=1)
#Run.LaunchBatch()

#%%
