import os
import subprocess
import shutil
import numpy as np
import libconf
import io
import impacts

filename = "input/gitrInput.cfg"
angles = np.empty(10)
a=0
for i in range(len(angles)):
    angles[i] = 0+a
    a = a + 360/(len(angles)-1)
Bmag = 2.2970
Btheta = 90.5289;
v_parallel = 24561; #24561;
bg = ["D"] #,"He"];
ZZ = [1] #,2];
amu = [2] #,4];
for k in range(0,len(bg)): #plasma species
    for i in range(0,len(angles)):
    
        with io.open(filename) as f:
            config = libconf.load(f)
            #config['backgroundPlasmaProfiles']['Z'] = ZZ[k]
            #config['backgroundPlasmaProfiles']['amu'] = amu[k]
            #config['backgroundPlasmaProfiles']['Bfield']['r'] = Bmag*np.cos(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
            #config['backgroundPlasmaProfiles']['Bfield']['z'] = Bmag*np.cos(np.deg2rad(90 - angles[i]))
            #config['backgroundPlasmaProfiles']['Bfield']['y'] = Bmag*np.sin(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
            config['backgroundPlasmaProfiles']['FlowVelocity']['flowVr'] = -v_parallel*np.cos(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
            config['backgroundPlasmaProfiles']['FlowVelocity']['flowVz'] = -v_parallel*np.cos(np.deg2rad(90 - angles[i]))
            config['backgroundPlasmaProfiles']['FlowVelocity']['flowVy'] = -v_parallel*np.sin(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
    
        with io.open(filename, 'w') as f:
            libconf.dump(config, f)
        
        impacts.setup_case(Zi = 6, amu_i = 12, charge_i = 4,
               Zb = ZZ[k], amu_b = amu[k], charge_b = 1,
               Ti = 19.334,
               ni = 4.68e19, Te = 18, ne = 4.68e19, Bmag = 2.297, phi = 90+angles[i],
               nP = 1e5, nT = 1e4, dt = 1.0e-9, height_factor = 4,
               filename = "input/gitrInput.cfg")
        subprocess.run("../build/GITR", shell=True, check=True)
        #shutil.copyfile("output/surface.nc","surface_C"+str(j)+"_loc_"+str(i) +".nc")
        shutil.copyfile("output/positions.nc","positions_"+str(i)+"_" +bg[k] + ".nc")
        os.remove("output/surface.nc")
        os.remove("output/positions.nc")
        os.remove("output/positions.m")
        os.remove("output/particleSource.nc")
