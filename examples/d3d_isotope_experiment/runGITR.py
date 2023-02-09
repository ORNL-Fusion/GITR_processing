import os
import subprocess
import shutil
import numpy as np
import libconf
import io

filename = "input/gitrInput.cfg"
angles = [1.2525, 3.0, 6.0, 9.0, 15.0, 30.0, 60.0, 90.0]
Bmag = 5; #0.7; #2.25; #2.2970
Btheta = 90; #90.5289;
v_parallel = 24561; #24561;
bg = ["D","He"];
ZZ = [1,2];
amu = [2,4];
siz = ["0p5mm","1mm"];
for k in range(0,len(bg)): #plasma species
    for j in range(0,len(siz)): #spot size
        for i in range(0,len(angles)):
        
            with io.open(filename) as f:
                config = libconf.load(f)
                config['particleSource']['ncFileString'] = "particle_source_"+siz[j] + ".nc"
                config['backgroundPlasmaProfiles']['Z'] = ZZ[k]
                config['backgroundPlasmaProfiles']['amu'] = amu[k]
                config['backgroundPlasmaProfiles']['Bfield']['r'] = Bmag*np.cos(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
                config['backgroundPlasmaProfiles']['Bfield']['z'] = Bmag*np.cos(np.deg2rad(90 - angles[i]))
                config['backgroundPlasmaProfiles']['Bfield']['y'] = Bmag*np.sin(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
                config['backgroundPlasmaProfiles']['FlowVelocity']['flowVr'] = -v_parallel*np.cos(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
                config['backgroundPlasmaProfiles']['FlowVelocity']['flowVz'] = -v_parallel*np.cos(np.deg2rad(90 - angles[i]))
                config['backgroundPlasmaProfiles']['FlowVelocity']['flowVy'] = -v_parallel*np.sin(np.deg2rad(Btheta))*np.sin(np.deg2rad(90 - angles[i]))
        
            with io.open(filename, 'w') as f:
                libconf.dump(config, f)
            subprocess.run("/home/tqd/Code/GITR/build/GITR", shell=True, check=True)
            #shutil.copyfile("output/surface.nc","surface_C"+str(j)+"_loc_"+str(i) +".nc")
            shutil.copyfile("output/positions.nc","positions_"+str(i)+"_" +bg[k]+"_"+siz[j]+ ".nc")
            os.remove("output/surface.nc")
            os.remove("output/positions.nc")
            os.remove("output/positions.m")
            os.remove("output/particleSource.nc")