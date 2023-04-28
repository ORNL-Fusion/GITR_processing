import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import shutil
import numpy as np
import solpsProcessing, makeGeom#, makeParticleSource

nP = int(5e2)

makeGeom.main(gitr_geometry_filename='gitr_geometry.cfg', \
            solps_rz = 'assets/solps_rz.txt', \
            gitr_rz = 'assets/gitr_rz.txt', \
            solps_mesh_extra='assets/mesh.extra', \
            profiles_filename = '../input/plasmaProfiles.nc', \
            rmrs_fine_file = 'assets/rmrs_fine.txt', \
            W_fine_file = 'assets/W_fine.txt', \
            numAddedPoints = 20, \
            solps_geom = 'assets/b2fgmtry')


with open('assets/W_fine.txt') as f:
    Wdata = f.readlines()
W_surface_indices = np.empty(len(Wdata))
for i in range(len(Wdata)):
    W_surface_indices[i] = int(Wdata[i])
    
os.remove('gitr_geometry.cfg0')
shutil.move('gitr_geometry.cfg', '../input/gitr_geometry.cfg')

solpsProcessing.readEquilibrium(equilibrium_filename = 'assets/west_54034_10p2s_mag.X4.equ', \
                    W_indices = W_surface_indices, \
                    solps_geom = 'assets/b2fgmtry', \
                    solps_mesh_extra = 'assets/mesh.extra', \
                    plot_variables = 0)


shutil.move('bField.nc', '../input/bField.nc')
'''
solpsProcessing.plot_surf_plasma_params(W_surf = W_surface_indices)

makeParticleSource.point_source(nP)

makeParticleSource.distributed_source(nP, surfW = W_surface_indices, \
                    geom = '../input/gitrGeometry.cfg', \
                    profiles_file = '../input/plasmaProfiles.nc', \
                    gitr_rz = 'assets/gitr_rz.txt', \
                    rmrs_fine_file = 'assets/rmrs_fine.txt', \
                    W_fine_file = 'assets/W_fine.txt', \
                    ftDFile = 'assets/ftridynBackgroundD.nc', \
                    ftCFile = 'assets/ftridynBackgroundC.nc', \
                    configuration = 'random', \
                    plot_variables = 0)

shutil.move('particleSource.nc', '../input/particleSource.nc')
'''