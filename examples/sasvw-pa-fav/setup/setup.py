import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import shutil
import numpy as np
import solpsProcessing, makeGeom, makeParticleSource

nP = int(1e3)
W_surface_indices = np.arange(10,24)

makeGeom.main(gitr_geometry_filename='gitrGeometry.cfg', \
                    solps_geomfile = 'assets/sas-vw_v004.ogr', \
                    solps_targfile = 'assets/b2fgmtry', \
                    profiles_file = '../input/plasmaProfiles.nc', \
                    surfW = W_surface_indices, \
                    solps_rz = 'assets/solps_rz.txt', \
                    gitr_rz = 'assets/gitr_rz.txt', \
                    rmrs_fine_file = 'assets/rmrs_fine.txt', \
                    W_fine_file = 'assets/W_fine.txt', \
                    numAddedPoints = 100, \
                    plot_variables = 0)

os.remove('gitrGeometry.cfg0')
shutil.move('gitrGeometry.cfg', '../input/gitrGeometry.cfg')

solpsProcessing.readEquilibrium(equilibrium_filename = 'assets/dg.equ', \
                    W_indices = W_surface_indices, \
                    solps_geom = 'assets/b2fgmtry', \
                    plot_variables = 0)

shutil.move('bField.nc', '../input/bField.nc')
'''
solpsProcessing.plot_surf_plasma_params(W_surf = W_surface_indices)

makeParticleSource.point_source(nP)
'''
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
