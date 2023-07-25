import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import shutil
import numpy as np
import solpsProcessing, makeGeom, makeParticleSource

nP = int(1e4)
W_indices = np.arange(11,22)

makeGeom.main(gitr_geometry_filename='gitrGeometry.cfg', \
                    solps_geomfile = 'assets/sas-vw_v005_mod.ogr', \
                    solps_targfile = 'assets/b2fgmtry', \
                    profiles_file = '../input/plasmaProfiles.nc', \
                    W_indices_profiles = W_indices, \
                    numAddedPoints = 100, \
                    plot_variables = 0)

os.remove('gitrGeometry.cfg0')
shutil.move('gitrGeometry.cfg', '../input/gitrGeometry.cfg')

solpsProcessing.readEquilibrium(equilibrium_filename = 'assets/dg.equ', \
                    W_indices = W_indices, \
                    solps_geom = 'assets/b2fgmtry', \
                    plot_variables = 1)

shutil.move('bField.nc', '../input/bField.nc')

solpsProcessing.plot_surf_plasma_params(W_surf = W_indices)
'''
makeParticleSource.point_source(nP)
'''
makeParticleSource.distributed_source(nP, surfW = W_indices, \
                    geom = '../input/gitrGeometry.cfg', \
                    profiles_file = '../input/plasmaProfiles.nc', \
                    ftDFile = 'assets/ftridynBackgroundD.nc', \
                    ftCFile = 'assets/ftridynBackgroundC.nc', \
                    configuration = 'random', \
                    plot_variables = 1)

shutil.move('particleSource.nc', '../input/particleSource.nc')
