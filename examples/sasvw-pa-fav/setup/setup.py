import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import shutil
import numpy as np
import solpsProcessing
import makeGeom
import makeParticleSource
'''
solpsProcessing.readEquilibrium(equilibrium_filename = 'assets/dg.equ', \
                                    solps_geom = 'assets/b2fgmtry', \
                                    solps_mesh_extra = None, \
                                    plot_variables = 1)

shutil.move('bField.nc', '../input/bField.nc')

#solpsProcessing.plot_surf_plasma_params()

'''
makeGeom.main(gitr_geometry_filename='gitrGeometry.cfg', \
                                    solps_geomfile = 'assets/sas-vw_v004.ogr', \
                                    solps_targfile = 'assets/b2fgmtry', \
                                    profiles_file = '../input/plasmaProfiles.nc', \
                                    surfW = np.arange(10,24), \
                                    solps_rz = 'assets/solps_rz.txt', \
                                    gitr_rz = 'assets/gitr_rz.txt', \
                                    rmrs_fine_file = 'assets/rmrs_fine.txt', \
                                    W_fine_file = 'assets/W_fine.txt', \
                                    numAddedPoints = 100)

os.remove('gitrGeometry.cfg0')
shutil.move('gitrGeometry.cfg', '../input/gitrGeometry.cfg')
'''

makeParticleSource.point_source(nP = int(5e2))

make_ParticleSource_sasvw.simple2D(nP = int(1e3), \
                                    geom = '../input/gitrGeometry.cfg', \
                                    targFile = 'assets/rightTargOutput', \
                                    wallFile = 'assets/gitr_rz.txt', \
                                    surf_coarse = 'assets/surf_coarse.txt', \
                                    surf_ind = 'assets/surf_ind.txt', \
                                    profilesFile = '../input/profiles.nc', \
                                    ftDFile = 'assets/ftridynBackgroundD.nc', \
                                    ftCFile = 'assets/ftridynBackgroundC.nc', \
                                    configuration = 'random', \
                                    plot_variables = 0, \
                                    r_W=rW, z_W=zW, rCoarse=rCoarse, zCoarse=zCoarse, addedPoints=addedPoints)

shutil.move('particleSource.nc', '../input/particleSource.nc')
'''
