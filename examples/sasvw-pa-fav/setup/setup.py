import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import shutil
import solpsProcessing
import makeGeom
import makeParticleSource

solpsProcessing.readEquilibrium(equilibrium_filename = 'assets/dg.equ', \
                                    solps_geom = 'assets/b2fgmtry', \
                                    solps_mesh_extra = None, \
                                    plot_variables = 1)

shutil.move('bField.nc', '../input/bField.nc')

solpsProcessing.plot_surf_plasma_params()

'''
r,z, rW,zW, rCoarse, zCoarse, addedPoints = makeGeom.V6e_v002(gitr_geometry_filename='gitrGeometry.cfg', \
                                    solps_geomfile = 'assets/geom-SASV6/SAS-V6e_v002.ogr', \
                                    solps_targfile = 'assets/b2fgmtry', \
                                    solps_rz = 'assets/solps_rz.txt', \
                                    gitr_rz = 'assets/gitr_rz.txt', \
                                    surf_coarse = 'assets/surf_coarse.txt', \
                                    surf_ind = 'assets/surf_ind.txt', \
                                    numAddedPoints = 100)

os.remove('gitrGeometry.cfg0')
shutil.move('gitrGeometry.cfg', '../input/gitrGeometry.cfg')


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
