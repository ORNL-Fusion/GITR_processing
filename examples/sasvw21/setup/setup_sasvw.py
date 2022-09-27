import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import shutil
import solps
import make_geom_sasv6
import make_ParticleSource_sasvw

r,z, rW,zW = make_geom_sasv6.V6e_v002(gitr_geometry_filename='gitrGeometry.cfg', \
                                    solps_geomfile = 'assets/geom-SASV6/SAS-V6e_v002.ogr', \
                                    solps_targfile = 'assets/b2fgmtry', \
                                    solps_rz = 'assets/geom-SASV6/solps_rz.txt', \
                                    gitr_rz = 'assets/geom-SASV6/gitr_rz.txt')

os.remove('gitrGeometry.cfg0')
shutil.move('gitrGeometry.cfg', '../input/gitrGeometry.cfg')
'''
solps.readEquilibrium(filename = 'assets/vertex_sasvw.eq', \
                                    solps_geom = 'assets/b2fgmtry', \
                                    solps_mesh_extra = None, \
                                    plot_variables = 0, \
                                    r_wall = r, z_wall = z)

shutil.move('bField.nc', '../input/bField.nc')

solps.process_solps_output_for_gitr(dakota_filename = 'assets/dakota', \
                                   nR = 500, nZ = 1000, plot_variables = 0, \
                                   b2fstate_filename = 'assets/b2fstate', \
                                   r_wall = r, z_wall = z)

shutil.move('profiles.nc', '../input/profiles.nc')

solps.make_solps_targ_coord_file(gitr_geom_filename = '../input/gitrGeometry.cfg', \
                                    solps_geom = 'assets/b2fgmtry', \
                                    coords_file = 'assets/right_target_coordinates.txt', \
                                    right_target_filename = 'assets/rightTargOutput')
'''
make_ParticleSource_sasvw.point_source(nP = int(2e2))
'''
make_ParticleSource_sasvw.simple2D(nP = int(3e2), \
                                    geom = '../input/gitrGeometry.cfg', \
                                    targFile = 'assets/rightTargOutput', \
                                    coordsFile = 'assets/right_target_coordinates.txt', \
                                    profilesFile = '../input/profiles.nc', \
                                    ftBFile = 'assets/ftridynBackground.nc', \
                                    configuration = 'midpoint', \
                                    plot_variables = 1, \
                                    r_W = rW, z_W = zW)
'''
shutil.move('particleSource.nc', '../input/particleSource.nc')
