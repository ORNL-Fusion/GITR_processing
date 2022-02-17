import sys
sys.path.insert(0, '../../../python/')

import shutil
import gitr
import solps
import gitrParticleSource

make_gitr_geometry_from_solps_sasvw(gitr_geometry_filename = 'gitrGeometry.cfg', \
                                   solps_geomfile = 'assets/vvfile')

shutil.move('gitrGeometry.cfg', '../input/gitrGeometry.cfg')

solps.process_solps_output_for_gitr(dakota_filename = 'assets/dakota', \
                                   nR = 500, nZ = 1000, plot_variables=0, \
                                   b2fstate_filename = 'assets/b2fstate')

shutil.move('profiles.nc', '../input/profiles.nc')

solps.readEquilibrium(filename = 'assets/west_54034_10p2s_mag.X4.equ', \
                                    solps_mesh_extra = 'assets/mesh.extra', \
                                    solps_geom = 'assets/b2fgmtry')

shutil.move('bField.nc', '../input/bField.nc')

solps.make_solps_targ_coord_file(gitr_geom_filename = 'gitrGeometry.cfg', \
                                    solps_geom = 'assets/b2fgmtry', \
                                    coords_file = 'assets/right_target_coordinates.txt', \
                                    right_target_filename = 'assets/rightTargOutput')

solps.make_solps_targ_file(solps_geom = 'assets/b2fgmtry', \
                                    b_field_file = 'assets/west_54034_10p2s_mag.X4.equ', \
                                    coords_file = 'assets/right_target_coordinates.txt', \
                                    right_target_filename= 'assets/rightTargOutput')

gitrParticleSource.particleSource2d_west(nParticles = int(1e3), \
                                    geom = 'gitrGeometry.cfg', \
                                    targFile = 'assets/rightTargOutput', \
                                    coordsFile = 'assets/right_target_coordinates.txt')

shutil.move('particleSource.nc', '../input/particleSource.nc')
