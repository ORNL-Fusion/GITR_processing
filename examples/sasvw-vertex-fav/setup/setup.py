import sys, os
sys.path.insert(0, os.path.abspath('../../../python/'))

import shutil
import numpy as np
import solpsProcessing, makeGeom, makeParticleSource

nP = int(1e6)
run_directory = '..'

W_indices = np.arange(16,25)
tile_shift_indices = [2,6]
Bangle_shift_indices = [3,6]

print_separator = '\n-------------------------------------------------\n'

print('\n',print_separator,'Making gitrGeometry.cfg',print_separator,'\n')
makeGeom.main(gitr_geometry_filename='gitrGeometry.cfg', \
                    solps_geomfile = 'assets/sas-vw_v005_mod.ogr', \
                    solps_targfile = 'assets/b2fgmtry', \
                    profiles_file = run_directory+'/input/plasmaProfiles.nc', \
                    W_indices_profiles = W_indices, \
                    tile_shift_indices = tile_shift_indices, \
                    numAddedPoints = 100, \
                    plot_variables = 0, \
                    show_plots = 0)

os.remove('gitrGeometry.cfg0')
shutil.move('gitrGeometry.cfg', run_directory+'/input/gitrGeometry.cfg')

print('\n',print_separator,'Making bField.nc',print_separator,'\n')
solpsProcessing.readEquilibrium(equilibrium_filename = 'assets/dg.equ', \
                    W_indices = W_indices, \
                    solps_geom = 'assets/b2fgmtry', \
                    flip_Bt = True, \
                    plot_variables = 0)

shutil.move('bField.nc', run_directory+'/input/bField.nc')
'''
solpsProcessing.plot_surf_plasma_params(W_surf = W_indices, \
                    tile_shift_indices = tile_shift_indices, \
                    Bangle_shift_indices = Bangle_shift_indices)

makeParticleSource.point_source(nP)
'''
print('\n',print_separator,'Making particleSource.nc',print_separator,'\n')
makeParticleSource.distributed_source(nP, surfW = W_indices, \
                    tile_shift_indices = tile_shift_indices, \
                    Bangle_shift_indices = Bangle_shift_indices, \
                    geom = run_directory+'/input/gitrGeometry.cfg', \
                    profiles_file = run_directory+'/input/plasmaProfiles.nc', \
                    ftDFile = 'assets/ftridynBackgroundD.nc', \
                    ftCFile = 'assets/ftridynBackgroundC.nc', \
                    configuration = 'random', \
                    plot_variables = 0)

shutil.move('particleSource.nc', run_directory+'/input/particleSource.nc')
