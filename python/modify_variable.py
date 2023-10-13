import shutil
import netCDF4

filepath = 'input/plasmaProfiles0.nc'
filepath_new = 'input/plasmaProfiles.nc'

shutil.copyfile(filepath, filepath_new)

dset = netCDF4.Dataset(filepath_new,'r+')
vr = dset['vr'][:]
vt = dset['vt'][:]
vz = dset['vz'][:]
dset['vr'][:] = 1.0*vr;
dset['vt'][:] = -1.0*vt;
dset['vz'][:] = 1.0*vz;
#vz[vz < 0] = -flow_v[k]
dset.close()
